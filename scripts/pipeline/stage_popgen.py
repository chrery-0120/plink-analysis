"""Stage 2: population structure, relatedness, LD, ROH and FST.

Everything here runs on the QC-cleaned fileset (``01_qc/data/clean``) and the
LD-pruned variant list (``01_qc/data/clean_pruned.prune.in``).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import utils

try:
    from . import plots
except Exception:  # noqa: BLE001
    plots = None  # type: ignore[assignment]

LD_BIN_EDGES_KB = [0, 1, 5, 10, 25, 50, 100, 250, 500, 1000, 1e9]


def _try_plot(func_name: str, *args, **kwargs) -> str | None:
    if plots is None:
        return None
    func = getattr(plots, func_name, None)
    if func is None:
        return None
    try:
        return str(func(*args, **kwargs))
    except Exception as exc:  # noqa: BLE001
        return f"__error__:{exc}"


def run(ctx) -> dict[str, Any]:
    cfg = ctx.config
    stage = ctx.stage_dir
    tables, plots_dir, data_dir = stage / "tables", stage / "plots", stage / "data"
    run_dir = ctx.run_dir
    chr_flags = list(cfg.get("_derived", {}).get("chr_flags") or [])
    clean = run_dir / "01_qc" / "data" / "clean"
    prune_in = run_dir / "01_qc" / "data" / "clean_pruned.prune.in"
    plink: utils.Plink = ctx.plink
    figures: list[str] = []
    plot_errors: list[str] = []
    warnings: list[str] = []
    tables_map: dict[str, str] = {}

    if not clean.with_suffix(".bed").exists():
        raise FileNotFoundError(f"缺少 QC 清理后的数据 {clean}.bed，请先运行 qc 阶段")
    if not prune_in.exists():
        raise FileNotFoundError(f"缺少 LD 剪枝列表 {prune_in}，请先运行 qc 阶段")

    fam = utils.read_fam(clean.with_suffix(".fam"))
    bim = utils.read_bim(clean.with_suffix(".bim"))
    n_samples = len(fam)

    # ------------------------------------------------------------------
    # PCA
    # ------------------------------------------------------------------
    n_pcs = int(cfg["popgen"].get("pca", 20))
    pca_prefix = tables / "pca"
    plink.run(["--bfile", str(clean), *chr_flags, "--autosome",
               "--extract", str(prune_in), "--pca", str(n_pcs)],
              stage="popgen", name="pca", out_prefix=pca_prefix,
              expect=[".eigenvec", ".eigenval"])
    eigenvec = utils.read_eigenvec(pca_prefix.with_suffix(".eigenvec"))
    eigenval = utils.read_eigenval(pca_prefix.with_suffix(".eigenval"))
    # PLINK 1.9 eigenvalues come from the variance-standardised GRM: they sum to
    # the sample count, so variance explained = eigenvalue / N.
    var_explained = eigenval / float(len(eigenvec)) if len(eigenvec) else np.zeros_like(eigenval)
    eigenvec.to_csv(tables / "pca.eigenvec", index=False, encoding="utf-8-sig")
    pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(len(eigenval))],
                  "eigenvalue": eigenval, "variance_explained": var_explained}
                 ).to_csv(tables / "pca.eigenval.csv", index=False, encoding="utf-8-sig")
    tables_map["pca_eigenvec"] = "tables/pca.eigenvec"
    tables_map["pca_eigenval"] = "tables/pca.eigenval.csv"

    groups: pd.Series | None = None
    group_present = bool(cfg.get("_derived", {}).get("group_present"))
    if group_present:
        group_path = run_dir / "00_intake" / "data" / "group_clean.txt"
        if group_path.exists():
            group_df = pd.read_csv(group_path, sep=r"\s+", header=None,
                                   names=["FID", "IID", "GROUP"], dtype=str)
            groups = group_df.set_index("IID")["GROUP"]

    # PCA outlier flag (SD rule on first two PCs)
    pca_outliers: list[str] = []
    if len(eigenvec) > 3:
        mask = pd.Series(False, index=eigenvec.index)
        for col in [c for c in eigenvec.columns if c.startswith("PC")][:2]:
            values = pd.to_numeric(eigenvec[col], errors="coerce")
            sd = float(values.std(ddof=1))
            if np.isfinite(sd) and sd > 0:
                mask |= (values - values.mean()).abs() > float(cfg["qc"]["pca_outlier_sd"]) * sd
        pca_outliers = eigenvec.loc[mask, "IID"].astype(str).tolist()

    def add_fig(name: str, path: str | None) -> None:
        if path is None:
            return
        if path.startswith("__error__:"):
            plot_errors.append(f"{name}: {path.split(':', 1)[1]}")
        else:
            figures.append(path)

    add_fig("pca_scatter", _try_plot("plot_pca_scatter", eigenvec, list(eigenval),
                                     plots_dir / "pca.png", groups=groups,
                                     n_samples=len(eigenvec)))
    add_fig("pca_scree", _try_plot("plot_pca_scree", list(eigenval),
                                   plots_dir / "pca_scree.png", n_samples=len(eigenvec)))

    # ------------------------------------------------------------------
    # kinship (GCTA-style relationship matrix) and IBD
    # ------------------------------------------------------------------
    kinship_prefix = tables / "kinship"
    plink.run(["--bfile", str(clean), *chr_flags, "--autosome",
               "--extract", str(prune_in), "--make-rel", "square"],
              stage="popgen", name="kinship", out_prefix=kinship_prefix,
              expect=[".rel"])
    kinship_summary: dict[str, Any] = {"available": False}
    relationship = None
    rel_ids: list[str] = []
    try:
        relationship = np.loadtxt(kinship_prefix.with_suffix(".rel"), dtype=float)
        rel_id_file = kinship_prefix.with_suffix(".rel.id")
        if rel_id_file.exists():
            rel_df = pd.read_csv(rel_id_file, sep=r"\s+", header=None, dtype=str)
            rel_ids = rel_df[1].astype(str).tolist()
        off_diag = relationship[~np.eye(relationship.shape[0], dtype=bool)]
        kinship_summary = {
            "available": True,
            "n": int(relationship.shape[0]),
            "mean_offdiag": round(float(np.nanmean(off_diag)), 5),
            "max_offdiag": round(float(np.nanmax(off_diag)), 5),
            "n_pairs_gt_025": int((off_diag > 0.25).sum() // 2),
            "n_pairs_gt_05": int((off_diag > 0.5).sum() // 2),
        }
        if len(rel_ids) == relationship.shape[0] and len(rel_ids) <= 300:
            add_fig("kinship_heatmap", _try_plot("plot_kinship_heatmap", relationship,
                                                 rel_ids, plots_dir / "kinship.png"))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"亲缘矩阵解析失败: {exc}")
        plot_errors.append(f"kinship: {exc}")

    genome_prefix = tables / "genome"
    n_pairs = n_samples * (n_samples - 1) // 2
    genome_args = ["--bfile", str(clean), *chr_flags, "--autosome",
                   "--extract", str(prune_in), "--genome"]
    if n_samples > 2000:
        genome_args += ["--min", "0.1"]
    plink.run(genome_args, stage="popgen", name="genome", out_prefix=genome_prefix,
              expect=[".genome"])
    genome = utils.read_csv(genome_prefix.with_suffix(".genome"), sep=r"\s+")
    if "PI_HAT" in genome.columns:
        genome["PI_HAT"] = pd.to_numeric(genome["PI_HAT"], errors="coerce")
        # keep the original PLINK .genome file untouched (tables/genome.genome);
        # write the curated pair table separately so downstream tools can read both
        genome.to_csv(tables / "genome_pairs.csv", sep="\t", index=False, encoding="utf-8-sig")
        tables_map["genome"] = "tables/genome.genome"
        tables_map["genome_pairs"] = "tables/genome_pairs.csv"
        add_fig("relatedness", _try_plot("plot_relatedness_hist",
                                         genome["PI_HAT"], plots_dir / "relatedness.png"))
        relatedness_summary = {
            "n_pairs": int(len(genome)),
            "mean_pihat": round(float(genome["PI_HAT"].mean()), 5),
            "max_pihat": round(float(genome["PI_HAT"].max()), 5),
            "n_pairs_gt_025": int((genome["PI_HAT"] >= 0.25).sum()),
            "n_pairs_gt_09": int((genome["PI_HAT"] >= 0.90).sum()),
        }
    else:
        relatedness_summary = {"n_pairs": 0}

    # ------------------------------------------------------------------
    # LD decay
    # ------------------------------------------------------------------
    ld_info: dict[str, Any] = {"available": False}
    rng = np.random.default_rng(int(cfg["calibration"].get("seed", 20261023)))
    max_ld_snps = 8000
    subset_idx = rng.choice(len(bim), size=min(max_ld_snps, len(bim)), replace=False)
    ld_subset = data_dir / "ld_subset.txt"
    bim.iloc[np.sort(subset_idx)]["SNP"].to_csv(ld_subset, index=False, header=False,
                                                encoding="utf-8")
    ld_prefix = tables / "ld"
    try:
        plink.run(["--bfile", str(clean), *chr_flags, "--autosome",
                   "--extract", str(ld_subset), "--r2",
                   "--ld-window", "100",
                   "--ld-window-kb", str(cfg["popgen"]["ld_window_kb"]),
                   "--ld-window-r2", "0"],
                  stage="popgen", name="ld", out_prefix=ld_prefix, expect=[".ld"])
        ld = pd.read_csv(ld_prefix.with_suffix(".ld"), sep=r"\s+")
        ld.columns = [str(c).strip() for c in ld.columns]
        if {"BP_A", "BP_B", "R2"} <= set(ld.columns):
            ld["distance_kb"] = (pd.to_numeric(ld["BP_B"], errors="coerce") -
                                 pd.to_numeric(ld["BP_A"], errors="coerce")) / 1000.0
            ld["R2"] = pd.to_numeric(ld["R2"], errors="coerce")
            ld = ld.loc[(ld["distance_kb"] >= 0) & ld["R2"].notna()]
            ld["bin"] = pd.cut(ld["distance_kb"], bins=LD_BIN_EDGES_KB, right=False)
            decay = (ld.groupby("bin", observed=True)
                     .agg(distance_low_kb=("distance_kb", "min"),
                          distance_high_kb=("distance_kb", "max"),
                          mean_r2=("R2", "mean"), n_pairs=("R2", "size"))
                     .reset_index(drop=True))
            decay.to_csv(tables / "ld_decay.csv", index=False, encoding="utf-8-sig")
            tables_map["ld_decay"] = "tables/ld_decay.csv"
            add_fig("ld_decay", _try_plot("plot_ld_decay", decay, plots_dir / "ld_decay.png"))
            r2_50 = ld.loc[(ld["distance_kb"] >= 40) & (ld["distance_kb"] <= 60), "R2"]
            ld_info = {
                "available": True, "n_pairs": int(len(ld)),
                "mean_r2": round(float(ld["R2"].mean()), 5),
                "r2_at_50kb": round(float(r2_50.mean()), 5) if len(r2_50) else None,
                "decay": [{"bin_low_kb": float(row.distance_low_kb),
                           "bin_high_kb": float(row.distance_high_kb),
                           "mean_r2": round(float(row.mean_r2), 5),
                           "n_pairs": int(row.n_pairs)}
                          for row in decay.itertuples()],
            }
            if ld_info.get("r2_at_50kb") and ld_info["r2_at_50kb"] > 0.5:
                warnings.append(f"50kb 处平均 r²={ld_info['r2_at_50kb']:.3f} 偏高，提示群体较小或存在结构")
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"LD 衰减计算失败: {exc}")

    # ------------------------------------------------------------------
    # ROH / FROH
    # ------------------------------------------------------------------
    roh_cfg = cfg["popgen"]["roh"]
    roh_prefix = tables / "roh"
    roh_info: dict[str, Any] = {"available": False}
    try:
        plink.run(["--bfile", str(clean), *chr_flags, "--autosome", "--homozyg",
                   "--homozyg-snp", str(roh_cfg["snp"]),
                   "--homozyg-kb", str(roh_cfg["kb"]),
                   "--homozyg-density", str(roh_cfg["density"]),
                   "--homozyg-gap", str(roh_cfg["gap"])],
                  stage="popgen", name="roh", out_prefix=roh_prefix,
                  expect=[".hom.indiv"])
        roh = pd.read_csv(roh_prefix.with_suffix(".hom.indiv"), sep=r"\s+")
        roh.columns = [str(c).strip() for c in roh.columns]
        genome_kb = _genome_length_kb(bim)
        roh["KB"] = pd.to_numeric(roh.get("KB"), errors="coerce").fillna(0.0)
        roh["NSEG"] = pd.to_numeric(roh.get("NSEG"), errors="coerce").fillna(0).astype(int)
        roh["FROH"] = roh["KB"] / genome_kb if genome_kb > 0 else 0.0
        roh_out = roh.rename(columns={"NSEG": "n_roh", "KB": "total_kb"})[
            ["FID", "IID", "n_roh", "total_kb", "FROH"]]
        roh_out.to_csv(tables / "roh_per_sample.csv", index=False, encoding="utf-8-sig")
        tables_map["roh_per_sample"] = "tables/roh_per_sample.csv"
        add_fig("roh_per_sample", _try_plot("plot_roh_per_sample", roh_out,
                                            plots_dir / "roh.png"))
        top = roh_out.sort_values("FROH", ascending=False).head(10)
        roh_info = {"available": True, "genome_length_kb": round(genome_kb, 1),
                    "mean_roh_per_indiv": round(float(roh_out["n_roh"].mean()), 3),
                    "mean_froh": round(float(roh_out["FROH"].mean()), 5),
                    "max_froh": round(float(roh_out["FROH"].max()), 5),
                    "top_samples": [{"iid": str(r.IID), "froh": round(float(r.FROH), 4),
                                     "n_roh": int(r.n_roh), "total_kb": round(float(r.total_kb), 1)}
                                    for r in top.itertuples()]}
        if roh_info["mean_froh"] > 0.2:
            warnings.append(f"平均 FROH={roh_info['mean_froh']:.3f} 偏高，提示近交程度较高")
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"ROH 计算失败: {exc}")

    # ------------------------------------------------------------------
    # FST (requires groups)
    # ------------------------------------------------------------------
    fst_info: dict[str, Any] = {"available": False}
    if group_present:
        within = run_dir / "00_intake" / "data" / "within.txt"
        fst_prefix = tables / "fst"
        try:
            plink.run(["--bfile", str(clean), *chr_flags, "--autosome",
                       "--fst", "--within", str(within)],
                      stage="popgen", name="fst", out_prefix=fst_prefix,
                      expect=[".fst"])
            fst = pd.read_csv(fst_prefix.with_suffix(".fst"), sep=r"\s+")
            fst.columns = [str(c).strip() for c in fst.columns]
            fst["FST"] = pd.to_numeric(fst["FST"], errors="coerce")
            fst.to_csv(tables / "fst.csv", index=False, encoding="utf-8-sig")
            tables_map["fst"] = "tables/fst.csv"
            add_fig("fst", _try_plot("plot_fst_manhattan", fst, plots_dir / "fst.png"))
            log_text = (fst_prefix.with_suffix(".log")).read_text(encoding="utf-8",
                                                                  errors="replace") \
                if fst_prefix.with_suffix(".log").exists() else ""
            mean_match = _re_search(r"Mean Fst estimate:\s*(-?[\d.eE-]+)", log_text)
            weighted_match = _re_search(r"Weighted Fst estimate:\s*(-?[\d.eE-]+)", log_text)
            top = fst.sort_values("FST", ascending=False).head(50)
            fst_info = {
                "available": True,
                "n_variants": int(len(fst)),
                "mean_fst": round(float(mean_match), 5) if mean_match is not None else
                round(float(fst["FST"].mean()), 5),
                "weighted_fst": round(float(weighted_match), 5) if weighted_match is not None else None,
                "levels": cfg.get("_derived", {}).get("group_levels", []),
                "top": [{"snp": str(r.SNP), "chr": str(r.CHR),
                         "bp": int(float(r.POS)) if str(r.POS).replace('.', '').isdigit() else None,
                         "fst": round(float(r.FST), 5)} for r in top.itertuples()],
            }
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"FST 计算失败（已跳过）: {exc}")
            fst_info = {"available": False, "reason": str(exc)}
    else:
        fst_info = {"available": False, "reason": "未提供分组文件"}

    # ------------------------------------------------------------------
    # observed heterozygosity per group
    # ------------------------------------------------------------------
    diversity: dict[str, Any] = {"available": False}
    het_prefix = tables / "het"
    try:
        plink.run(["--bfile", str(clean), *chr_flags, "--autosome", "--het"],
                  stage="popgen", name="het", out_prefix=het_prefix, expect=[".het"])
        het = utils.read_het(het_prefix.with_suffix(".het"))
        hom = pd.to_numeric(het.get("O(HOM)"), errors="coerce")
        nnm = pd.to_numeric(het.get("N(NM)"), errors="coerce")
        het["observed_het"] = 1 - hom / nnm
        if groups is not None:
            het["GROUP"] = het["IID"].astype(str).map(groups).fillna("NA")
            per_group = het.groupby("GROUP")["observed_het"].agg(["mean", "std", "count"]).reset_index()
            per_group.to_csv(tables / "diversity_by_group.csv", index=False, encoding="utf-8-sig")
            tables_map["diversity_by_group"] = "tables/diversity_by_group.csv"
            diversity = {"available": True,
                         "mean_het": round(float(het["observed_het"].mean()), 5),
                         "by_group": [{"group": str(r.GROUP),
                                       "mean_het": round(float(r.mean), 5),
                                       "n": int(r.count)} for r in per_group.itertuples()]}
        else:
            diversity = {"available": True,
                         "mean_het": round(float(het["observed_het"].mean()), 5),
                         "by_group": []}
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"杂合度统计失败: {exc}")

    # ------------------------------------------------------------------
    # metrics
    # ------------------------------------------------------------------
    gates: list[dict[str, Any]] = []

    def gate(name: str, value: float, warn: float, fail: float, higher_is_better: bool,
             detail: str) -> None:
        if higher_is_better:
            status = "pass" if value >= warn else ("warn" if value >= fail else "fail")
        else:
            status = "pass" if value <= warn else ("warn" if value <= fail else "fail")
        gates.append({"name": name, "status": status, "value": round(float(value), 5),
                      "threshold": {"warn": warn, "fail": fail}, "detail": detail})

    if relatedness_summary.get("mean_pihat") is not None:
        gate("mean_relatedness", relatedness_summary["mean_pihat"], 0.15, 0.30, False,
             f"平均 PI_HAT={relatedness_summary['mean_pihat']}；高亲缘会增加 GWAS 假阳性风险")
    if kinship_summary.get("mean_offdiag") is not None:
        gate("mean_kinship", kinship_summary["mean_offdiag"], 0.10, 0.20, False,
             f"平均亲缘系数={kinship_summary['mean_offdiag']}")
    gate("pca_outliers", utils.safe_div(len(pca_outliers), max(n_samples, 1)), 0.02, 0.05,
         False, f"PCA 离群样本 {len(pca_outliers)}")
    if roh_info.get("available"):
        gate("mean_froh", roh_info["mean_froh"], 0.15, 0.30, False,
             f"平均 FROH={roh_info['mean_froh']}")

    overall = "fail" if any(g["status"] == "fail" for g in gates) else (
        "warn" if any(g["status"] == "warn" for g in gates) else "pass")

    metrics = {
        "stage": "popgen",
        "status": overall,
        "n_samples": int(n_samples),
        "n_variants": int(len(bim)),
        "pca": {"n_pcs": int(len(eigenval)),
                "var_explained": [round(float(v), 5) for v in var_explained[:10]],
                "n_outliers": int(len(pca_outliers)), "outliers": pca_outliers[:50]},
        "kinship": kinship_summary,
        "relatedness": relatedness_summary,
        "ld": ld_info,
        "roh": roh_info,
        "fst": fst_info,
        "diversity": diversity,
        "gates": gates,
        "tables": {k: v for k, v in tables_map.items() if v},
        "figures": figures,
        "plot_errors": plot_errors,
        "warnings": warnings + [g["detail"] for g in gates if g["status"] == "warn"],
        "errors": [g["detail"] for g in gates if g["status"] == "fail"],
    }
    metrics["summary"] = _summary(metrics)
    utils.write_json(stage / "metrics.json", metrics)
    ctx.logger.info("群体结构分析完成，状态 %s", overall)
    return metrics


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _genome_length_kb(bim: pd.DataFrame) -> float:
    total_bp = 0.0
    for _, frame in bim.groupby("CHR"):
        bp = pd.to_numeric(frame["BP"], errors="coerce").dropna()
        if len(bp):
            total_bp += float(bp.max() - min(0.0, float(bp.min())))
    return total_bp / 1000.0


def _re_search(pattern: str, text: str) -> float | None:
    import re
    match = re.search(pattern, text)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _summary(metrics: dict[str, Any]) -> str:
    pca = metrics["pca"]
    rel = metrics["relatedness"]
    ld = metrics.get("ld", {})
    roh = metrics.get("roh", {})
    fst = metrics.get("fst", {})
    lines = [
        f"群体结构：PCA 前 {min(5, len(pca['var_explained']))} 个主成分解释 "
        f"{sum(pca['var_explained'][:5]):.1%} 方差；离群样本 {pca['n_outliers']} 个。",
    ]
    if rel.get("mean_pihat") is not None:
        lines.append(f"亲缘：PI_HAT 均值 {rel['mean_pihat']}，最大值 {rel.get('max_pihat')}，"
                     f"≥0.25 的对数 {rel.get('n_pairs_gt_025', 0)}。")
    if ld.get("available"):
        lines.append(f"LD 衰减：50kb 处平均 r²={ld.get('r2_at_50kb')}（{ld['n_pairs']} 对）。")
    if roh.get("available"):
        lines.append(f"ROH：平均每个个体 {roh['mean_roh_per_indiv']} 段，平均 FROH="
                     f"{roh['mean_froh']}（最大值 {roh['max_froh']}）。")
    if fst.get("available"):
        lines.append(f"FST：全局均值 {fst.get('mean_fst')}（加权 {fst.get('weighted_fst')}）。")
    else:
        lines.append(f"FST：未计算（{fst.get('reason', '无分组')}）。")
    failed = [g["name"] for g in metrics["gates"] if g["status"] == "fail"]
    warned = [g["name"] for g in metrics["gates"] if g["status"] == "warn"]
    if failed:
        lines.append("群体结构警告项：" + "、".join(failed) + "。")
    elif warned:
        lines.append("群体结构提示项：" + "、".join(warned) + "。")
    return "\n".join(lines)
