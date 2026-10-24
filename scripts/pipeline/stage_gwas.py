"""Stage 3: genome-wide association analysis.

Runs PLINK --linear (quantitative) or --logistic (binary) with covariates,
then adds allele frequencies, BH q-values, lambda-GC, lead SNPs, LD clumping
and Manhattan/QQ figures.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import utils

try:
    from . import plots
except Exception:  # noqa: BLE001
    plots = None  # type: ignore[assignment]


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
    plink: utils.Plink = ctx.plink
    figures: list[str] = []
    plot_errors: list[str] = []
    warnings: list[str] = []
    tables_map: dict[str, str] = {}

    if not cfg.get("_derived", {}).get("phenotype_present"):
        metrics = _skipped("未提供表型文件，跳过 GWAS（可补充 --pheno 后重跑）")
        utils.write_json(stage / "metrics.json", metrics)
        return metrics

    clean = run_dir / "01_qc" / "data" / "clean"
    pheno_path = run_dir / "00_intake" / "data" / "pheno_clean.txt"
    pca_covar_path = run_dir / "01_qc" / "data" / "pca_covar.txt"
    if not clean.with_suffix(".bed").exists():
        raise FileNotFoundError(f"缺少 QC 清理后的数据: {clean}.bed，请先运行 qc 阶段")
    if not pheno_path.exists():
        raise FileNotFoundError(f"缺少表型文件: {pheno_path}，请先运行 intake 阶段")

    fam = utils.read_fam(clean.with_suffix(".fam"))
    pheno = pd.read_csv(pheno_path, sep=r"\s+", header=None, names=["FID", "IID", "PHENO"],
                        dtype={"FID": str, "IID": str})
    pheno["PHENO"] = pd.to_numeric(pheno["PHENO"], errors="coerce")
    pheno_type = cfg.get("_derived", {}).get("phenotype_type", "quantitative")
    model = str(cfg["gwas"].get("model", "auto")).lower()
    if model in {"auto", ""}:
        model = "logistic" if pheno_type == "binary" else "linear"
    if model == "logistic" and pheno_type != "binary":
        warnings.append("配置要求 logistic，但表型不是二元性状；已改用 linear")
        model = "linear"
    if model == "linear" and pheno_type == "binary":
        warnings.append("二元表型使用 linear 模型不符合规范；已改用 logistic")
        model = "logistic"

    # ---- sample set ----------------------------------------------------
    analysis = fam[["FID", "IID"]].merge(pheno, on=["FID", "IID"], how="left")
    analysis = analysis.loc[analysis["PHENO"].notna() & (analysis["PHENO"] != -9)].copy()

    covar_cols: list[str] = []
    covar_frame = analysis[["FID", "IID"]].copy()
    user_covar_path = run_dir / "00_intake" / "data" / "covar_clean.txt"
    user_covar_columns = list(cfg.get("_derived", {}).get("covariate_columns") or [])
    if user_covar_path.exists() and user_covar_columns:
        user_covar = pd.read_csv(user_covar_path, sep=r"\s+", header=None,
                                 names=["FID", "IID", *user_covar_columns],
                                 dtype={"FID": str, "IID": str})
        covar_frame = covar_frame.merge(user_covar, on=["FID", "IID"], how="left")
        covar_cols.extend([c for c in user_covar_columns if c in covar_frame.columns])
    n_pcs = int(cfg["gwas"].get("n_pcs", 5))
    if pca_covar_path.exists() and n_pcs > 0:
        pca_cols = [f"PC{i + 1}" for i in range(n_pcs)]
        pca = pd.read_csv(pca_covar_path, sep=r"\s+", header=None,
                          names=["FID", "IID", *pca_cols], dtype={"FID": str, "IID": str})
        pca = pca[["FID", "IID", *[c for c in pca_cols if c in pca.columns]]]
        covar_frame = covar_frame.merge(pca, on=["FID", "IID"], how="left")
        covar_cols.extend([c for c in pca_cols if c in covar_frame.columns])

    analysis = analysis.merge(covar_frame[["FID", "IID", *covar_cols]], on=["FID", "IID"], how="left")
    for col in covar_cols:
        analysis[col] = pd.to_numeric(analysis[col], errors="coerce")
        if analysis[col].isna().any():
            analysis[col] = analysis[col].fillna(analysis[col].mean())
    keep_df = analysis[["FID", "IID"]]
    keep_path = data_dir / "analysis_samples.txt"
    keep_df.to_csv(keep_path, sep="\t", header=False, index=False, encoding="utf-8")

    covar_path = data_dir / "covar.txt"
    if covar_cols:
        analysis[["FID", "IID", *covar_cols]].to_csv(covar_path, sep="\t", header=False,
                                                     index=False, encoding="utf-8")
    else:
        covar_path = None  # type: ignore[assignment]

    n_samples = len(analysis)
    if n_samples < 20:
        warnings.append(f"进入关联分析的样本只有 {n_samples} 个，结果不稳定")

    # ---- run PLINK association -----------------------------------------
    assoc_prefix = tables / "gwas"
    args = ["--bfile", str(clean), *chr_flags, "--autosome",
            "--pheno", str(pheno_path), "--keep", str(keep_path),
            "--allow-no-sex", "--keep-allele-order",
            "--ci", str(cfg["gwas"].get("ci", 0.95)), "--adjust"]
    if covar_path is not None:
        args += ["--covar", str(covar_path), "--covar-number", f"1-{len(covar_cols)}"]
    if model == "logistic":
        args.append("--logistic")
        assoc_suffix = ".assoc.logistic"
    else:
        args.append("--linear")
        assoc_suffix = ".assoc.linear"
    plink.run(args, stage="gwas", name=f"gwas_{model}", out_prefix=assoc_prefix,
              expect=[assoc_suffix, assoc_suffix + ".adjusted"])
    assoc_file = assoc_prefix.with_suffix(assoc_suffix)
    adjusted_file = Path(str(assoc_file) + ".adjusted")

    assoc = utils.read_assoc(assoc_file)
    assoc.columns = [str(c).strip() for c in assoc.columns]
    assoc = assoc.loc[assoc["TEST"].astype(str).str.upper() == "ADD"].copy()
    for col in ("BETA", "OR", "SE", "STAT", "P", "BP", "NMISS"):
        if col in assoc.columns:
            assoc[col] = pd.to_numeric(assoc[col], errors="coerce")
    if model == "logistic":
        assoc["BETA"] = np.log(pd.to_numeric(assoc["OR"], errors="coerce"))
    assoc["SE"] = pd.to_numeric(assoc.get("SE"), errors="coerce")

    # allele frequencies from the cleaned data
    frq_path = run_dir / "01_qc" / "tables" / "final_stats.frq"
    bim = utils.read_bim(clean.with_suffix(".bim"))
    if frq_path.exists():
        freq = pd.read_csv(frq_path, sep=r"\s+")
        freq.columns = [str(c).strip().lstrip("#") for c in freq.columns]
        freq["MAF"] = pd.to_numeric(freq["MAF"], errors="coerce")
        freq = freq.merge(bim[["SNP", "A1"]].rename(columns={"A1": "BIM_A1"}), on="SNP", how="left")
        freq["AF"] = np.where(freq["A1"].astype(str) == freq["BIM_A1"].astype(str),
                              freq["MAF"], 1 - freq["MAF"])
        assoc = assoc.merge(freq[["SNP", "AF", "MAF"]].rename(columns={"MAF": "MAF_CLEAN"}),
                            on="SNP", how="left")
    else:
        warnings.append("未找到 .frq 文件，结果表缺少等位基因频率")

    assoc["Q"] = utils.bh_fdr(assoc["P"].fillna(1).to_numpy())
    assoc["P_PLOT"] = assoc["P"].clip(lower=1e-300)
    n_tested = int(len(assoc))
    bonferroni = 0.05 / max(n_tested, 1)
    suggestive = 1.0 / max(n_tested, 1)
    lambda_val = utils.lambda_gc(assoc["P_PLOT"].dropna().tolist())

    assoc_sorted = assoc.sort_values("P")
    lead = assoc_sorted.head(10).copy()
    lead_rows = [{
        "snp": str(row.SNP), "chr": str(row.CHR), "bp": _int_or_none(row.BP),
        "a1": str(row.A1), "af": _round_or_none(getattr(row, "AF", None)),
        "beta": _round_or_none(row.BETA), "se": _round_or_none(row.SE),
        "p": _sig_or_none(row.P), "q": _sig_or_none(row.Q),
        "p_text": _format_p(row.P), "n": _int_or_none(row.NMISS),
    } for row in lead.itertuples()]

    results_cols = [c for c in ("CHR", "SNP", "BP", "A1", "AF", "BETA", "SE", "STAT", "P", "Q", "NMISS")
                    if c in assoc.columns]
    results = assoc[results_cols].copy()
    results.to_csv(tables / "gwas_results.csv", index=False, encoding="utf-8-sig")
    tables_map["gwas_results"] = "tables/gwas_results.csv"
    assoc_sorted.head(100).to_csv(tables / "gwas_top100.csv", index=False, encoding="utf-8-sig")
    tables_map["gwas_top100"] = "tables/gwas_top100.csv"

    # lambda comparison with PLINK's own estimate
    plink_lambda = _parse_lambda_from_log(Path(str(assoc_prefix) + ".log"))

    # ---- clump ---------------------------------------------------------
    clump_info: dict[str, Any] = {"available": False, "n_clumps": 0, "clumps": []}
    clump_cfg = cfg["gwas"].get("clump", {})
    clump_prefix = tables / "clump"
    p1 = float(clump_cfg.get("p1", 5e-8))
    has_significant = bool((assoc["P"] <= bonferroni).any())
    if not has_significant:
        p1 = max(p1, float(suggestive))
        warnings.append(f"无全基因组显著位点，clump 阈值放宽为 {p1:.2e}（仅作信号整理）")
    try:
        plink.run(["--bfile", str(clean), *chr_flags, "--autosome",
                   "--clump", str(assoc_file),
                   "--clump-snp-field", "SNP", "--clump-field", "P",
                   "--clump-p1", f"{p1:g}", "--clump-p2", f"{float(clump_cfg.get('p2', 1e-4)):g}",
                   "--clump-r2", str(clump_cfg.get("r2", 0.2)),
                   "--clump-kb", str(clump_cfg.get("kb", 1000))],
                  stage="gwas", name="clump", out_prefix=clump_prefix,
                  check=False)
        clumped_path = clump_prefix.with_suffix(".clumped")
        clumps: list[dict[str, Any]] = []
        if clumped_path.exists():
            lines = [line for line in clumped_path.read_text(encoding="utf-8", errors="replace").splitlines()
                     if line.strip() and not line.startswith("CHR")]
            for line in lines:
                parts = line.split()
                if len(parts) < 5:
                    continue
                clumps.append({"chr": parts[0], "snp": parts[1], "bp": _int_or_none(parts[2]),
                               "p": _float_or_none(parts[3]), "n_snps": _int_or_none(parts[4]),
                               "span_kb": _float_or_none(parts[5]) if len(parts) > 5 else None})
            if clumps:
                pd.DataFrame(clumps).to_csv(tables / "clumps.csv", index=False,
                                            encoding="utf-8-sig")
                tables_map["clumps"] = "tables/clumps.csv"
            clump_info = {"available": True, "n_clumps": len(clumps),
                          "p1": p1, "clumps": clumps[:50]}
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"clump 失败（不影响主结果）: {exc}")

    # ---- thresholds / metrics -----------------------------------------
    n_significant = int((assoc["P"] <= bonferroni).sum())
    n_fdr05 = int((assoc["Q"] <= 0.05).sum())
    n_suggestive = int((assoc["P"] <= suggestive).sum())
    n_cases = cfg.get("_derived", {}).get("phenotype_type") == "binary" and \
        int((analysis["PHENO"] == 2).sum()) or None
    n_controls = cfg.get("_derived", {}).get("phenotype_type") == "binary" and \
        int((analysis["PHENO"] == 1).sum()) or None

    gates: list[dict[str, Any]] = []

    def gate(name: str, status: str, detail: str, value: Any = None) -> None:
        gates.append({"name": name, "status": status, "value": value, "detail": detail})

    if n_tested == 0:
        gates.append({"name": "variants_tested", "status": "fail", "value": 0,
                      "detail": "没有任何位点进入关联分析"})
    if np.isfinite(lambda_val):
        if 0.9 <= lambda_val <= 1.1:
            gate("lambda_gc", "pass", f"λGC={lambda_val:.3f}，在可接受范围 [0.9, 1.1]", round(lambda_val, 3))
        elif 0.8 <= lambda_val <= 1.3:
            gate("lambda_gc", "warn", f"λGC={lambda_val:.3f} 偏离 1（0.8-1.3），提示结构/亲缘或检验膨胀", round(lambda_val, 3))
        else:
            gate("lambda_gc", "fail", f"λGC={lambda_val:.3f} 严重偏离 1，结果不可信", round(lambda_val, 3))
    if pheno_type == "binary" and n_cases and n_controls:
        if min(n_cases, n_controls) < 10:
            gate("case_control_balance", "warn", f"病例 {n_cases} / 对照 {n_controls}，样本量偏小")
        else:
            gate("case_control_balance", "pass", f"病例 {n_cases} / 对照 {n_controls}")
    if n_samples < 100:
        gate("sample_size", "warn", f"关联分析样本量 {n_samples} 偏小，统计效力有限")
    else:
        gate("sample_size", "pass", f"关联分析样本量 {n_samples}")

    overall = "fail" if any(g["status"] == "fail" for g in gates) else (
        "warn" if any(g["status"] == "warn" for g in gates) else "pass")

    add_fig = figures.append
    if assoc["P_PLOT"].notna().any():
        plot_frame = assoc.copy()
        plot_frame["P"] = plot_frame["P_PLOT"]
        manhattan = _try_plot("plot_manhattan", plot_frame, plots_dir / "manhattan.png",
                              bonferroni=bonferroni)
        qq = _try_plot("plot_qq", plot_frame["P"], plots_dir / "qq.png", lambda_gc=lambda_val)
        for name, path in (("manhattan", manhattan), ("qq", qq)):
            if path and path.startswith("__error__:"):
                plot_errors.append(f"{name}: {path.split(':', 1)[1]}")
            elif path:
                add_fig(path)

    metrics = {
        "stage": "gwas",
        "status": overall,
        "model": model,
        "phenotype_type": pheno_type,
        "n_samples": int(n_samples),
        "n_tested": int(n_tested),
        "n_cases": n_cases, "n_controls": n_controls,
        "covariates": [c for c in covar_cols if not c.startswith("PC")] +
                      [c for c in covar_cols if c.startswith("PC")],
        "covar_file": str(covar_path) if covar_path else None,
        "assoc_file": str(assoc_file),
        "adjusted_file": str(adjusted_file),
        "lambda_gc": round(float(lambda_val), 4) if np.isfinite(lambda_val) else None,
        "lambda_plink": plink_lambda,
        "lambda_source": "p-values (median chi-square)",
        "thresholds": {"bonferroni": bonferroni, "suggestive": suggestive,
                       "fdr_0.05": 0.05 if n_fdr05 else None},
        "n_significant": n_significant, "n_suggestive": n_suggestive, "n_fdr05": n_fdr05,
        "lead_snps": lead_rows,
        "clumps": clump_info,
        "inflation": {"mean_chi2": _mean_chi2(assoc["STAT"]), "lambda": round(float(lambda_val), 4)
                      if np.isfinite(lambda_val) else None},
        "p_distribution": {"min": _round_or_none(assoc["P"].min()),
                           "min_text": _format_p(assoc["P"].min()),
                           "median": _round_or_none(assoc["P"].median()),
                           "n_p_lt_5e8": int((assoc["P"] < 5e-8).sum())},
        "gates": gates,
        "tables": tables_map,
        "figures": figures,
        "plot_errors": plot_errors,
        "warnings": warnings + [g["detail"] for g in gates if g["status"] == "warn"],
        "errors": [g["detail"] for g in gates if g["status"] == "fail"],
    }
    metrics["summary"] = _summary(metrics)
    utils.write_json(stage / "metrics.json", metrics)
    ctx.logger.info("GWAS 完成：%s，测试 %d 个位点，显著 %d，λGC=%s，状态 %s",
                    model, n_tested, n_significant,
                    metrics["lambda_gc"], overall)
    return metrics


# ---------------------------------------------------------------------------
def _skipped(reason: str) -> dict[str, Any]:
    return {"stage": "gwas", "status": "skipped", "summary": f"GWAS 未执行：{reason}",
            "model": None, "n_samples": 0, "n_tested": 0, "covariates": [],
            "lambda_gc": None, "n_significant": 0, "lead_snps": [], "gates": [],
            "tables": {}, "figures": [], "warnings": [reason], "errors": []}


def _parse_lambda_from_log(path: Path) -> float | None:
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"lambda[^=]*=\s*([0-9]+(?:\.[0-9]+)?)", text, re.IGNORECASE)
    return float(match.group(1)) if match else None


def _format_p(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(number):
        return "NA"
    if number == 0:
        return "<1e-300"
    return f"{number:.3e}"


def _mean_chi2(stats: pd.Series) -> float | None:
    values = pd.to_numeric(stats, errors="coerce").dropna()
    if values.empty:
        return None
    return round(float((values ** 2).mean()), 4)


def _int_or_none(value: Any) -> int | None:
    try:
        number = float(value)
        return int(number) if np.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _float_or_none(value: Any) -> float | None:
    try:
        number = float(value)
        return round(number, 6) if np.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _round_or_none(value: Any, digits: int = 6) -> float | None:
    try:
        number = float(value)
        return round(number, digits) if np.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _sig_or_none(value: Any, digits: int = 6) -> float | None:
    """Round to significant digits so tiny p-values (1e-12) survive JSON."""
    try:
        number = float(value)
        if not np.isfinite(number):
            return None
        return float(f"{number:.{digits}g}")
    except (TypeError, ValueError):
        return None


def _summary(metrics: dict[str, Any]) -> str:
    lines = [
        f"GWAS（{metrics['model']}）完成：{metrics['n_samples']} 个样本、{metrics['n_tested']} 个位点，"
        f"协变量 {metrics['covariates'] or '无'}。",
        f"多重检验：Bonferroni 阈值 {metrics['thresholds']['bonferroni']:.2e}，"
        f"显著位点 {metrics['n_significant']} 个，FDR<0.05 位点 {metrics['n_fdr05']} 个，"
        f"建议阈值 {metrics['thresholds']['suggestive']:.2e} 下 {metrics['n_suggestive']} 个。",
        f"基因组膨胀因子 λGC={metrics['lambda_gc']}。",
    ]
    if metrics["lead_snps"]:
        top = metrics["lead_snps"][0]
        p_text = top.get("p_text") or (f"{top['p']:.2e}" if top.get("p") is not None else "NA")
        lines.append(f"最强信号：{top['snp']}（{top['chr']}:{top['bp']}，P={p_text}，"
                     f"效应={top['beta']}）。")
    if metrics["clumps"].get("available"):
        lines.append(f"LD clump：{metrics['clumps']['n_clumps']} 个独立信号区间。")
    failed = [g["name"] for g in metrics["gates"] if g["status"] == "fail"]
    warned = [g["name"] for g in metrics["gates"] if g["status"] == "warn"]
    if failed:
        lines.append("警告：关联结果存在严重问题（" + "、".join(failed) + "），请勿直接下结论。")
    elif warned:
        lines.append("提示：" + "、".join(warned) + "（详见 03_gwas/report.md）。")
    return "\n".join(lines)
