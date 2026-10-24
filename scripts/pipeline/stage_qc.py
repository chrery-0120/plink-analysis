"""Stage 1: strict QC of samples and variants.

Order of operations (mirrors PLINK's own filter semantics):

1. compute variant-level statistics on the raw fileset (autosomes)
2. remove variants failing missingness / MAF / HWE thresholds
3. on the variant-filtered data compute sample missingness, heterozygosity,
   sex check, LD pruning, IBD and PCA
4. remove samples failing missingness / heterozygosity / sex / duplicate /
   (optionally) relatedness / PCA outlier rules
5. recompute all statistics on the cleaned fileset and emit tables, figures,
   gates and metrics
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import utils

try:  # figures are optional so QC still works if the reporting module is absent
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
        path = func(*args, **kwargs)
        return str(path)
    except Exception as exc:  # noqa: BLE001 - a failed plot must not fail QC
        return f"__error__:{exc}"


def run(ctx) -> dict[str, Any]:
    cfg = ctx.config
    qc = cfg["qc"]
    stage = ctx.stage_dir
    tables, plots_dir, data_dir = stage / "tables", stage / "plots", stage / "data"
    base = utils.resolve_base_prefix(ctx.run_dir)
    chr_flags = list(cfg.get("_derived", {}).get("chr_flags") or [])
    plink: utils.Plink = ctx.plink
    figures: list[str] = []
    plot_errors: list[str] = []

    # ------------------------------------------------------------------
    # 1. variant-level statistics on the raw data
    # ------------------------------------------------------------------
    pheno_path = ctx.run_dir / "00_intake" / "data" / "pheno_clean.txt"
    binary_pheno = (cfg.get("_derived", {}).get("phenotype_type") == "binary"
                    and pheno_path.exists())
    stats_prefix = tables / "raw_stats"
    stats_args = ["--bfile", str(base), *chr_flags, "--autosome",
                  "--freq", "--missing", "--hardy", "midp", "--het"]
    if binary_pheno:
        # PLINK --hwe semantics: for case/control data HWE is tested in controls
        # only, so stratify --hardy by phenotype instead of testing all samples.
        stats_args += ["--pheno", str(pheno_path)]
    plink.run(stats_args, stage="qc", name="raw_stats", out_prefix=stats_prefix,
              expect=[".frq", ".lmiss", ".imiss", ".hwe", ".het"])

    frq = _read_plink_table(stats_prefix.with_suffix(".frq"))
    lmiss = _read_plink_table(stats_prefix.with_suffix(".lmiss"))
    imiss = _read_plink_table(stats_prefix.with_suffix(".imiss"))
    hwe = _read_plink_table(stats_prefix.with_suffix(".hwe"))
    het = _read_plink_table(stats_prefix.with_suffix(".het"))

    hwe_all = _hwe_for_filter(hwe, binary=binary_pheno)
    variant_table = (frq[["CHR", "SNP", "A1", "A2", "MAF"]]
                     .merge(lmiss[["SNP", "F_MISS"]], on="SNP", how="left")
                     .merge(hwe_all[["SNP", "P"]].rename(columns={"P": "HWE_P"}), on="SNP", how="left"))
    variant_table["MAF"] = pd.to_numeric(variant_table["MAF"], errors="coerce")
    variant_table["F_MISS"] = pd.to_numeric(variant_table["F_MISS"], errors="coerce")
    variant_table["HWE_P"] = pd.to_numeric(variant_table["HWE_P"], errors="coerce")

    # ------------------------------------------------------------------
    # 2. variant filters
    # ------------------------------------------------------------------
    fail_missing = variant_table["F_MISS"] > float(qc["geno"])
    fail_maf = variant_table["MAF"] < float(qc["maf"])
    fail_hwe = variant_table["HWE_P"] < float(qc["hwe"])
    exclude = variant_table.loc[fail_missing | fail_maf | fail_hwe, "SNP"].astype(str)
    exclude_path = data_dir / "exclude_variants.txt"
    exclude.to_csv(exclude_path, index=False, header=False, encoding="utf-8")

    qc_vars = data_dir / "qc_variants"
    plink.run(["--bfile", str(base), *chr_flags, "--autosome",
               "--exclude", str(exclude_path), "--make-bed", "--keep-allele-order"],
              stage="qc", name="apply_variant_filters", out_prefix=qc_vars,
              expect=[".bed", ".bim", ".fam"])
    n_variants_after = len(utils.read_bim(qc_vars.with_suffix(".bim")))

    # ------------------------------------------------------------------
    # 3. sample-level statistics on the variant-filtered data
    # ------------------------------------------------------------------
    sample_prefix = tables / "sample_stats"
    plink.run(["--bfile", str(qc_vars), *chr_flags, "--autosome",
               "--missing", "--het", "--freq"],
              stage="qc", name="sample_stats", out_prefix=sample_prefix,
              expect=[".imiss", ".het", ".frq"])
    imiss2 = _read_plink_table(sample_prefix.with_suffix(".imiss"))
    het2 = _read_plink_table(sample_prefix.with_suffix(".het"))

    # LD pruning (needed for PCA / relatedness / LD estimates)
    prune_prefix = data_dir / "pruned"
    plink.run(["--bfile", str(qc_vars), *chr_flags, "--autosome",
               "--indep-pairwise", str(qc["ld_prune"]["window"]),
               str(qc["ld_prune"]["step"]), str(qc["ld_prune"]["r2"])],
              stage="qc", name="ld_prune", out_prefix=prune_prefix,
              expect=[".prune.in", ".prune.out"])
    prune_in = prune_prefix.with_suffix(".prune.in")
    n_pruned = _count_lines(prune_in)

    # PCA on pruned variants (used for outlier detection and GWAS covariates)
    pca_prefix = tables / "pca"
    n_pcs = max(int(cfg["popgen"]["pca"]), int(cfg["gwas"].get("n_pcs", 5)) + 2)
    plink.run(["--bfile", str(qc_vars), *chr_flags, "--autosome",
               "--extract", str(prune_in), "--pca", str(n_pcs)],
              stage="qc", name="pca", out_prefix=pca_prefix,
              expect=[".eigenvec", ".eigenval"])
    eigenvec = utils.read_eigenvec(pca_prefix.with_suffix(".eigenvec"))
    eigenval = utils.read_eigenval(pca_prefix.with_suffix(".eigenval"))
    var_explained = _variance_explained(eigenval, n_samples=len(imiss2))

    # IBD / relatedness on pruned variants
    genome_prefix = tables / "genome"
    genome_args = ["--bfile", str(qc_vars), *chr_flags, "--autosome",
                   "--extract", str(prune_in), "--genome"]
    if len(imiss2) > 2000:
        genome_args += ["--min", "0.1"]
    plink.run(genome_args, stage="qc", name="genome", out_prefix=genome_prefix,
              expect=[".genome"])
    genome = utils.read_csv(genome_prefix.with_suffix(".genome"), sep=r"\s+")
    if "PI_HAT" in genome.columns:
        genome["PI_HAT"] = pd.to_numeric(genome["PI_HAT"], errors="coerce")

    # sex check (XY species only, and only when the X chromosome exists)
    sex_prefix = tables / "sexcheck"
    sex_checked = False
    sex_problems: list[dict[str, Any]] = []
    species = cfg.get("species", {})
    chromosomes = set(map(str, cfg.get("_derived", {}).get("chromosomes", [])))
    if (str(species.get("sex_system", "")).upper() == "XY" and "X" in chromosomes
            and str(qc.get("sex_check", "auto")) != "false"):
        try:
            plink.run(["--bfile", str(qc_vars), *chr_flags, "--check-sex"],
                      stage="qc", name="check_sex", out_prefix=sex_prefix,
                      expect=[".sexcheck"])
            sex_df = _read_plink_table(sex_prefix.with_suffix(".sexcheck"))
            sex_df["F"] = pd.to_numeric(sex_df.get("F"), errors="coerce")
            if "STATUS" in sex_df.columns:
                sex_df["STATUS"] = sex_df["STATUS"].astype(str)
                sex_problems = sex_df.loc[sex_df["STATUS"].str.upper().str.contains("PROBLEM"),
                                          ["FID", "IID", "PEDSEX", "SNPSEX", "F"]].to_dict("records")
            sex_checked = True
        except utils.PlinkError as exc:
            ctx.logger.warning("性别检查跳过：%s", exc)
            plot_errors.append(f"sex_check skipped: {exc}")

    # ------------------------------------------------------------------
    # 4. sample filters
    # ------------------------------------------------------------------
    fam = utils.read_fam(qc_vars.with_suffix(".fam"))
    imiss2 = imiss2.merge(fam[["FID", "IID"]], on=["FID", "IID"], how="right")
    imiss2["F_MISS"] = pd.to_numeric(imiss2["F_MISS"], errors="coerce").fillna(0.0)
    het2 = het2.merge(fam[["FID", "IID"]], on=["FID", "IID"], how="right")
    het2["F"] = pd.to_numeric(het2["F"], errors="coerce")

    remove_reasons: dict[tuple[str, str], list[str]] = {}
    def mark(frame: pd.DataFrame, reason: str) -> None:
        for _, row in frame.iterrows():
            remove_reasons.setdefault((str(row["FID"]), str(row["IID"])), []).append(reason)

    mark(imiss2.loc[imiss2["F_MISS"] > float(qc["mind"])], "missingness")

    f_mean, f_sd = float(het2["F"].mean()), float(het2["F"].std(ddof=1))
    if np.isfinite(f_sd) and f_sd > 0:
        low, high = f_mean - float(qc["het_sd"]) * f_sd, f_mean + float(qc["het_sd"]) * f_sd
        mark(het2.loc[(het2["F"] < low) | (het2["F"] > high)], "heterozygosity_outlier")
    else:
        low, high = float("nan"), float("nan")

    if sex_checked and sex_problems:
        mark(pd.DataFrame(sex_problems), "sex_discordant")

    # duplicates / related pairs: keep the sample with lower missingness
    miss_lookup = {(str(r.FID), str(r.IID)): float(r.F_MISS) for r in imiss2.itertuples()}
    related_pairs: list[dict[str, Any]] = []
    dup_remove: set[tuple[str, str]] = set()
    if "PI_HAT" in genome.columns:
        self_pairs = ((genome["FID1"].astype(str) == genome["FID2"].astype(str)) &
                      (genome["IID1"].astype(str) == genome["IID2"].astype(str)))
        genome = genome.loc[~self_pairs]
        for row in genome.itertuples():
            key1 = (str(row.FID1), str(row.IID1))
            key2 = (str(row.FID2), str(row.IID2))
            pihat = float(getattr(row, "PI_HAT", np.nan))
            if not np.isfinite(pihat):
                continue
            if pihat >= float(qc["related_pihat"]):
                related_pairs.append({"FID1": key1[0], "IID1": key1[1],
                                      "FID2": key2[0], "IID2": key2[1], "PI_HAT": pihat})
            if pihat >= float(qc["duplicate_pihat"]):
                drop = key1 if miss_lookup.get(key1, 0) > miss_lookup.get(key2, 0) else key2
                dup_remove.add(drop)
    for key in dup_remove:
        remove_reasons.setdefault(key, []).append("duplicate")

    related_remove: set[tuple[str, str]] = set()
    if qc.get("remove_related"):
        for pair in sorted(related_pairs, key=lambda p: -p["PI_HAT"]):
            k1 = (pair["FID1"], pair["IID1"])
            k2 = (pair["FID2"], pair["IID2"])
            if k1 in related_remove or k2 in related_remove:
                continue
            drop = k1 if miss_lookup.get(k1, 0) > miss_lookup.get(k2, 0) else k2
            related_remove.add(drop)
        for key in related_remove:
            remove_reasons.setdefault(key, []).append("related")

    # PCA outliers (flag; removal only when configured)
    pca_outliers: list[str] = []
    pca_flags: list[str] = []
    if len(eigenvec) > 3 and len(eigenval) > 1:
        pc_cols = [c for c in eigenvec.columns if c.startswith("PC")][:2]
        outlier_mask = pd.Series(False, index=eigenvec.index)
        for col in pc_cols:
            values = pd.to_numeric(eigenvec[col], errors="coerce")
            sd = float(values.std(ddof=1))
            if np.isfinite(sd) and sd > 0:
                outlier_mask |= (values - values.mean()).abs() > float(qc["pca_outlier_sd"]) * sd
        pca_outliers = eigenvec.loc[outlier_mask, "IID"].astype(str).tolist()
        pca_flags = pca_outliers
        if qc.get("remove_pca_outliers"):
            for iid in pca_outliers:
                match = eigenvec.loc[eigenvec["IID"].astype(str) == iid]
                if not match.empty:
                    key = (str(match.iloc[0]["FID"]), iid)
                    remove_reasons.setdefault(key, []).append("pca_outlier")

    remove_df = pd.DataFrame([
        {"FID": key[0], "IID": key[1], "reasons": ";".join(sorted(set(reasons)))}
        for key, reasons in remove_reasons.items()
    ])
    if remove_df.empty:
        remove_df = pd.DataFrame(columns=["FID", "IID", "reasons"])
    remove_df.to_csv(data_dir / "remove_samples.txt", sep="\t", header=False, index=False,
                     encoding="utf-8")
    utils.write_csv(remove_df, tables / "removed_samples.csv")

    clean = data_dir / "clean"
    plink.run(["--bfile", str(qc_vars), *chr_flags,
               "--remove", str(data_dir / "remove_samples.txt"),
               "--make-bed", "--keep-allele-order"],
              stage="qc", name="apply_sample_filters", out_prefix=clean,
              expect=[".bed", ".bim", ".fam"])

    # ------------------------------------------------------------------
    # 5. final statistics on the cleaned data
    # ------------------------------------------------------------------
    final_prefix = tables / "final_stats"
    final_args = ["--bfile", str(clean), *chr_flags, "--autosome",
                  "--freq", "--missing", "--hardy", "midp", "--het"]
    if binary_pheno:
        final_args += ["--pheno", str(pheno_path)]
    plink.run(final_args, stage="qc", name="final_stats", out_prefix=final_prefix,
              expect=[".frq", ".lmiss", ".imiss", ".hwe", ".het"])
    final_frq = _read_plink_table(final_prefix.with_suffix(".frq"))
    final_lmiss = _read_plink_table(final_prefix.with_suffix(".lmiss"))
    final_imiss = _read_plink_table(final_prefix.with_suffix(".imiss"))
    final_hwe = _hwe_for_filter(_read_plink_table(final_prefix.with_suffix(".hwe")),
                                binary=binary_pheno)
    final_het = _read_plink_table(final_prefix.with_suffix(".het"))

    # recompute pruning + PCA on clean data for downstream analyses
    clean_prune = data_dir / "clean_pruned"
    plink.run(["--bfile", str(clean), *chr_flags, "--autosome",
               "--indep-pairwise", str(qc["ld_prune"]["window"]),
               str(qc["ld_prune"]["step"]), str(qc["ld_prune"]["r2"])],
              stage="qc", name="clean_ld_prune", out_prefix=clean_prune,
              expect=[".prune.in", ".prune.out"])
    clean_prune_in = clean_prune.with_suffix(".prune.in")
    clean_pca = tables / "clean_pca"
    plink.run(["--bfile", str(clean), *chr_flags, "--autosome",
               "--extract", str(clean_prune_in), "--pca", str(n_pcs)],
              stage="qc", name="clean_pca", out_prefix=clean_pca,
              expect=[".eigenvec", ".eigenval"])
    clean_eigenvec = utils.read_eigenvec(clean_pca.with_suffix(".eigenvec"))
    clean_eigenval = utils.read_eigenval(clean_pca.with_suffix(".eigenval"))
    clean_var_explained = _variance_explained(clean_eigenval, n_samples=len(clean_eigenvec))

    # PC covariates for GWAS (PC1..k, no header, PLINK --covar format)
    n_cov_pcs = int(cfg["gwas"].get("n_pcs", 5))
    pc_cols = [f"PC{i + 1}" for i in range(min(n_cov_pcs, clean_eigenvec.shape[1] - 2))]
    pca_covar = clean_eigenvec[["FID", "IID", *pc_cols]].copy()
    pca_covar.to_csv(data_dir / "pca_covar.txt", sep="\t", header=False, index=False,
                     encoding="utf-8")

    # ------------------------------------------------------------------
    # 6. tables / figures
    # ------------------------------------------------------------------
    variant_table.to_csv(tables / "variant_qc_raw.csv", index=False, encoding="utf-8-sig")
    final_variant = (final_frq[["CHR", "SNP", "A1", "A2", "MAF"]]
                     .merge(final_lmiss[["SNP", "F_MISS"]], on="SNP", how="left")
                     .merge(final_hwe[["SNP", "P"]].rename(columns={"P": "HWE_P"}), on="SNP", how="left"))
    final_variant["MAF"] = pd.to_numeric(final_variant["MAF"], errors="coerce")
    final_variant["F_MISS"] = pd.to_numeric(final_variant["F_MISS"], errors="coerce")
    final_variant["HWE_P"] = pd.to_numeric(final_variant["HWE_P"], errors="coerce")
    final_variant.to_csv(tables / "variant_qc_clean.csv", index=False, encoding="utf-8-sig")
    final_imiss.to_csv(tables / "sample_qc_clean.csv", index=False, encoding="utf-8-sig")
    final_het.to_csv(tables / "het_clean.csv", index=False, encoding="utf-8-sig")
    if not remove_df.empty:
        remove_df.to_csv(tables / "removed_samples.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([{"SNP": str(s)} for s in exclude]).to_csv(
        tables / "removed_variants.csv", index=False, encoding="utf-8-sig")
    if related_pairs:
        pd.DataFrame(related_pairs).to_csv(tables / "related_pairs.csv", index=False,
                                           encoding="utf-8-sig")
    clean_eigenvec.to_csv(tables / "pca.eigenvec", index=False, encoding="utf-8-sig")
    pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(len(clean_eigenval))],
                  "eigenvalue": clean_eigenval,
                  "variance_explained": clean_var_explained}
                 ).to_csv(tables / "pca.eigenval.csv", index=False, encoding="utf-8-sig")

    # figures
    def add_fig(name: str, path: str | None) -> None:
        if path is None:
            return
        if path.startswith("__error__:"):
            plot_errors.append(f"{name}: {path.split(':', 1)[1]}")
        else:
            figures.append(path)

    add_fig("maf_spectrum", _try_plot("plot_maf_spectrum", final_frq["MAF"].astype(float),
                                      plots_dir / "maf_spectrum.png"))
    add_fig("missingness", _try_plot("plot_missingness", final_imiss["F_MISS"].astype(float),
                                     final_lmiss["F_MISS"].astype(float),
                                     plots_dir / "missingness.png"))
    add_fig("het_distribution", _try_plot("plot_het_distribution", final_het,
                                          plots_dir / "heterozygosity.png"))
    add_fig("hwe_qq", _try_plot("plot_hwe_qq", final_hwe["P"].astype(float),
                                plots_dir / "hwe_qq.png"))
    add_fig("pca_scatter", _try_plot("plot_pca_scatter", clean_eigenvec, list(clean_eigenval),
                                     plots_dir / "pca.png", n_samples=len(clean_eigenvec)))
    add_fig("pca_scree", _try_plot("plot_pca_scree", list(clean_eigenval),
                                   plots_dir / "pca_scree.png", n_samples=len(clean_eigenvec)))
    add_fig("relatedness", _try_plot("plot_relatedness_hist",
                                     pd.to_numeric(genome.get("PI_HAT", pd.Series(dtype=float)),
                                                   errors="coerce"),
                                     plots_dir / "relatedness.png"))

    # ------------------------------------------------------------------
    # 7. gates + metrics
    # ------------------------------------------------------------------
    n_samples_in, n_variants_in = len(fam), len(variant_table)
    n_samples_out, n_variants_out = len(final_imiss), len(final_frq)
    sample_ret = utils.safe_div(n_samples_out, n_samples_in)
    variant_ret = utils.safe_div(n_variants_out, n_variants_in)
    hwe_fail_frac = utils.safe_div(int((final_variant["HWE_P"] < 1e-6).sum()), max(n_variants_out, 1))
    mean_miss = float(pd.to_numeric(final_lmiss["F_MISS"], errors="coerce").mean())
    mean_maf = float(pd.to_numeric(final_frq["MAF"], errors="coerce").mean())
    sex_disc = len(sex_problems)
    n_related = len(related_pairs)
    remove_counts = remove_df["reasons"].str.split(";").explode().value_counts().to_dict() if not remove_df.empty else {}

    gates: list[dict[str, Any]] = []

    def gate(name: str, value: float, warn: float, fail: float, higher_is_better: bool,
             detail: str) -> None:
        if higher_is_better:
            status = "pass" if value >= warn else ("warn" if value >= fail else "fail")
        else:
            status = "pass" if value <= warn else ("warn" if value <= fail else "fail")
        gates.append({"name": name, "status": status, "value": round(float(value), 5),
                      "threshold": {"warn": warn, "fail": fail}, "detail": detail})

    gate("sample_retention", sample_ret, 0.90, 0.60, True,
         f"保留 {n_samples_out}/{n_samples_in} 个样本")
    gate("variant_retention", variant_ret, 0.80, 0.30, True,
         f"保留 {n_variants_out}/{n_variants_in} 个位点（常染色体）")
    gate("genotype_call_rate", 1 - mean_miss, 0.97, 0.95, True,
         f"清理后平均位点缺失率 {mean_miss:.4f}")
    gate("hwe_fail_fraction", hwe_fail_frac, 0.01, 0.05, False,
         f"HWE p<1e-6 的位点比例 {hwe_fail_frac:.4f}")
    general_het_outliers = int(remove_counts.get("heterozygosity_outlier", 0))
    gate("het_outlier_fraction", utils.safe_div(general_het_outliers, n_samples_in), 0.02, 0.05,
         False, f"杂合度离群样本 {general_het_outliers}")
    if sex_checked:
        gate("sex_discordance", utils.safe_div(sex_disc, n_samples_in), 0.01, 0.03, False,
             f"性别不一致样本 {sex_disc}")
    gates.append({"name": "duplicate_samples", "status": "pass", "value": int(len(dup_remove)),
                  "threshold": {"duplicate_pihat": float(qc["duplicate_pihat"])},
                  "detail": (f"已识别并剔除 {len(dup_remove)} 个重复样本"
                             f"（PI_HAT≥{qc['duplicate_pihat']}；数据质量问题，已处理）")})
    gate("pca_outliers", utils.safe_div(len(pca_outliers), n_samples_out), 0.02, 0.05, False,
         f"PCA 离群样本 {len(pca_outliers)}（默认仅标记不删除）")

    overall = "fail" if any(g["status"] == "fail" for g in gates) else (
        "warn" if any(g["status"] == "warn" for g in gates) else "pass")

    metrics = {
        "stage": "qc",
        "status": overall,
        "n_samples_in": int(n_samples_in), "n_samples_out": int(n_samples_out),
        "n_variants_in": int(n_variants_in), "n_variants_out": int(n_variants_out),
        "removed": {
            "samples": {k: int(v) for k, v in remove_counts.items()},
            "samples_total": int(len(remove_df)),
            "variants": {
                "missingness": int(fail_missing.sum()),
                "maf": int(fail_maf.sum()),
                "hwe": int(fail_hwe.sum()),
                "total": int(len(exclude)),
            },
        },
        "missing_rate": {"sample_mean": round(mean_miss, 6),
                         "variant_mean": round(float(pd.to_numeric(final_lmiss["F_MISS"],
                                                                   errors="coerce").mean()), 6)},
        "maf": {"mean": round(mean_maf, 5),
                "n_below_001": int((pd.to_numeric(final_frq["MAF"], errors="coerce") < 0.01).sum())},
        "hwe": {"n_p_lt_1e-6": int((final_variant["HWE_P"] < 1e-6).sum()),
                "fraction": round(hwe_fail_frac, 5)},
        "het": {"mean": round(float(pd.to_numeric(final_het["F"], errors="coerce").mean()), 5),
                "sd": round(float(pd.to_numeric(final_het["F"], errors="coerce").std(ddof=1)), 5),
                "sd_thresholds": [None if not np.isfinite(low) else round(low, 4),
                                  None if not np.isfinite(high) else round(high, 4)]},
        "sex": {"checked": sex_checked, "n_discordant": int(sex_disc),
                "discordant": sex_problems[:20],
                "thresholds": {"F_male_max": 0.2, "F_female_min": 0.8}},
        "relatedness": {"n_duplicates": int(len(dup_remove)),
                        "n_related_pairs": int(n_related),
                        "max_pihat": round(float(genome["PI_HAT"].max()), 4) if n_related and "PI_HAT" in genome else None,
                        "mean_pihat": round(float(pd.to_numeric(genome.get("PI_HAT", pd.Series(dtype=float)),
                                                                errors="coerce").mean()), 4)},
        "pca": {"n_pcs": int(len(clean_eigenval)),
                "var_explained": [round(float(v), 5) for v in clean_var_explained[:10]],
                "n_outliers": int(len(pca_outliers)), "outliers": pca_outliers[:50],
                "ld_pruned_variants": int(n_pruned)},
        "pruned": {"n_variants": int(_count_lines(clean_prune_in))},
        "gates": gates,
        "tables": {
            "variant_qc_raw": "tables/variant_qc_raw.csv",
            "variant_qc_clean": "tables/variant_qc_clean.csv",
            "sample_qc_clean": "tables/sample_qc_clean.csv",
            "het_clean": "tables/het_clean.csv",
            "removed_samples": "tables/removed_samples.csv" if not remove_df.empty else None,
            "removed_variants": "tables/removed_variants.csv",
            "related_pairs": "tables/related_pairs.csv" if related_pairs else None,
            "pca_covar": "data/pca_covar.txt",
        },
        "figures": figures,
        "plot_errors": plot_errors,
        "warnings": [g["detail"] for g in gates if g["status"] == "warn"],
        "errors": [g["detail"] for g in gates if g["status"] == "fail"],
        "summary": _summary(n_samples_in, n_samples_out, n_variants_in, n_variants_out,
                            remove_counts, mean_maf, mean_miss, hwe_fail_frac,
                            sex_checked, sex_disc, n_related, len(pca_outliers), gates),
    }
    utils.write_json(stage / "metrics.json", metrics)
    ctx.logger.info("QC 完成：%d→%d 样本，%d→%d 位点，状态 %s",
                    n_samples_in, n_samples_out, n_variants_in, n_variants_out, overall)
    return metrics


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _read_plink_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", dtype=str)
    df.columns = [str(c).strip().lstrip("#") for c in df.columns]
    return df


def _count_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip())


def _hwe_for_filter(hwe: pd.DataFrame, *, binary: bool = False) -> pd.DataFrame:
    """Pick the HWE rows used for filtering.

    For case/control data the standard is to test HWE in controls only; PLINK's
    ``--hardy`` then emits TEST=ALL/AFF/UNAFF rows, and we keep UNAFF.
    """
    if "TEST" in hwe.columns:
        tests = hwe["TEST"].astype(str).str.upper()
        preferred = "UNAFF" if binary else "ALL"
        rows = hwe.loc[tests == preferred]
        if not rows.empty:
            return rows
        rows = hwe.loc[tests == "ALL"]
        if not rows.empty:
            return rows
    return hwe


def _variance_explained(eigenval: np.ndarray, n_samples: int | None = None) -> np.ndarray:
    """PLINK 1.9 eigenvectors come from the variance-standardised GRM, whose
    eigenvalues sum to the *sample count* (trace of the GRM), not to the sum of
    the retained leading eigenvalues.  Proportion of variance explained is
    therefore eigenvalue / N."""
    total = float(n_samples) if n_samples else float(np.sum(eigenval))
    if total <= 0:
        return np.zeros_like(eigenval)
    return eigenval / total


def _summary(n_samples_in: int, n_samples_out: int, n_variants_in: int, n_variants_out: int,
             remove_counts: dict[str, int], mean_maf: float, mean_miss: float,
             hwe_fail_frac: float, sex_checked: bool, sex_disc: int, n_related: int,
             n_pca_outliers: int, gates: list[dict[str, Any]]) -> str:
    lines = [
        f"质控完成：样本 {n_samples_in}→{n_samples_out}（保留 {n_samples_out / max(n_samples_in, 1):.1%}），"
        f"常染色体位点 {n_variants_in}→{n_variants_out}（保留 {n_variants_out / max(n_variants_in, 1):.1%}）。",
        f"清理后平均 MAF={mean_maf:.3f}，平均位点缺失率={mean_miss:.4f}，"
        f"HWE p<1e-6 比例={hwe_fail_frac:.4f}。",
    ]
    if remove_counts:
        parts = "、".join(f"{k} {v}" for k, v in sorted(remove_counts.items()))
        lines.append(f"样本剔除原因：{parts}。")
    if sex_checked:
        lines.append(f"性别检查：不一致 {sex_disc} 个样本。")
    lines.append(f"高亲缘对（PI_HAT≥0.25）：{n_related} 对；PCA 离群样本：{n_pca_outliers} 个。")
    failed = [g["name"] for g in gates if g["status"] == "fail"]
    warned = [g["name"] for g in gates if g["status"] == "warn"]
    if failed:
        lines.append("质控未通过项：" + "、".join(failed) + "。建议调整阈值或检查数据来源后重跑。")
    elif warned:
        lines.append("质控警告项：" + "、".join(warned) + "（详见 01_qc/report.md）。")
    else:
        lines.append("所有质控门限均通过。")
    return "\n".join(lines)
