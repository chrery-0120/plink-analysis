"""Self-test for the report-eng deliverables (plots.py / report.py / report.css).

Run from anywhere with::

    python scripts/pipeline/_report_selftest.py

The script builds a temporary, fully fabricated PLINK run directory (all five
stages, fake metrics per DESIGN.md section 8, real CSV tables and real PNG
figures produced by plots.py), then exercises:

* ``plots.*`` on realistic synthetic data and on degenerate input (missing
  columns / empty / all-NaN) -- placeholder figures must never raise;
* ``report.stage_report`` for all five stages (plus one skipped-stage run);
* ``report.final_report`` -- Markdown, self-contained HTML with base64 images,
  and machine-readable summary.json;
* ``report.render_markdown_to_html`` on its own with tables, code, lists,
  blockquotes, links and images.

Every check prints PASS/FAIL; the temporary run directory is removed at the
end.  The script itself is kept (it is the only file created in pipeline/ that
is not part of the public module API).
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import plots  # noqa: E402
import report  # noqa: E402

CHECKS: list[tuple[str, bool, str]] = []


def check(name: str, condition: object, detail: str = "") -> bool:
    ok = bool(condition)
    CHECKS.append((name, ok, detail))
    marker = "PASS" if ok else "FAIL"
    suffix = f"  [{detail}]" if detail else ""
    print(f"{marker}  {name}{suffix}", flush=True)
    return ok


def expect_raises_free(name: str, func, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception as exc:  # pragma: no cover - this is the failure path
        CHECKS.append((name, False, f"raised {type(exc).__name__}: {exc}"))
        print(f"FAIL  {name}  [{type(exc).__name__}: {exc}]", flush=True)
        traceback.print_exc()
        return None


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
CONFIG: dict = {
    "schema_version": 1,
    "run_id": "20261023_2200_selftest",
    "created_at": "2026-10-23T22:00:00+08:00",
    "input": {
        "genotype": {"type": "bed", "prefix": "demo_data/pig", "path": None},
        "phenotype": {
            "type": "quantitative",
            "path": "demo_data/pig.pheno",
            "missing": ["-9", "NA", "0"],
        },
        "covariates": {
            "path": "demo_data/pig.covar",
            "columns": ["sex", "age"],
            "categorical": ["sex"],
        },
        "groups": {"path": "demo_data/pig.group", "name": "breed"},
        "pedigree": {"path": None},
    },
    "species": {
        "name": "pig",
        "chr_set": None,
        "autosome_count": 18,
        "sex_system": "XY",
        "nonautosomal": ["X", "Y", "MT"],
    },
    "qc": {
        "mind": 0.10,
        "geno": 0.10,
        "maf": 0.01,
        "hwe": 1e-06,
        "het_sd": 3.0,
        "duplicate_pihat": 0.90,
        "related_pihat": 0.25,
        "sex_check": "auto",
        "pca_outlier_sd": 6.0,
        "ld_prune": {"window": 50, "step": 5, "r2": 0.2},
        "remove_related": False,
        "remove_pca_outliers": False,
    },
    "popgen": {
        "pca": 20,
        "ld_window_kb": 1000,
        "ld_bin_kb": 50,
        "roh": {"snp": 50, "kb": 1000, "density": 50, "gap": 1000, "het_missing": 0.0},
    },
    "gwas": {
        "model": "auto",
        "covariates": "auto",
        "n_pcs": 5,
        "ci": 0.95,
        "maf": 0.01,
        "geno": 0.10,
        "hwe": 1e-06,
        "clump": {"p1": 5e-08, "p2": 1e-04, "r2": 0.2, "kb": 1000},
        "genomic_control": "none",
    },
    "calibration": {
        "enabled": True,
        "n_replicates": 10,
        "n_causal": 20,
        "h2": 0.5,
        "n_snps_subset": 20000,
        "alpha": 5e-08,
    },
    "runtime": {"threads": 4, "memory_mb": 4000, "timeout_s": 7200, "keep_intermediate": True},
    "report": {"language": "zh-CN", "title": "伪造数据自测项目"},
}


def _write_csv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8")


def _make_metadata(rng: np.random.Generator) -> dict[str, dict]:
    n = 400
    lead_snps = []
    for i in range(25):
        lead_snps.append(
            {
                "snp": f"rs{100000 + i}",
                "chr": str(i % 18 + 1),
                "bp": int(rng.integers(10_000, 2_000_000)),
                "a1": str(rng.choice(["A", "C", "G", "T"])),
                "a2": str(rng.choice(["A", "C", "G", "T"])),
                "af": float(rng.uniform(0.05, 0.5)),
                "beta": float(rng.normal(0, 0.4)),
                "se": float(abs(rng.normal(0.05, 0.01))),
                "p": float(10 ** rng.uniform(-12, -5)),
                "n": n,
            }
        )
    lead_snps.sort(key=lambda item: item["p"])
    # PLINK prints p = 0 for underflowed p-values; stage_gwas supplies p_text.
    for index, item in enumerate(lead_snps):
        if index == 0:
            item["p"] = 0.0
            item["p_text"] = "<1e-300"
        else:
            item["p_text"] = f"{item['p']:.3e}"

    intake = {
        "genotype_type": "bed",
        "n_samples": n,
        "n_variants": 8000,
        "chromosomes": [str(i) for i in range(1, 19)] + ["X"],
        "species_inferred": {
            "name": "sus_scrofa",
            "chr_set": "18",
            "sex_system": "XY",
            "autosomes": 18,
        },
        "files": [
            {"path": "demo_data/pig.bed", "size": 3_276_800, "sha256": "a" * 64},
            {"path": "demo_data/pig.bim", "size": 1_048_576, "sha256": "b" * 64},
            {"path": "demo_data/pig.fam", "size": 12_800, "sha256": "c" * 64},
        ],
        "preflight": {
            "status": "pass",
            "checks": [
                {"name": "files_exist", "status": "pass", "detail": "3 个文件均可读"},
                {"name": "sample_id_unique", "status": "pass", "detail": "400 个 IID 唯一"},
                {"name": "phenotype_parse", "status": "pass", "detail": "400 条记录，2 个缺失"},
                {"name": "chromosome_names", "status": "warn", "detail": "检测到 chr 前缀，已规范化"},
            ],
        },
        "phenotype": {"present": True, "type": "quantitative", "n": n, "missing": 2},
        "covariates": {"present": True, "columns": ["sex", "age"]},
        "groups": {"present": True, "levels": [{"name": "Duroc", "n": 150}, {"name": "Landrace", "n": 250}]},
    }
    qc = {
        "n_samples_in": 400,
        "n_samples_out": 390,
        "n_variants_in": 8000,
        "n_variants_out": 6500,
        "removed": {
            "samples": {
                "missingness": 2,
                "het_outlier": 1,
                "sex_discordant": 0,
                "duplicate": 1,
                "related": 3,
                "pca_outlier": 0,
            },
            "variants": {"missingness": 300, "maf": 800, "hwe": 400},
        },
        "missing_rate": {"sample_mean": 0.01, "variant_mean": 0.005},
        "maf": {"mean": 0.24, "n_below_001": 0},
        "hwe": {"n_p_lt_1e-6": 12, "lambda": 0.99},
        "het": {"mean": -0.01, "sd": 0.02},
        "sex": {"checked": True, "discordant": [], "thresholds": {"F": 0.2, "M": 0.8}},
        "relatedness": {"n_duplicates": 1, "n_related": 3, "max_pihat": 0.98, "mean_pihat": 0.02},
        "pca": {"eigenvalues": [1.2, 0.8], "var_explained": [0.05, 0.03], "n_outliers": 0},
        "pruned": {"n_variants": 1200},
        "gates": [
            {"name": "sample_missingness", "status": "pass", "value": 0.01, "threshold": 0.10},
            {"name": "variant_missingness", "status": "pass", "value": 0.005, "threshold": 0.10},
            {"name": "maf", "status": "pass", "value": 0.24, "threshold": 0.01},
            {"name": "hwe", "status": "warn", "value": 12, "threshold": 20},
            {"name": "heterozygosity", "status": "pass", "value": -0.01, "threshold": 3.0},
            {"name": "relatedness", "status": "pass", "value": 3, "threshold": 5},
            {"name": "sex_check", "status": "pass", "value": 0, "threshold": 0},
            {"name": "sample_retention", "status": "pass", "value": 0.975, "threshold": 0.90},
        ],
        "overall": "warn",
    }
    popgen = {
        "pca": {
            "n_pcs": 20,
            "var_explained": [0.052, 0.031, 0.021, 0.017, 0.014, 0.012, 0.011, 0.010, 0.009, 0.008],
            "n_outliers": 2,
            "outliers": ["ID0007", "ID0123"],
        },
        "relatedness": {"n_pairs_gt_025": 5, "max": 0.51},
        "ld": {
            "n_pairs": 123456,
            "bins": [
                {"distance_bp": 50_000, "mean_r2": 0.31, "n_pairs": 12_345},
                {"distance_bp": 100_000, "mean_r2": 0.24, "n_pairs": 11_111},
                {"distance_bp": 150_000, "mean_r2": 0.19, "n_pairs": 10_000},
                {"distance_bp": 200_000, "mean_r2": 0.15, "n_pairs": 9_000},
                {"distance_bp": 250_000, "mean_r2": 0.12, "n_pairs": 8_800},
            ],
        },
        "roh": {
            "n_roh": 412,
            "mean_roh_per_indiv": 2.1,
            "mean_froh": 0.05,
            "top_samples": [
                {"iid": f"ID{i:04d}", "froh": 0.21 - i * 0.01, "n_roh": 8 - i, "total_kb": 12000 - i * 500}
                for i in range(12)
            ],
        },
        "fst": {
            "global": 0.08,
            "n_variants": 6500,
            "top": [
                {"snp": f"rs{i}", "chr": str(i % 18 + 1), "bp": 123 + i, "fst": 0.4 - i * 0.01}
                for i in range(12)
            ],
        },
        "diversity": {"mean_het": 0.28, "n_groups": 2},
        "overall": "pass",
    }
    gwas = {
        "model": "linear",
        "phenotype_type": "quantitative",
        "n_samples": 390,
        "n_tested": 6400,
        "covariates": ["PC1", "PC2", "PC3", "PC4", "PC5", "sex"],
        "lambda_gc": 1.01,
        "lambda_source": "p-values",
        "thresholds": {"bonferroni": 7.8e-06, "fdr_0.05": None, "suggestive": 1.56e-04},
        "n_significant": 3,
        "n_fdr05": 12,
        "n_suggestive": 40,
        "lead_snps": lead_snps,
        "p_distribution": {"min": 0.0, "min_text": "<1e-300", "median": 0.125, "n_p_lt_5e8": 3},
        "inflation": {"mean_chi2": 1.0, "lambda": 1.01},
        "overall": "pass",
    }
    validation = {
        "checks": [
            {
                "name": name,
                "status": "pass" if index not in {5, 9} else "warn",
                "detail": detail,
                "value": value,
                "threshold": threshold,
            }
            for index, (name, detail, value, threshold) in enumerate(
                [
                    ("sample_count", "fam 行数与 manifest 一致", 400, 400),
                    ("variant_count", "bim 行数与 manifest 一致", 8000, 8000),
                    ("maf_recompute", "20 个 SNP 最大绝对误差", 1e-07, 1e-05),
                    ("missingness_recompute", "缺失率最大绝对误差", 2e-07, 1e-05),
                    ("duplicate_ids", "重复 IID 数量", 0, 0),
                    ("monomorphic_sites", "单态位点数量", 3, 10),
                    ("phenotype_encoding", "缺失/异常编码数量", 2, 10),
                    ("pca_eigenvalues", "特征值和偏差", 1e-06, 1e-04),
                    ("kinship_range", "PI_HAT 取值范围", 0.98, 1.0),
                    ("lambda_recompute", "λGC 绝对误差", 0.004, 0.01),
                    ("beta_se_statsmodels", "beta 最大绝对误差", 3e-05, 1e-04),
                    ("hwe_extreme_fraction", "HWE p<1e-6 比例", 0.0018, 0.01),
                    ("allele_frequency_range", "AF 范围", 0.5, 1.0),
                    ("report_consistency", "metrics 与表格计数一致", 1, 1),
                ]
            )
        ],
        "overall": "pass",
    }
    calibration = {
        "null": {
            "n_replicates": 10,
            "lambda_mean": 1.0,
            "lambda_sd": 0.02,
            "type1": {"alpha_0.05": 0.049, "alpha_1e-3": 0.0011, "alpha_5e-8": 0.0},
        },
        "power": {
            "h2": 0.5,
            "n_causal": 20,
            "alpha_5e-8": {"mean": 0.72, "sd": 0.08},
            "alpha_1e-4": {"mean": 0.95},
        },
        "effect_bias": {"mean_beta_ratio": 1.02, "se_calibration": 0.98},
        "concordance": {"top_causal_in_top100": 0.95, "rank_auc": 0.97},
        "validation": {
            "n_checks": 14,
            "n_pass": 12,
            "n_warn": 2,
            "n_fail": 0,
            "overall": "pass",
        },
        "overall": "pass",
    }
    return {
        "intake": intake,
        "qc": qc,
        "popgen": popgen,
        "gwas": gwas,
        "calibration": calibration,
        "validation": validation,
    }


def _make_stage_plots(stage_dir: Path, stage: str, rng: np.random.Generator) -> None:
    plot_dir = stage_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    n = 400
    if stage == "intake":
        plots.plot_phenotype_distribution(pd.Series(rng.normal(100, 15, n)), plot_dir / "phenotype_distribution.png")
        return
    if stage == "qc":
        plots.plot_maf_spectrum(pd.Series(np.clip(rng.beta(1.2, 4.0, 6000), 0.001, 0.5)), plot_dir / "maf_spectrum.png")
        plots.plot_missingness(
            pd.Series(rng.beta(0.3, 30, n)),
            pd.Series(rng.beta(0.2, 40, 6000)),
            plot_dir / "missingness.png",
            mind=CONFIG["qc"]["mind"],
            geno=CONFIG["qc"]["geno"],
        )
        plots.plot_het_distribution(
            pd.DataFrame({"FID": [f"F{i}" for i in range(n)], "IID": [f"ID{i:04d}" for i in range(n)], "F": rng.normal(-0.01, 0.02, n)}),
            plot_dir / "het_distribution.png",
        )
        plots.plot_hwe_qq(pd.Series(rng.uniform(0, 1, 6000)), plot_dir / "hwe_qq.png")
        plots.plot_relatedness_hist(
            pd.Series(np.concatenate([rng.beta(0.4, 12, 2000), rng.uniform(0.2, 1.0, 5)])),
            plot_dir / "relatedness_hist.png",
        )
        return
    if stage == "popgen":
        eig = np.array([1.2, 0.8, 0.5, 0.3, 0.2, 0.1, 0.08, 0.06, 0.05, 0.04])
        eigenvec = pd.DataFrame(
            {
                "FID": [f"F{i}" for i in range(n)],
                "IID": [f"ID{i:04d}" for i in range(n)],
                **{f"PC{k + 1}": rng.normal(size=n) for k in range(5)},
            }
        )
        groups = pd.DataFrame(
            {"IID": eigenvec["IID"], "breed": rng.choice(["Duroc", "Landrace", "Pietrain"], n)}
        )
        plots.plot_pca_scatter(eigenvec, eig, plot_dir / "pca_scatter.png", groups=groups)
        plots.plot_pca_scree(eig, plot_dir / "pca_scree.png", n_samples=len(eigenvec))
        distances = np.arange(50_000, 1_000_001, 50_000)
        plots.plot_ld_decay(
            pd.DataFrame({"distance_bp": distances, "mean_r2": np.exp(-distances / 250_000)}),
            plot_dir / "ld_decay.png",
        )
        plots.plot_roh_per_sample(
            pd.DataFrame({"IID": [f"ID{i:04d}" for i in range(120)], "FROH": rng.beta(1.5, 8, 120)}),
            plot_dir / "roh_per_sample.png",
        )
        plots.plot_fst_manhattan(
            pd.DataFrame(
                {
                    "CHR": rng.integers(1, 19, 3000).astype(str),
                    "BP": rng.integers(1, 2_000_000, 3000),
                    "FST": rng.beta(0.5, 6, 3000),
                }
            ),
            plot_dir / "fst_manhattan.png",
        )
        kinship = rng.normal(0, 0.1, (40, 40))
        kinship = (kinship + kinship.T) / 2
        np.fill_diagonal(kinship, 0.5)
        plots.plot_kinship_heatmap(
            pd.DataFrame(kinship, index=[f"ID{i:04d}" for i in range(40)], columns=[f"ID{i:04d}" for i in range(40)]),
            [f"ID{i:04d}" for i in range(40)],
            plot_dir / "kinship_heatmap.png",
        )
        return
    if stage == "gwas":
        snps = 6000
        p_values = rng.uniform(0, 1, snps) ** 2
        p_values[:5] = 10 ** rng.uniform(-12, -8, 5)
        plots.plot_manhattan(
            pd.DataFrame(
                {
                    "SNP": [f"rs{i}" for i in range(snps)],
                    "CHR": rng.integers(1, 19, snps).astype(str),
                    "BP": rng.integers(1, 2_000_000, snps),
                    "P": p_values,
                }
            ),
            plot_dir / "manhattan.png",
            bonferroni=7.8e-06,
            suggestive=1.56e-04,
        )
        p_qq = np.sort(rng.uniform(0, 1, 3000) ** 1.02)
        plots.plot_qq(pd.Series(p_qq), plot_dir / "qq.png")
        return
    if stage == "calibration":
        plots.plot_calibration(rng.normal(1.0, 0.02, 10), plot_dir / "calibration.png")
        plots.plot_power_curve(
            pd.DataFrame(
                {
                    "alpha": [0.05, 1e-03, 1e-04, 1e-05, 5e-08],
                    "power": [0.99, 0.92, 0.78, 0.55, 0.21],
                }
            ),
            plot_dir / "power_curve.png",
        )
        plots.plot_effect_bias(
            pd.DataFrame(
                {
                    "true_beta": rng.normal(0, 0.2, 100),
                    "est_beta": rng.normal(0, 0.2, 100),
                }
            ),
            plot_dir / "effect_bias.png",
        )


def _make_tables(stage_dir: Path, stage: str, rng: np.random.Generator) -> dict[str, str]:
    tables_dir = stage_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    if stage == "intake":
        frame = pd.DataFrame(
            {
                "check": ["files_exist", "sample_id_unique", "phenotype_parse", "chromosome_names"],
                "status": ["PASS", "PASS", "PASS", "WARN"],
                "detail": ["3 个文件均可读", "400 个 IID 唯一", "400 条记录，2 个缺失", "检测到 chr 前缀"],
            }
        )
        _write_csv(tables_dir / "preflight.csv", frame)
        return {"preflight": "tables/preflight.csv"}
    if stage == "qc":
        sample_qc = pd.DataFrame(
            {
                "IID": [f"ID{i:04d}" for i in range(30)],
                "MISSING_RATE": np.round(rng.uniform(0, 0.08, 30), 5),
                "F": np.round(rng.normal(-0.01, 0.02, 30), 4),
                "OUTLIER": ["no"] * 30,
            }
        )
        _write_csv(tables_dir / "sample_qc.csv", sample_qc)
        summary = pd.DataFrame(
            {
                "metric": ["n_samples_in", "n_samples_out", "n_variants_in", "n_variants_out", "mean_maf"],
                "value": [400, 390, 8000, 6500, 0.24],
            }
        )
        _write_csv(tables_dir / "qc_summary.csv", summary)
        return {"sample_qc": "tables/sample_qc.csv", "qc_summary": "tables/qc_summary.csv"}
    if stage == "popgen":
        distances = np.arange(50_000, 650_000, 50_000)
        frame = pd.DataFrame(
            {
                "distance_bp": distances,
                "mean_r2": np.round(np.exp(-distances / 250_000), 4),
                "n_pairs": rng.integers(5000, 20000, distances.size),
            }
        )
        _write_csv(tables_dir / "ld_bins.csv", frame)
        return {"ld_bins": "tables/ld_bins.csv"}
    if stage == "gwas":
        p_values = ["0.0"] + [f"{v:.3e}" for v in 10 ** rng.uniform(-12, -5, 24)]
        frame = pd.DataFrame(
            {
                "SNP": [f"rs{100000 + i}" for i in range(25)],
                "CHR": [str(i % 18 + 1) for i in range(25)],
                "BP": rng.integers(10_000, 2_000_000, 25),
                "A1": rng.choice(["A", "C", "G", "T"], 25),
                "A2": rng.choice(["A", "C", "G", "T"], 25),
                "AF": np.round(rng.uniform(0.05, 0.5, 25), 3),
                "BETA": np.round(rng.normal(0, 0.4, 25), 4),
                "SE": np.round(np.abs(rng.normal(0.05, 0.01, 25)), 4),
                "P": p_values,
                "P_TEXT": ["<1e-300"] + p_values[1:],
            }
        )
        _write_csv(tables_dir / "lead_snps.csv", frame)
        return {"lead_snps": "tables/lead_snps.csv"}
    if stage == "calibration":
        frame = pd.DataFrame(
            {
                "check": [f"check_{i:02d}" for i in range(14)],
                "status": ["PASS"] * 14,
                "value": [f"{rng.uniform(0, 1):.3e}" for _ in range(14)],
                "threshold": ["1e-05"] * 14,
                "detail": ["独立重算一致"] * 14,
            }
        )
        _write_csv(tables_dir / "validation_checks.csv", frame)
        return {"validation_checks": "tables/validation_checks.csv"}
    return {}


def _make_run(root: Path) -> tuple[Path, dict, dict, dict]:
    """Create a complete fake run directory; return (run_dir, config, state, metrics_by_stage)."""
    rng = np.random.default_rng(20261023)
    run_dir = root / "plink_runs" / CONFIG["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "logs").mkdir(exist_ok=True)
    (run_dir / "tmp").mkdir(exist_ok=True)
    metrics = _make_metadata(rng)
    tables_map: dict[str, dict[str, str]] = {}
    state: dict = {"run_id": CONFIG["run_id"], "updated_at": "2026-10-23T22:10:00+08:00", "stages": {}}
    for index, stage in enumerate(["intake", "qc", "popgen", "gwas", "calibration"]):
        stage_dir = run_dir / report.STAGE_DIRS[stage]
        stage_dir.mkdir(parents=True, exist_ok=True)
        _make_stage_plots(stage_dir, stage, rng)
        tables_map[stage] = _make_tables(stage_dir, stage, rng)
        stage_metrics = metrics[stage]
        if stage == "qc":
            wrapper = {
                "stage": "qc",
                "status": "warn",
                "summary": (
                    "质控完成：400 个样本中保留 390 个（97.5%），8000 个位点中保留 6500 个（81.3%）。"
                    "HWE 偏分离位点略多，已按阈值过滤；整体质量满足下游分析要求，但需关注 HWE 警告。"
                ),
                "metrics": stage_metrics,
                "tables": tables_map[stage],
                "figures": sorted(f"plots/{p.name}" for p in (stage_dir / "plots").glob("*.png")),
                "warnings": ["HWE 检验 p < 1e-6 的位点数 12，接近预设上限 20。"],
                "errors": [],
                "elapsed_s": 12.3,
            }
            (stage_dir / "metrics.json").write_text(
                json.dumps(wrapper, ensure_ascii=False, indent=2), encoding="utf-8"
            )
        else:
            payload = dict(stage_metrics)
            payload["tables"] = tables_map[stage]
            payload["figures"] = sorted(f"plots/{p.name}" for p in (stage_dir / "plots").glob("*.png"))
            payload["status"] = "pass"
            payload["warnings"] = []
            payload["errors"] = []
            payload["summary"] = f"{report._STAGE_META[stage]['name']}已完成，关键指标见报告。"
            (stage_dir / "metrics.json").write_text(
                json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
            )
        if stage == "calibration":
            (stage_dir / "validation.json").write_text(
                json.dumps(metrics["validation"], ensure_ascii=False, indent=2), encoding="utf-8"
            )
        state["stages"][stage] = {
            "status": "warn" if stage == "qc" else "pass",
            "started": "2026-10-23T22:00:00+08:00",
            "finished": "2026-10-23T22:10:00+08:00",
            "elapsed_s": float(5 + index * 7),
            "metrics_path": f"{report.STAGE_DIRS[stage]}/metrics.json",
            "report_path": f"{report.STAGE_DIRS[stage]}/report.md",
            "error": None,
        }
    (run_dir / "run.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")
    (run_dir / "state.json").write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8")
    (run_dir / "logs" / "pipeline.log").write_text("fake pipeline log\n", encoding="utf-8")
    return run_dir, CONFIG, state, metrics


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def _test_plots(root: Path) -> None:
    print("\n--- plots.py ---", flush=True)
    plot_dir = root / "plot_checks"
    plot_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    n = 300

    calls = {
        "maf_spectrum.png": lambda p: plots.plot_maf_spectrum(
            pd.Series(np.clip(rng.beta(1.3, 4, 5000), 0.001, 0.5)), p
        ),
        "missingness.png": lambda p: plots.plot_missingness(
            pd.Series(rng.beta(0.3, 30, n)), pd.Series(rng.beta(0.2, 40, 5000)), p
        ),
        "het_distribution.png": lambda p: plots.plot_het_distribution(
            pd.DataFrame({"F": rng.normal(-0.01, 0.02, n), "IID": [f"ID{i}" for i in range(n)]}), p
        ),
        "hwe_qq.png": lambda p: plots.plot_hwe_qq(pd.Series(rng.uniform(0, 1, 5000)), p),
        "pca_scatter.png": lambda p: plots.plot_pca_scatter(
            pd.DataFrame(
                {
                    "IID": [f"ID{i}" for i in range(n)],
                    **{f"PC{k + 1}": rng.normal(size=n) for k in range(4)},
                }
            ),
            [1.2, 0.8, 0.5, 0.3],
            p,
            groups=pd.Series(rng.choice(["A", "B"], n), index=[f"ID{i}" for i in range(n)]),
        ),
        "pca_scree.png": lambda p: plots.plot_pca_scree([1.2, 0.8, 0.5, 0.3, 0.2], p),
        "pca_scree_n.png": lambda p: plots.plot_pca_scree([6.36, 1.94, 1.86, 1.85, 1.81], p, n_samples=398),
        "relatedness_hist.png": lambda p: plots.plot_relatedness_hist(pd.Series(rng.beta(0.4, 12, 1500)), p),
        "kinship_heatmap.png": lambda p: plots.plot_kinship_heatmap(
            pd.DataFrame(rng.normal(0, 0.1, (20, 20))), [f"ID{i}" for i in range(20)], p
        ),
        "ld_decay.png": lambda p: plots.plot_ld_decay(
            pd.DataFrame({"distance_bp": np.arange(50_000, 600_000, 50_000), "mean_r2": np.linspace(0.4, 0.05, 11)}), p
        ),
        "roh_per_sample.png": lambda p: plots.plot_roh_per_sample(
            pd.DataFrame({"IID": [f"ID{i}" for i in range(80)], "FROH": rng.beta(1.5, 8, 80)}), p
        ),
        "fst_manhattan.png": lambda p: plots.plot_fst_manhattan(
            pd.DataFrame(
                {
                    "CHR": rng.integers(1, 19, 2000).astype(str),
                    "BP": rng.integers(1, 2_000_000, 2000),
                    "FST": rng.beta(0.5, 6, 2000),
                }
            ),
            p,
        ),
        "manhattan.png": lambda p: plots.plot_manhattan(
            pd.DataFrame(
                {
                    "CHR": rng.integers(1, 19, 2000).astype(str),
                    "BP": rng.integers(1, 2_000_000, 2000),
                    "P": rng.uniform(0, 1, 2000) ** 3,
                    "SNP": [f"rs{i}" for i in range(2000)],
                }
            ),
            p,
            bonferroni=7.8e-06,
        ),
        "qq.png": lambda p: plots.plot_qq(pd.Series(np.sort(rng.uniform(0, 1, 3000) ** 1.05)), p),
        "calibration.png": lambda p: plots.plot_calibration(rng.normal(1.0, 0.02, 10), p),
        "power_curve.png": lambda p: plots.plot_power_curve(
            pd.DataFrame({"alpha": [0.05, 1e-3, 1e-4, 5e-8], "power": [0.99, 0.9, 0.5, 0.1]}), p
        ),
        "effect_bias.png": lambda p: plots.plot_effect_bias(
            pd.DataFrame({"true_beta": rng.normal(0, 0.2, 80), "est_beta": rng.normal(0, 0.2, 80)}), p
        ),
        "phenotype_quantitative.png": lambda p: plots.plot_phenotype_distribution(pd.Series(rng.normal(100, 15, n)), p),
        "phenotype_binary.png": lambda p: plots.plot_phenotype_distribution(
            pd.Series(rng.integers(0, 2, n)), p, kind="binary"
        ),
    }
    real_ok = 0
    for filename, func in calls.items():
        target = plot_dir / filename
        result = expect_raises_free(f"plots {filename} does not raise", func, target)
        if result is None:
            continue
        size = target.stat().st_size if target.exists() else 0
        check(f"plots {filename} PNG > 5KB", size > 5000, f"{size / 1024:.1f} KB")
        if size > 5000:
            real_ok += 1
    check(f"plots: all {len(calls)} plot kinds produced a usable PNG", real_ok == len(calls), f"{real_ok}/{len(calls)}")

    # --- PCA variance normalisation regression ---------------------------
    # PLINK 1.9 --pca eigenvalues come from the variance-standardised GRM and
    # sum to ~N, so PC_i% = eigenvalue_i / N * 100 (NOT eigenvalue_i / sum(eig)).
    print("\n--- PCA variance normalisation regression ---", flush=True)
    synthetic = plots.pca_variance_percent([6.36279, 1.93992, 1.85881], n_samples=398)
    expected_pc1 = 6.36279 / 398 * 100
    check(
        "pca_variance_percent(PC1) == eigenvalue / N * 100",
        synthetic.size > 0 and abs(float(synthetic[0]) - expected_pc1) < 1e-9,
        f"{float(synthetic[0]):.4f}% (expected {expected_pc1:.4f}%)",
    )
    check(
        "pca_variance_percent() without N returns empty (no silent sum-normalisation)",
        plots.pca_variance_percent([6.36, 1.94]).size == 0,
    )
    demo_eigenval = Path("plink_runs/demo3/02_popgen/tables/pca.eigenval")
    demo_eigenvec = Path("plink_runs/demo3/02_popgen/tables/pca.eigenvec")
    if demo_eigenval.is_file() and demo_eigenvec.is_file():
        demo_eig = [float(x) for x in demo_eigenval.read_text(encoding="utf-8").split()]
        demo_lines = [line for line in demo_eigenvec.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
        first_line = demo_lines[0].strip() if demo_lines else ""
        if "," in first_line:
            has_header = first_line.split(",")[0].strip().upper() in {"FID", "IID", "#FID"}
        else:
            has_header = bool(first_line) and first_line.split()[0].upper() in {"FID", "IID", "#FID"}
        n_demo = len(demo_lines) - (1 if has_header else 0)
        demo_pct = plots.pca_variance_percent(demo_eig, n_samples=n_demo)
        expected_demo = demo_eig[0] / n_demo * 100.0
        buggy_pct = demo_eig[0] / sum(demo_eig) * 100.0
        check(
            "demo PCA: PC1 uses eigenvalue/N (not eigenvalue/sum)",
            demo_pct.size > 0
            and abs(float(demo_pct[0]) - expected_demo) < 1e-9
            and float(demo_pct[0]) < 5.0
            and buggy_pct > 5.0,
            f"N={n_demo}, PC1={float(demo_pct[0]):.3f}%, old-bug would be {buggy_pct:.2f}%",
        )
        demo_out = plot_dir / "pca_scree_demo_n.png"
        result = expect_raises_free(
            "plots scree with demo eigenvalues + n_samples",
            lambda p: plots.plot_pca_scree(demo_eig, p, n_samples=n_demo),
            demo_out,
        )
        if result is not None:
            size = demo_out.stat().st_size if demo_out.exists() else 0
            check("demo scree PNG > 5KB", size > 5000, f"{size / 1024:.1f} KB")
    else:
        print("SKIP demo scree check (demo3 run not found; synthetic check covers the formula)", flush=True)

    # degenerate inputs -> placeholder PNG, never raise
    degenerate = {
        "empty_series_maf.png": lambda p: plots.plot_maf_spectrum(pd.Series([], dtype=float), p),
        "empty_series_qq.png": lambda p: plots.plot_qq(pd.Series([], dtype=float), p),
        "all_nan_manhattan.png": lambda p: plots.plot_manhattan(
            pd.DataFrame({"CHR": ["1"], "BP": [1], "P": [np.nan]}), p
        ),
        "missing_columns_ld.png": lambda p: plots.plot_ld_decay(pd.DataFrame({"foo": [1]}), p),
        "none_heatmap.png": lambda p: plots.plot_kinship_heatmap(None, None, p),
        "empty_phenotype.png": lambda p: plots.plot_phenotype_distribution(pd.Series([], dtype=float), p),
    }
    placeholder_ok = 0
    for filename, func in degenerate.items():
        target = plot_dir / filename
        result = expect_raises_free(f"plots degenerate {filename} does not raise", func, target)
        if result is None:
            continue
        size = target.stat().st_size if target.exists() else 0
        check(f"plots degenerate {filename} writes placeholder", size > 1000, f"{size / 1024:.1f} KB")
        if size > 1000:
            placeholder_ok += 1
    check("plots: degenerate inputs handled gracefully", placeholder_ok == len(degenerate), f"{placeholder_ok}/{len(degenerate)}")


def _test_stage_reports(run_dir: Path, state: dict, metrics: dict) -> None:
    print("\n--- report.stage_report ---", flush=True)
    expected_sections = ["一、阶段目的", "二、方法与参数", "三、结果", "四、质控判定", "五、异常与处理", "六、对后续分析的影响", "七、产物清单"]
    for stage in report.STAGE_ORDER:
        stage_dir = run_dir / report.STAGE_DIRS[stage]
        if stage == "qc":
            payload = json.loads((stage_dir / "metrics.json").read_text(encoding="utf-8"))
        else:
            payload = dict(metrics[stage])
            payload.update(
                {
                    "status": "pass",
                    "summary": f"{report._STAGE_META[stage]['name']}已完成。",
                }
            )
        summary_override = "自定义阶段概要（自测注入文本）。" if stage == "popgen" else None
        report_path = expect_raises_free(
            f"stage_report({stage}) executes",
            report.stage_report,
            run_dir,
            stage,
            payload,
            CONFIG,
            summary_override,
        )
        if report_path is None:
            continue
        check(f"stage_report({stage}) returns Path", isinstance(report_path, Path), str(report_path))
        exists = report_path.exists()
        size = report_path.stat().st_size if exists else 0
        check(f"stage_report({stage}) report.md exists", exists)
        check(f"stage_report({stage}) report.md > 1KB", size > 1024, f"{size / 1024:.1f} KB")
        text = report_path.read_text(encoding="utf-8") if exists else ""
        for section in expected_sections:
            check(f"stage_report({stage}) contains '{section}'", section in text)
        if stage == "intake":
            check("stage_report(intake) renders preflight checks", "files_exist" in text)
        if stage == "qc":
            check("stage_report(qc) renders WARN gate", "WARN" in text and "sample_missingness" in text)
            check("stage_report(qc) renders table", "sample_qc" in text or "IID" in text)
        if stage == "popgen":
            check("stage_report(popgen) honours summary_override", "自定义阶段概要" in text)
            check("stage_report(popgen) references figure", "plots/" in text)
        if stage == "gwas":
            check("stage_report(gwas) renders lead SNP detail", "Lead SNP" in text and "rs100000" in text)
            check(
                "stage_report(gwas) prefers p_text over underflowed p=0",
                "<1e-300" in text and "p = 0" not in text and "P = 0" not in text,
            )
            check("stage_report(gwas) renders P value distribution", "最小 P 值" in text and "<1e-300" in text)
        if stage == "calibration":
            check("stage_report(calibration) renders validation", "独立重算" in text and "校验" in text)


def _test_final_report(run_dir: Path, state: dict) -> None:
    print("\n--- report.final_report ---", flush=True)
    result = expect_raises_free("final_report executes", report.final_report, run_dir, CONFIG, state)
    if result is None:
        return
    check("final_report returns dict", isinstance(result, dict))
    for key in ("md", "html", "summary"):
        check(f"final_report result has '{key}'", key in result and isinstance(result[key], str))
    if not all(key in result for key in ("md", "html", "summary")):
        return
    md_path = Path(result["md"])
    html_path = Path(result["html"])
    summary_path = Path(result["summary"])
    for label, path in (("md", md_path), ("html", html_path), ("summary", summary_path)):
        check(f"final {label} exists", path.exists(), str(path))
    md_size = md_path.stat().st_size if md_path.exists() else 0
    html_size = html_path.stat().st_size if html_path.exists() else 0
    check("final_report.md > 8KB", md_size > 8 * 1024, f"{md_size / 1024:.1f} KB")
    check("final_report.html > 100KB (base64 figures)", html_size > 100 * 1024, f"{html_size / 1024:.1f} KB")
    md_text = md_path.read_text(encoding="utf-8") if md_path.exists() else ""
    html_text = html_path.read_text(encoding="utf-8") if html_path.exists() else ""
    for phrase in ("执行摘要", "质量控制", "Lead SNP", "复现命令", "局限性", "文件与产物清单", "λGC"):
        check(f"final_report.md contains '{phrase}'", phrase in md_text)
    for phrase in ("data:image/png;base64,", "<style>", "</html>", "st-pass", "st-warn", "table-wrap"):
        check(f"final_report.html contains '{phrase}'", phrase in html_text)
    check("final_report.html embeds all figures (no missing-image marker)", "图片缺失" not in html_text)
    check(
        "final_report uses p_text for underflowed p values",
        "<1e-300" in md_text and "1e-300" in html_text,
    )
    check(
        "final_report never prints p = 0",
        "p = 0" not in md_text and "P = 0" not in md_text and "p = 0" not in html_text,
    )
    check("final_report includes P value distribution", "P 值分布" in md_text and "最小 P 值" in md_text)
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception as exc:
        check("summary.json parses", False, str(exc))
        return
    check("summary.json parses", True)
    check("summary.json has run_id", summary.get("run_id") == CONFIG["run_id"])
    check("summary.json overall_status == 'warn'", summary.get("overall_status") == "warn", str(summary.get("overall_status")))
    check("summary.json data.n_samples == 400", summary.get("data", {}).get("n_samples") == 400)
    check(
        "summary.json gwas lambda_gc == 1.01",
        abs(float(summary.get("gwas", {}).get("lambda_gc", 0)) - 1.01) < 1e-9,
    )
    check("summary.json figures embedded", isinstance(summary.get("figures"), list) and len(summary["figures"]) >= 10)
    check("summary.json calibration null recorded", summary.get("calibration", {}).get("null", {}).get("lambda_mean") == 1.0)


def _test_skipped_and_empty(root: Path) -> None:
    print("\n--- skipped stage / empty run ---", flush=True)
    skipped_run = root / "skipped_run"
    stage_dir = skipped_run / report.STAGE_DIRS["gwas"]
    stage_dir.mkdir(parents=True, exist_ok=True)
    skipped_payload = {"stage": "gwas", "status": "skipped", "summary": "本次运行跳过了 GWAS 阶段。"}
    report_path = expect_raises_free(
        "stage_report(gwas, skipped) executes",
        report.stage_report,
        skipped_run,
        "gwas",
        skipped_payload,
        CONFIG,
    )
    if report_path is not None:
        text = report_path.read_text(encoding="utf-8")
        check("skipped stage report exists", report_path.exists())
        check("skipped stage report > 1KB", report_path.stat().st_size > 1024)
        check("skipped stage report marks not executed", "未执行" in text and "SKIPPED" in text)

    # calibration skipped: final report must say "本次未执行（原因：…）" and not 0/0
    cal_dir = skipped_run / report.STAGE_DIRS["calibration"]
    cal_dir.mkdir(parents=True, exist_ok=True)
    (cal_dir / "metrics.json").write_text(
        json.dumps(
            {
                "stage": "calibration",
                "status": "skipped",
                "summary": "本机缺少 simulate.py，仿真校准未执行。",
                "null": {},
                "validation": {"n_checks": 0, "n_pass": 0, "n_warn": 0, "n_fail": 0},
            },
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    result = expect_raises_free(
        "final_report on run with skipped calibration", report.final_report, skipped_run, {"run_id": "skipped_run"}, None
    )
    if result is not None:
        for key in ("md", "html", "summary"):
            path = Path(result[key])
            check(f"skipped-run final {key} exists", path.exists(), str(path))
        md_text = Path(result["md"]).read_text(encoding="utf-8")
        check("skipped calibration explicitly says 本次未执行", "本次未执行" in md_text)
        check("skipped calibration gives a reason", "原因" in md_text and "simulate.py" in md_text)
        check("skipped calibration does not print 0/0", "通过 0/0" not in md_text and "0/0 项" not in md_text)
        check("skipped calibration is flagged as unverified", "未经独立验证" in md_text)
        check("limitations mention missing simulation calibration", "未完成仿真校准" in md_text)
        try:
            summary_json = json.loads(Path(result["summary"]).read_text(encoding="utf-8"))
        except Exception as exc:
            check("skipped-run summary.json parses", False, str(exc))
            summary_json = {}
        check(
            "summary.json records calibration skip reason",
            "simulate.py" in str(summary_json.get("calibration", {}).get("skip_reason")),
            str(summary_json.get("calibration", {}).get("skip_reason")),
        )
        check(
            "summary.json records calibration evidence flags",
            summary_json.get("calibration", {}).get("evidence", {}).get("null") is False,
        )

    # calibration warn with zero evidence: must not read as "通过 0/0"
    noevid_run = root / "noevid_run"
    noevid_cal = noevid_run / report.STAGE_DIRS["calibration"]
    noevid_cal.mkdir(parents=True, exist_ok=True)
    (noevid_cal / "metrics.json").write_text(
        json.dumps(
            {
                "stage": "calibration",
                "status": "warn",
                "summary": "零模型校准：未获得结果。",
                "null": {},
                "validation": {"n_checks": 0, "n_pass": 0, "n_warn": 0, "n_fail": 0},
            },
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    result = expect_raises_free(
        "final_report on run with warn calibration but no evidence",
        report.final_report,
        noevid_run,
        {"run_id": "noevid_run"},
        None,
    )
    if result is not None:
        md_text = Path(result["md"]).read_text(encoding="utf-8")
        check("no-evidence calibration reports no valid results", "未获得有效" in md_text)
        check("no-evidence calibration avoids 0/0 wording", "通过 0/0" not in md_text)
        check("no-evidence calibration notes unverified p values", "未经独立验证" in md_text)

    # completely empty run directory: report is still generated
    empty_run = root / "completely_empty_run"
    empty_run.mkdir(parents=True, exist_ok=True)
    result = expect_raises_free(
        "final_report on empty run directory", report.final_report, empty_run, {"run_id": "empty_run"}, None
    )
    if result is not None:
        for key in ("md", "html", "summary"):
            path = Path(result[key])
            check(f"empty-run final {key} exists", path.exists(), str(path))
        md_text = Path(result["md"]).read_text(encoding="utf-8")
        check("empty-run final report notes missing stages", "未执行" in md_text or "未产生" in md_text)
        check("empty-run final HTML is valid shell", "</html>" in Path(result["html"]).read_text(encoding="utf-8"))


def _test_p_value_helpers() -> None:
    print("\n--- P-value display helpers ---", flush=True)
    check("_p_value_text(0) -> '<1e-300'", report._p_value_text(0.0) == "<1e-300")
    check("_p_value_text(0.0) never returns '0'", report._p_value_text(0.0) != "0")
    check(
        "_p_display falls back to numeric p and marks underflow",
        report._p_display({"p": 0.0}) == "<1e-300",
    )
    check(
        "_p_display prefers p_text when present",
        report._p_display({"p": 0.0, "p_text": "<1e-300"}) == "<1e-300",
    )
    check(
        "_p_display keeps precise p_text",
        report._p_display({"p": 0.0, "p_text": "1.273e-12"}) == "1.273e-12",
    )
    check(
        "_min_p_display falls back to numeric min with underflow marker",
        report._min_p_display({"min": 0.0}) == "<1e-300",
    )
    headers, rows = report._collapse_p_text_columns(["SNP", "P"], [["rs1", "0"], ["rs2", "1.0e-08"]])
    check(
        "CSV numeric P=0 displayed as '<1e-300'",
        rows[0][1] == "<1e-300" and rows[1][1] == "1.0e-08",
        str(rows),
    )
    headers, rows = report._collapse_p_text_columns(
        ["SNP", "P", "P_TEXT"], [["rs1", "0", "<1e-300"], ["rs2", "1e-08", "1.000e-08"]]
    )
    check(
        "CSV P_TEXT column replaces numeric P column",
        headers == ["SNP", "P 值"] and rows[0][1] == "<1e-300" and rows[1][1] == "1.000e-08",
        f"{headers} {rows}",
    )


def _test_renderer(root: Path) -> None:
    print("\n--- report.render_markdown_to_html ---", flush=True)
    base = root / "renderer"
    figures_dir = base / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    plots.plot_maf_spectrum(pd.Series(np.clip(np.random.default_rng(3).beta(1.2, 4, 2000), 0.001, 0.5)), figures_dir / "test.png")
    markdown = """# 标题一

这是一个**粗体**、*斜体*、`行内代码` 与 [链接](https://example.com) 的段落。

| 项目 | 结果 |
| :--- | ---: |
| A | ✅ PASS |
| B | ⚠️ WARN |

- 列表项一
- 列表项二

> 引用文本。

```python
print("hello")
```

![测试图](figures/test.png)

![缺失图](figures/not_found.png)
"""
    out = base / "test.html"
    result = expect_raises_free("render_markdown_to_html executes", report.render_markdown_to_html, markdown, out, title="自测标题")
    if result is None:
        return
    check("renderer returns Path", isinstance(result, Path))
    check("renderer output exists", out.exists())
    html = out.read_text(encoding="utf-8") if out.exists() else ""
    for fragment in (
        "<h1>标题一</h1>",
        "<strong>粗体</strong>",
        "<em>斜体</em>",
        "<code>行内代码</code>",
        '<a href="https://example.com">链接</a>',
        "<table>",
        "st-pass",
        "st-warn",
        "<ul>",
        "<blockquote>",
        '<pre><code class="language-python">',
        "data:image/png;base64,",
        "<figcaption>测试图</figcaption>",
        "图片缺失",
        "<!DOCTYPE html>",
        "<title>自测标题</title>",
    ):
        check(f"renderer output contains '{fragment[:32]}'", fragment in html)


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    root = Path(tempfile.mkdtemp(prefix="plink_report_selftest_"))
    print(f"selftest sandbox: {root}", flush=True)
    try:
        run_dir, _, state, metrics = _make_run(root)
        _test_plots(root)
        _test_stage_reports(run_dir, state, metrics)
        _test_final_report(run_dir, state)
        _test_skipped_and_empty(root)
        _test_p_value_helpers()
        _test_renderer(root)
    except Exception:  # pragma: no cover - unexpected crash
        traceback.print_exc()
        CHECKS.append(("selftest completed without unexpected exception", False, "crashed"))
    finally:
        shutil.rmtree(root, ignore_errors=True)
        print(f"\nselftest sandbox removed: {root}", flush=True)

    passed = sum(1 for _name, ok, _detail in CHECKS if ok)
    failed = [(name, detail) for name, ok, detail in CHECKS if not ok]
    print("\n" + "=" * 72)
    print(f"CHECKS: {passed}/{len(CHECKS)} passed")
    if failed:
        print("FAILED CHECKS:")
        for name, detail in failed:
            print(f"  - {name}" + (f"  [{detail}]" if detail else ""))
        print("RESULT: FAIL")
        return 1
    print("RESULT: ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
