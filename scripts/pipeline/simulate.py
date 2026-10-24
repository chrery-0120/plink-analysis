#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""仿真校准：零模型（Type-I error / λGC）+ 已知 QTL（power / 效应量 / 排序一致性）。

为什么需要
----------
GWAS 的 p 值只有在"零假设下均匀分布"时才可信。本模块用**与主 GWAS 完全相同的
PLINK 命令与协变量**（``--linear/--logistic --covar ... --keep-allele-order``）在
两类仿真数据上重复 R 次，直接度量：

1. **零模型**：表型与基因型无关（随机正态 / 随机 1/2 病例标签），
   统计每次的 λGC 与 α=0.05 / 1e-3 / 5e-8 的经验假阳性率（FPR）。
   验收：λ 均值 ∈ [0.95, 1.05]，α=0.05 的 FPR ∈ [0.02, 0.08]，α=5e-8 的 FPR≈0。
2. **已知 QTL**：用 truth 中的因果 SNP 按 h²=0.3 / 0.5 生成表型（因果效应 +
   多基因背景 + 环境噪声），统计全基因组显著（5e-8）与 1e-4 阈值的检出率、
   效应量偏差（β_est/β_true）、SE 标定（SE_est vs β_est 跨重复的 SD）、
   因果 SNP 排序 AUC 与 top100 命中率。
3. **可复现**：种子来自 ``config.calibration.seed``（默认 20261023）。

调用约定
--------
``run_calibration(ctx_or_run_dir) -> dict``

* calibration_runner 的调用方式（推荐）::

      from .simulate import run_calibration
      metrics = run_calibration(ctx)     # ctx: PipelineContext

* 单独对已有 run 目录执行（自动读 run.json/state.json、定位 PLINK/truth）::

      python -m pipeline.simulate D:/plink_runs/val_qtl

产物（``<run_dir>/04_calibration/``）::

    metrics.json                        # 结构化结果（calibration_runner 会合并）
    tables/calibration_null.csv         # 每次重复的 λ 与 FPR
    tables/calibration_power.csv        # 每次重复的 power
    tables/calibration_effect_bias.csv  # 因果 SNP 的真实/估计效应
    tables/calibration_concordance.csv  # 排序一致性
    plots/calibration_null_lambda.png
    plots/calibration_power.png
    plots/calibration_effect_bias.png

容差与验收阈值在 ``GATES`` 中集中定义，可用 ``run.json`` 的
``calibration.gates`` 覆盖。
"""
from __future__ import annotations

import datetime as _dt
import logging
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:  # pipeline 包内
    from . import utils
    from .validate import (
        _detect_header, _prefix_of, _read_plink_table, _recode, _resolve,
        discover_artifacts,
    )
except ImportError:  # pragma: no cover - 兼容直接执行
    _here = Path(__file__).resolve().parent
    sys.path.insert(0, str(_here.parent))
    sys.path.insert(0, str(_here))
    import utils  # type: ignore
    from validate import (  # type: ignore
        _detect_header, _prefix_of, _read_plink_table, _recode, _resolve,
        discover_artifacts,
    )

# ---------------------------------------------------------------------------
# 常量 / 验收阈值
# ---------------------------------------------------------------------------
ALPHAS = (0.05, 1e-3, 1e-4, 5e-8)
ALPHA_KEYS = {0.05: "alpha_0.05", 1e-3: "alpha_1e-3", 1e-4: "alpha_1e-4", 5e-8: "alpha_5e-8"}
ALPHA_GW = 5e-8
DEFAULT_N_CAUSAL = 20
DEFAULT_N_BG = 1500
DEFAULT_H2 = (0.3, 0.5)
DEFAULT_SEED = 20261023
DEFAULT_MINUTES = 9.0            # 校准总时长预算（selftest 15 分钟里留余量）

GATES: dict[str, Any] = {
    "lambda_pass": (0.95, 1.05),     # 零模型 λ 均值验收区间（Lead/methods 共识）
    "lambda_warn": (0.90, 1.10),
    "fpr_alpha_0.05": (0.02, 0.08),
    "fpr_warn": (0.01, 0.12),
    "rank_auc_min": 0.65,
    "top100_min": 0.15,
    "max_warn": 6,                   # 允许的 warn 项数上限（超过则整体 warn→fail 由调用方决定）
}


# ---------------------------------------------------------------------------
# 数据准备
# ---------------------------------------------------------------------------
def _read_pheno_any(path: Path) -> pd.DataFrame:
    """读 FID/IID/PHENO（自动识别表头），保留 -9/NA 用于判断可用样本。"""
    has_header = _detect_header(path)
    df = pd.read_csv(path, sep=r"\s+", header=0 if has_header else None,
                     encoding="utf-8-sig")
    if has_header:
        df.columns = [str(c).lstrip("#").strip() for c in df.columns]
        df = df.iloc[:, :3].copy()
    else:
        df = df.iloc[:, :3].copy()
    df.columns = ["FID", "IID", "PHENO"]
    df["FID"], df["IID"] = df["FID"].astype(str), df["IID"].astype(str)
    df["PHENO_RAW"] = df["PHENO"]
    df["PHENO"] = pd.to_numeric(df["PHENO"], errors="coerce")
    return df


def _sample_table(t) -> pd.DataFrame:
    """clean.fam 的 FID/IID（分析样本的基础集合）。"""
    clean = _prefix_of(t.artifacts.get("clean"))
    fam = utils.read_fam(Path(str(clean) + ".fam"))
    fam["FID"], fam["IID"] = fam["FID"].astype(str), fam["IID"].astype(str)
    return fam[["FID", "IID"]].copy()


def _covar_info(t, found) -> tuple[Path | None, list[str]]:
    covar_p = found["covar"]["path"]
    if covar_p is None or not Path(covar_p).exists():
        return None, []
    cols = list((t.metrics.get("gwas") or {}).get("covariates") or [])
    if not cols:
        try:
            first = Path(covar_p).read_text(encoding="utf-8-sig",
                                            errors="replace").splitlines()[0].split()
            cols = [f"COV{i + 1}" for i in range(max(0, len(first) - 2))]
        except Exception:
            cols = []
    return Path(covar_p), cols


def _find_truth(t, found) -> tuple[Path | None, str]:
    """定位 truth.csv；顺序：config.calibration.truth_path → 表型/基因型同目录 → run 目录。"""
    candidates: list[tuple[Path, str]] = []
    cfg = t.config or {}
    tp = (cfg.get("calibration") or {}).get("truth_path")
    if tp:
        candidates.append((Path(str(tp)), "calibration.truth_path"))
    for key, src in (
        (((t.metrics.get("intake") or {}).get("phenotype") or {}).get("file"),
         "intake.metrics.phenotype.file"),
        (((cfg.get("input") or {}).get("phenotype") or {}).get("path"),
         "config.input.phenotype.path"),
        ((((cfg.get("input") or {}).get("genotype") or {}).get("prefix")),
         "config.input.genotype.prefix"),
        (str(t.artifacts.get("base")), "00_intake base_prefix"),
        (str(found.get("clean", {}).get("path")), "01_qc clean"),
    ):
        if not key:
            continue
        try:
            candidates.append((Path(str(key)).parent / "truth.csv", src + "/../truth.csv"))
        except Exception:
            continue
    candidates += [
        (t.run_dir / "truth.csv", "run_dir/truth.csv"),
        (t.run_dir.parent / "truth.csv", "run_dir/../truth.csv"),
    ]
    for path, src in candidates:
        try:
            if path.is_file():
                head = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()[0]
                if "SNP" in head:
                    return path, src
        except Exception:
            continue
    return None, "未找到 truth.csv"


def _load_truth(path: Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    try:
        truth = pd.read_csv(path, encoding="utf-8-sig")
    except Exception:
        return None
    truth.columns = [str(c).strip() for c in truth.columns]
    if "SNP" not in truth.columns:
        return None
    for col in ("BETA_UNIT", "BETA_QT", "BETA_LIAB", "AF_A", "AF_B", "MAF"):
        if col in truth.columns:
            truth[col] = pd.to_numeric(truth[col], errors="coerce")
    return truth


def _code_truth_flag(series: pd.Series) -> np.ndarray:
    """把 IS_CAUSAL 之类的列规范成 bool 数组。"""
    vals = series.astype(str).str.strip().str.lower()
    return vals.isin({"1", "true", "t", "yes", "y", "causal"}).to_numpy()


def _pick_snp_set(t, found, n_subset: int, seed: int) -> tuple[list[str], str]:
    """选择校准用的 SNP 集合：优先 LD 剪枝集（控制 LD 冗余），过大则确定性抽样。"""
    clean = _prefix_of(t.artifacts.get("clean"))
    bim = utils.read_bim(Path(str(clean) + ".bim"))
    all_snps = bim["SNP"].astype(str).tolist()
    prune_p = found.get("pruned", {}).get("path")
    if prune_p is not None and Path(prune_p).exists():
        snps = [s for s in Path(prune_p).read_text(encoding="utf-8-sig",
                                                   errors="replace").split() if s]
        src = "LD 剪枝位点"
        snps = [s for s in snps if s in set(all_snps)]
    else:
        snps, src = all_snps, "全部 QC 后位点"
    if len(snps) > n_subset:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(len(snps), n_subset, replace=False))
        snps = [snps[i] for i in idx]
        src += f"（确定性抽取 {n_subset} 个）"
    return snps, src


def _prepare_effects(t, truth: pd.DataFrame | None, snps: list[str], raw,
                     n_causal: int, n_bg: int, seed: int,
                     phenotype_type: str) -> dict[str, Any]:
    """确定因果位点与单位遗传方差尺度上的真实效应量（β_unit）。

    注意：truth.csv 的 β 是相对 truth 中 A1 的效应；若当前 bfile 的 A1 被翻转，
    这里按 ``raw.alleles``（= 当前 .bim A1）翻转符号。
    """
    rng = np.random.default_rng(seed + 7)
    snp_index = {s: j for j, s in enumerate(raw.snps)}
    x = raw.matrix()
    p = np.clip(np.nanmean(x, axis=0) / 2.0, 0.01, 0.99)
    maf = np.minimum(p, 1 - p)
    used = [s for s in snps if s in snp_index]
    if not used:
        raise ValueError("校准 SNP 集合与 .raw 无交集")
    px = np.array([p[snp_index[s]] for s in used])
    maf_used = np.array([maf[snp_index[s]] for s in used])
    info: dict[str, Any] = {"source": "random", "truth_path": None, "n_bg": 0}

    causal: list[str] = []
    beta_c = np.zeros(0)
    bg: list[str] = []
    beta_bg = np.zeros(0)
    if truth is not None:
        tmap = truth.set_index(truth["SNP"].astype(str))
        used_set = set(used)
        allele_in_truth = (dict(zip(truth["SNP"].astype(str), truth["A1"].astype(str)))
                           if "A1" in truth.columns else {})

        def _beta_of(snp: str) -> float:
            val = float(tmap.loc[snp, "BETA_UNIT"])
            if allele_in_truth and allele_in_truth.get(snp) not in (None, "") \
                    and raw.alleles.get(snp) not in (None, "") \
                    and allele_in_truth[snp] != raw.alleles[snp]:
                val = -val
            return val

        if "IS_CAUSAL" in truth.columns and "BETA_UNIT" in truth.columns:
            is_c = _code_truth_flag(truth["IS_CAUSAL"])
            causal = [s for s, flag in zip(truth["SNP"].astype(str), is_c)
                      if flag and s in used_set]
        if causal:
            causal = sorted(causal, key=lambda s: -abs(_beta_of(s)))
            causal = causal[:n_causal]
            beta_c = np.array([_beta_of(s) for s in causal])
            info.update(source="truth.csv", n_truth_causal=len(causal))
            if "IS_POLY" in truth.columns:
                is_bg = _code_truth_flag(truth["IS_POLY"])
                bg = [s for s, flag in zip(truth["SNP"].astype(str), is_bg)
                      if flag and s in used_set and s not in set(causal)]
                if bg:
                    beta_bg = np.array([_beta_of(s) for s in bg])
                    info["n_bg"] = len(bg)
    if not causal:
        # 无 truth：随机选高 MAF 位点当因果位点，效应 ~ N(0,1) 后归一化
        cand = [s for s in used if maf[snp_index[s]] >= 0.05]
        if len(cand) < n_causal:
            cand = used
        picks = rng.choice(np.array(cand), min(n_causal, len(cand)), replace=False)
        causal = [str(s) for s in picks]
        b = rng.normal(size=len(causal))
        q = np.array([px[snp_index[s]] for s in causal])
        b *= np.sqrt(1.0 / max(np.sum(2 * q * (1 - q) * b ** 2), 1e-12))
        beta_c = b
        info["source"] = "random_effects"
        rest = [s for s in used if s not in set(causal)]
        if rest and n_bg > 0:
            bg = [str(s) for s in rng.choice(np.array(rest), min(n_bg, len(rest)),
                                             replace=False)]
            bb = rng.normal(size=len(bg))
            q = np.array([px[snp_index[s]] for s in bg])
            denom = np.sqrt(max(np.sum(2 * q * (1 - q) * bb ** 2), 1e-12))
            beta_bg = bb * np.sqrt(0.15) / denom
            info["n_bg"] = len(bg)
    info.update(causal=causal, beta_c=beta_c, bg=bg, beta_bg=beta_bg,
                snp_px={s: float(px[snp_index[s]]) for s in used},
                snp_maf={s: float(maf_used[i]) for i, s in enumerate(used)},
                phenotype_type=phenotype_type)
    return info


def _genetic_and_true_beta(raw, eff: dict[str, Any], rows: np.ndarray,
                           h2: float) -> tuple[np.ndarray, dict[str, float], float]:
    """在指定样本行上构造标准化的单位遗传值 G_std 与真实每等位基因效应 β_true。"""
    cols = {s: j for j, s in enumerate(raw.snps)}
    x = raw.matrix()[rows, :]

    def g_of(names: list[str], betas: np.ndarray) -> np.ndarray:
        if not names:
            return np.zeros(x.shape[0])
        keep = [(i, s) for i, s in enumerate(names) if s in cols]
        if not keep:
            return np.zeros(x.shape[0])
        idx = [cols[s] for _, s in keep]
        b = np.array([betas[i] for i, _ in keep], dtype=float)
        xx = x[:, idx]
        pj = np.clip(np.nanmean(xx, axis=0) / 2.0, 1e-6, 1 - 1e-6)
        xx = np.where(np.isnan(xx), 2.0 * pj[None, :], xx)
        return xx @ b

    g = g_of(eff["causal"], eff["beta_c"]) + g_of(eff["bg"], eff["beta_bg"])
    sd = float(np.std(g, ddof=0))
    if sd <= 0:
        raise ValueError("仿真遗传值方差为 0（因果位点可能全是单态）")
    g_std = (g - g.mean()) / sd
    a = math.sqrt(h2 / (1.0 - h2))
    beta_true = {s: float(eff["beta_c"][i] * a / sd)
                 for i, s in enumerate(eff["causal"])}
    return a * g_std, beta_true, a


# ---------------------------------------------------------------------------
# PLINK 关联运行
# ---------------------------------------------------------------------------
def _assoc_parse(t, prefix: Path, model: str) -> pd.DataFrame:
    """解析本轮关联结果（只保留 TEST=ADD），统一成 CHR/SNP/BP/A1/BETA/SE/P。"""
    for suffix in (".assoc.linear", ".assoc.logistic"):
        path = Path(str(prefix) + suffix)
        if path.exists():
            break
    else:
        raise FileNotFoundError(f"缺少关联结果: {prefix}.assoc.(linear|logistic)")
    df = _read_plink_table(path)
    if "TEST" in df.columns:
        df = df.loc[df["TEST"].astype(str).str.upper() == "ADD"].copy()
    for col in ("BETA", "OR", "SE", "L95", "U95", "P", "NMISS", "STAT"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "P" not in df.columns:
        raise ValueError(f"关联结果缺少 P 列: {path}")
    if "BETA" in df.columns:
        df["BETA_PLINK"] = df["BETA"]
    elif "OR" in df.columns:
        orv = df["OR"].astype(float)
        df["BETA_PLINK"] = np.where(orv > 0, np.log(orv), np.nan)
    else:
        df["BETA_PLINK"] = np.nan
    if "SE" not in df.columns and {"L95", "U95"}.issubset(df.columns):
        lo, hi = df["L95"].astype(float), df["U95"].astype(float)
        if model == "logistic":
            df["SE"] = (np.log(hi.where(hi > 0)) - np.log(lo.where(lo > 0))) / (2 * 1.959964)
        else:
            df["SE"] = (hi - lo) / (2 * 1.959964)
    if "SE" not in df.columns:
        df["SE"] = np.nan
    return df[["CHR", "SNP", "BP", "A1", "BETA_PLINK", "SE", "P", "NMISS"]].copy()


def _run_assoc(t, prefix: Path, keep_p: Path, pheno_p: Path, covar_p: Path | None,
               covar_cols: list[str], model: str, name: str) -> pd.DataFrame:
    clean = _prefix_of(t.artifacts.get("clean"))
    args = ["--bfile", str(clean), "--keep", str(keep_p), "--pheno", str(pheno_p)]
    if covar_p is not None and covar_cols:
        args += ["--covar", str(covar_p), "--covar-number", f"1-{len(covar_cols)}"]
    args += ["--keep-allele-order", "--allow-no-sex", "--ci", "0.95",
             "--linear" if model == "linear" else "--logistic"]
    pk = t.plink
    if pk is None:
        raise FileNotFoundError("未找到 PLINK，无法运行校准仿真")
    pk.run(args, stage="calibration", name=name, out_prefix=prefix, check=True)
    return _assoc_parse(t, prefix, model)


# ---------------------------------------------------------------------------
# 绘图（优先用 pipeline.plots，失败则本地兜底）
# ---------------------------------------------------------------------------
def _plot(t, kind: str, data, out_path: Path) -> Path | None:
    try:
        from . import plots as plots_mod  # type: ignore
    except Exception:
        try:
            import plots as plots_mod  # type: ignore
        except Exception:
            plots_mod = None
    if plots_mod is not None:
        try:
            if kind == "null":
                return Path(plots_mod.plot_calibration(data, out_path))
            if kind == "power":
                return Path(plots_mod.plot_power_curve(data, out_path))
            if kind == "effect":
                return Path(plots_mod.plot_effect_bias(data, out_path))
        except Exception as exc:  # noqa: BLE001
            t.logger.warning("plots.%s 失败，改用内置绘图: %s", kind, exc)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=150)
    try:
        if kind == "null":
            arr = np.asarray(list(data), dtype=float)
            ax.hist(arr, bins=max(5, min(20, len(arr))), color="#4C72B0", alpha=0.85)
            ax.axvline(1.0, color="k", ls="--", lw=1)
            ax.set_xlabel("lambda GC")
            ax.set_ylabel("replicates")
            ax.set_title("Null-model lambda GC")
        elif kind == "power":
            df = pd.DataFrame(data)
            ax.plot(df["alpha"], df["power"], marker="o", color="#C44E52")
            ax.set_xscale("log")
            ax.set_xlabel("alpha")
            ax.set_ylabel("empirical power")
            ax.set_title("Known-QTL power")
        elif kind == "effect":
            df = pd.DataFrame(data)
            ax.scatter(df["true_beta"], df["est_beta"], s=14, alpha=0.7)
            lim = [min(df["true_beta"].min(), df["est_beta"].min()),
                   max(df["true_beta"].max(), df["est_beta"].max())]
            ax.plot(lim, lim, "k--", lw=1)
            ax.set_xlabel("true beta (/allele)")
            ax.set_ylabel("estimated beta")
            ax.set_title("Effect-size calibration")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(out_path)
        return out_path
    finally:
        plt.close(fig)


# ---------------------------------------------------------------------------
# 统计
# ---------------------------------------------------------------------------
def _auc(causal_p: np.ndarray, other_p: np.ndarray) -> float:
    """因果 SNP 的 p 值比随机 SNP 更小的概率（Mann-Whitney AUC，越大越好）。"""
    if causal_p.size == 0 or other_p.size == 0:
        return float("nan")
    # 用 -log10(p) 的秩和计算 AUC
    x = -np.log10(np.clip(causal_p, 1e-300, 1))
    y = -np.log10(np.clip(other_p, 1e-300, 1))
    greater = (x[:, None] > y[None, :]).sum()
    ties = (x[:, None] == y[None, :]).sum()
    return float((greater + 0.5 * ties) / (x.size * y.size))


def _summ(vals: list[float]) -> dict[str, float]:
    arr = np.asarray([v for v in vals if v is not None and np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return {"mean": float("nan"), "sd": float("nan"), "min": float("nan"),
                "max": float("nan"), "n": 0}
    return {"mean": float(arr.mean()), "sd": float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
            "min": float(arr.min()), "max": float(arr.max()), "n": int(arr.size)}


def _write_pheno(path: Path, fid: np.ndarray, iid: np.ndarray, values: np.ndarray) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("FID\tIID\tPHENO\n")
        for f, i, v in zip(fid, iid, values):
            fh.write(f"{f}\t{i}\t{float(v):.8f}\n")


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def run_calibration(ctx_or_run_dir: Any, *, logger: logging.Logger | None = None) -> dict[str, Any]:
    """执行仿真校准并写 metrics/tables/plots，返回 metrics dict（见模块 docstring）。"""
    started = time.time()
    logger = logger or logging.getLogger("calibration")
    t = _resolve(ctx_or_run_dir, logger)
    found = discover_artifacts(t)
    cfg = t.config or {}
    cal_cfg = dict((cfg.get("calibration") or {}))
    seed = int(cal_cfg.get("seed", DEFAULT_SEED) or DEFAULT_SEED)
    n_reps = max(3, int(cal_cfg.get("n_replicates", 10) or 10))
    n_causal = max(5, int(cal_cfg.get("n_causal", DEFAULT_N_CAUSAL) or DEFAULT_N_CAUSAL))
    n_subset = max(500, int(cal_cfg.get("n_snps_subset", 20000) or 20000))
    h2_primary = float(cal_cfg.get("h2", 0.5) or 0.5)
    h2_list = [h for h in dict.fromkeys([0.3, h2_primary]) if 0.05 <= h <= 0.95]
    budget_s = float(cal_cfg.get("max_minutes", DEFAULT_MINUTES) or DEFAULT_MINUTES) * 60.0
    gates = dict(GATES)
    for k, v in (cal_cfg.get("gates") or {}).items():
        if k in gates:
            gates[k] = v

    warnings: list[str] = []
    errors: list[str] = []
    pheno_p = found.get("pheno", {}).get("path")
    if pheno_p is None:
        metrics = _finish(t, started, {
            "stage": "calibration", "status": "skipped",
            "summary": "未提供表型文件，跳过仿真校准。",
            "null": None, "power": None, "effect_bias": None, "concordance": None,
            "warnings": ["缺少 00_intake/data/pheno_clean.txt"], "errors": [],
        })
        return metrics
    clean = _prefix_of(t.artifacts.get("clean"))
    pheno_type = ((cfg.get("_derived") or {}).get("phenotype_type")
                  or (t.metrics.get("gwas") or {}).get("phenotype_type")
                  or "quantitative")
    model = str(cal_cfg.get("model") or "auto")
    if model == "auto":
        model = "logistic" if pheno_type == "binary" else "linear"
    if pheno_type == "binary" and model != "logistic":
        warnings.append("表型为二元但 calibration.model=linear，已按 linear 运行（结果解释需谨慎）")
    logger.info("仿真校准开始：model=%s, h2=%s, R=%d, seed=%d", model, h2_list, n_reps, seed)

    # --- 样本 / 协变量 / 表型 ------------------------------------------------
    # 保持 clean.fam 的行序（.raw 的样本行序与 .fam 一致），用布尔掩码取子集，
    # 这样生成的表型 y 与 .raw 剂量矩阵的行严格一一对应。
    samples = _sample_table(t)
    pheno_all = _read_pheno_any(Path(pheno_p))
    covariates, covar_cols = _covar_info(t, found)
    pheno_usable = pheno_all.loc[pheno_all["PHENO"].notna() &
                                 (pheno_all["PHENO_RAW"].astype(str).str.strip() != "-9")]
    ph_keys = set(zip(pheno_usable["FID"].astype(str), pheno_usable["IID"].astype(str)))
    cv_keys: set[tuple[str, str]] = set()
    if covariates is not None and covar_cols:
        cov = pd.read_csv(covariates, sep=r"\s+", header=None, encoding="utf-8-sig",
                          dtype={0: str, 1: str})
        cv_keys = set(zip(cov[0].astype(str), cov[1].astype(str)))
    mask = np.array([
        ((f, i) in ph_keys) and ((not cv_keys) or ((f, i) in cv_keys))
        for f, i in zip(samples["FID"].astype(str), samples["IID"].astype(str))
    ])
    keep = samples.loc[mask].reset_index(drop=True)
    if len(keep) < 50:
        raise ValueError(f"可用于校准的样本只有 {len(keep)} 个（需要 >=50）")
    keep_p = t.data_dir / "calib_keep.txt"
    keep_p.write_text("".join(f"{r.FID}\t{r.IID}\n" for r in keep.itertuples()),
                      encoding="utf-8")
    fid = keep["FID"].to_numpy()
    iid = keep["IID"].to_numpy()
    n_samples = len(keep)
    sample_rows = np.flatnonzero(mask)       # .fam / .raw 中的行号

    # --- SNP 集合 / 剂量 / 因果位点 ------------------------------------------
    snps, snp_source = _pick_snp_set(t, found, n_subset, seed)
    raw = _recode(t, clean, snps, "calib")
    if len(raw.df) != len(samples):
        # 理论上 --recode A 的样本行 = clean.fam，防御性校验
        warnings.append(f".raw 样本数 {len(raw.df)} 与 clean.fam {len(samples)} 不一致，"
                        "校准将按 IID 重新对齐")
        row_lookup = {sid: k for k, sid in enumerate(raw.df["IID"].astype(str))}
        sample_rows = np.array([row_lookup[i] for i in keep["IID"].astype(str)])
    truth_path, truth_src = _find_truth(t, found)
    eff = _prepare_effects(t, _load_truth(truth_path), snps, raw,
                           n_causal, DEFAULT_N_BG, seed, pheno_type)
    eff["truth_path"] = str(truth_path) if truth_path else None
    if truth_path is None:
        warnings.append("未找到 truth.csv，已改用随机效应模拟因果位点（仍可校准 power/λ，"
                        "但\"已知 QTL\"不是数据集的真实因果位点）")
    if eff.get("n_bg"):
        logger.info("已知 QTL 校准：%d 个因果位点 + %d 个多基因背景位点（%s）",
                    len(eff["causal"]), eff["n_bg"], eff["source"])

    # --- 运行预算：先测一次单轮耗时再决定重复次数 -----------------------------
    def _run_round(prefix_name: str, values: np.ndarray, model_name: str) -> pd.DataFrame:
        pheno_file = t.data_dir / f"{prefix_name}.pheno"
        _write_pheno(pheno_file, fid, iid, values)
        return _run_assoc(t, t.data_dir / prefix_name, keep_p, pheno_file,
                          covariates, covar_cols, model_name, prefix_name)

    # 零模型第一轮
    rng = np.random.default_rng(seed)
    t_round0 = time.time()
    if model == "logistic":
        prev = float(np.mean(pheno_usable["PHENO"] == 2)) if len(pheno_usable) else 0.2
        prev = float(np.clip(prev, 0.05, 0.5))
        null_vals = np.where(rng.random(n_samples) < prev, 2, 1).astype(float)
    else:
        null_vals = rng.normal(size=n_samples)
    assoc0 = _run_round("calib_null_rep1", null_vals, model)
    single_s = max(time.time() - t_round0, 0.05)
    # 预留：null R + len(h2_list)*R 轮
    total_rounds_planned = n_reps * (1 + len(h2_list))
    budget_rounds = int(max(6, budget_s / single_s))
    reps_eff = n_reps
    if total_rounds_planned > budget_rounds:
        reps_eff = max(3, budget_rounds // (1 + len(h2_list)))
        warnings.append(f"时间预算 {budget_s / 60:.0f} 分钟：单轮约 {single_s:.2f}s，"
                        f"重复次数由 {n_reps} 降为 {reps_eff}（可用 calibration.max_minutes 调整）")
        logger.warning("校准重复次数降为 %d（预算 %.0fs）", reps_eff, budget_s)

    # --- 零模型 -------------------------------------------------------------
    null_rows: list[dict[str, Any]] = []
    type1 = {ALPHA_KEYS[a]: [] for a in ALPHAS}
    null_lambdas: list[float] = []
    for rep in range(reps_eff):
        if rep == 0:
            assoc = assoc0
        else:
            vals = (np.where(rng.random(n_samples) < prev, 2, 1).astype(float)
                    if model == "logistic" else rng.normal(size=n_samples))
            assoc = _run_round(f"calib_null_rep{rep + 1}", vals, model)
        p = assoc["P"].to_numpy(dtype=float)
        p = p[np.isfinite(p)]
        lam = utils.lambda_gc(p)
        null_lambdas.append(lam)
        row = {"replicate": rep + 1, "lambda_gc": round(lam, 6),
               "n_tests": int(p.size)}
        for alpha, akey in ALPHA_KEYS.items():
            frac = float(np.mean(p < alpha)) if p.size else float("nan")
            type1[akey].append(frac)
            row[f"fpr_{akey}"] = round(frac, 8)
        null_rows.append(row)
    null = {
        "n_replicates": len(null_rows),
        "lambda_mean": _summ(null_lambdas)["mean"],
        "lambda_sd": _summ(null_lambdas)["sd"],
        "lambda_min": _summ(null_lambdas)["min"],
        "lambda_max": _summ(null_lambdas)["max"],
        "type1": {k: _summ(v)["mean"] for k, v in type1.items()},
        "type1_sd": {k: _summ(v)["sd"] for k, v in type1.items()},
        "expected": {"alpha_0.05": 0.05, "alpha_1e-3": 1e-3,
                     "alpha_1e-4": 1e-4, "alpha_5e-8": ALPHA_GW},
        "definition": ("每次重复随机生成与基因型无关的表型，用主 GWAS 相同命令与协变量检测；"
                       "λGC 由 P 值中位数 χ² 计算，FPR = 该 α 下显著位点占检验位点的比例"),
    }
    lambda_pass = gates["lambda_pass"]

    # --- 已知 QTL -----------------------------------------------------------
    power_rows: list[dict[str, Any]] = []
    effect_rows: list[dict[str, Any]] = []
    concordance_rows: list[dict[str, Any]] = []
    power_by_h2: dict[str, Any] = {}
    power_frac: dict[float, list[float]] = {}
    effect_by_h2: dict[str, Any] = {}
    lambda_qtl: dict[str, float] = {}
    causal_set = set(eff["causal"])
    for h2 in h2_list:
        power_5e8: list[float] = []
        power_1e4: list[float] = []
        any_hit: list[float] = []
        aucs: list[float] = []
        top100: list[float] = []
        top1: list[float] = []
        lambdas: list[float] = []
        beta_pairs: list[dict[str, Any]] = []
        sign_obs: list[float] = []          # logistic：估计 log-OR 与 liability β 的符号一致率
        for rep in range(reps_eff):
            g_std, beta_true, _a = _genetic_and_true_beta(raw, eff, sample_rows, h2)
            eps = rng.normal(size=n_samples)
            if model == "logistic":
                liab = g_std + eps
                case_rate = prev if "prev" in dir() else 0.2
                thr = float(np.quantile(liab, 1 - float(case_rate)))
                y = np.where(liab > thr, 2, 1).astype(float)
            else:
                y = g_std + eps
            assoc = _run_round(f"calib_h2_{h2:g}_rep{rep + 1}", y, model)
            p = assoc["P"].to_numpy(dtype=float)
            lam = utils.lambda_gc(p[np.isfinite(p)])
            lambdas.append(lam)
            smap = assoc.set_index(assoc["SNP"].astype(str))
            causal_rows = smap.loc[[s for s in eff["causal"] if s in smap.index]]
            if isinstance(causal_rows, pd.Series):
                causal_rows = causal_rows.to_frame().T
            cp = causal_rows["P"].to_numpy(dtype=float)
            cp = cp[np.isfinite(cp)]
            fracs = {a: (float(np.mean(cp < a)) if cp.size else float("nan"))
                     for a in ALPHAS}
            frac5, frac1e4 = fracs[ALPHA_GW], fracs[1e-4]
            hit = float(np.any(cp < ALPHA_GW)) if cp.size else 0.0
            power_5e8.append(frac5)
            power_1e4.append(frac1e4)
            any_hit.append(hit)
            for a in ALPHAS:
                power_frac.setdefault(a, []).append(fracs[a])
            # 排序一致性
            other_mask = ~assoc["SNP"].astype(str).isin(causal_set)
            other_p = assoc.loc[other_mask, "P"].to_numpy(dtype=float)
            aucs.append(_auc(cp, other_p[np.isfinite(other_p)]))
            order = np.argsort(p)
            top_set = set(assoc["SNP"].astype(str).to_numpy()[order[:100]].tolist())
            top100.append(float(len(causal_set & top_set) / max(len(causal_set), 1)))
            top1.append(float(str(assoc["SNP"].astype(str).to_numpy()[order[0]]) in causal_set))
            concordance_rows.append({"h2": h2, "replicate": rep + 1,
                                     "rank_auc": round(aucs[-1], 5),
                                     "top_causal_in_top100": round(top100[-1], 5),
                                     "top1_is_causal": int(top1[-1]),
                                     "lambda_gc": round(lam, 5)})
            # 效应量对照：线性模型比 β 值；logistic 只比方向（log-OR 与 liability β
            # 的尺度不可直接换算，但符号应当一致）
            if model == "logistic":
                for snp in eff["causal"]:
                    if snp not in smap.index:
                        continue
                    row = smap.loc[snp]
                    if isinstance(row, pd.DataFrame):
                        row = row.iloc[0]
                    b_est = float(row["BETA_PLINK"])
                    b_true = beta_true.get(snp, float("nan"))
                    if np.isfinite(b_est) and np.isfinite(b_true) and abs(b_true) > 1e-8:
                        sign_obs.append(1.0 if (b_est > 0) == (b_true > 0) else 0.0)
            if model == "linear":
                for snp in eff["causal"]:
                    if snp not in smap.index:
                        continue
                    row = smap.loc[snp]
                    if isinstance(row, pd.DataFrame):
                        row = row.iloc[0]
                    b_hat = float(row["BETA_PLINK"])
                    se_hat = float(row["SE"])
                    if np.isfinite(b_hat):
                        beta_pairs.append({"h2": h2, "replicate": rep + 1, "snp": snp,
                                           "beta_true": beta_true[snp],
                                           "beta_est": b_hat, "se_est": se_hat})
            power_rows.append({"h2": h2, "replicate": rep + 1,
                               **{f"power_{a:g}": round(fracs[a], 6) for a in ALPHAS},
                               "any_hit_5e-8": int(hit),
                               "lambda_gc": round(lam, 6),
                               "n_causal_checked": int(cp.size)})
        lam_s = _summ(lambdas)
        lambda_qtl[f"{h2:g}"] = lam_s["mean"]
        p5, p1e4 = _summ(power_5e8), _summ(power_1e4)
        power_by_h2[f"{h2:g}"] = {
            "h2": h2, "n_causal": int(len(eff["causal"])),
            "n_replicates": len(power_rows) and reps_eff,
            "alpha_5e-8": {"mean": p5["mean"], "sd": p5["sd"],
                           "any_hit_rate": _summ(any_hit)["mean"]},
            "alpha_1e-4": {"mean": p1e4["mean"], "sd": p1e4["sd"]},
            "alpha_0.05": {"mean": _summ(power_frac.get(0.05, []))["mean"],
                           "sd": _summ(power_frac.get(0.05, []))["sd"]},
            "lambda_mean": lam_s["mean"],
            "definition": ("power = 因果 SNP 中达到阈值的比例（每次重复内平均）；"
                           "any_hit_rate = 至少检出 1 个因果 SNP 的重复比例"),
        }
        # 二元性状：方向一致性（效应量尺度不可比，改比符号）
        if model == "logistic" and sign_obs:
            effect_by_h2[f"{h2:g}"] = {
                "h2": h2,
                "mean_beta_ratio": None,
                "median_beta_ratio": None,
                "se_calibration": None,
                "sign_concordance": float(np.mean(sign_obs)),
                "n_obs": int(len(sign_obs)),
                "definition": ("logistic：log-OR 与 liability-scale β 尺度不同，只比较方向；"
                               "sign_concordance = 估计 log-OR 与真实效应符号一致的因果位点比例 "
                               "（理想≈1）"),
            }
        # 效应量偏差
        if beta_pairs:
            dfe = pd.DataFrame(beta_pairs)
            # 每个 SNP 跨重复的均值
            per_snp = dfe.groupby("snp").agg(beta_true=("beta_true", "mean"),
                                             beta_est=("beta_est", "mean"),
                                             se_est=("se_est", "mean"),
                                             n_reps=("beta_est", "size"),
                                             beta_sd=("beta_est", "std")).reset_index()
            per_snp = per_snp.loc[per_snp["beta_true"].abs() > 1e-8].copy()
            per_snp["ratio"] = per_snp["beta_est"] / per_snp["beta_true"]
            per_snp["h2"] = h2
            effect_rows.extend(per_snp.to_dict("records"))
            ratios = per_snp["ratio"].to_numpy(dtype=float)
            ratios = ratios[np.isfinite(ratios)]
            # 汇总比值用回归斜率 Σ(β_est·β_true)/Σ(β_true²)：
            # 符号感知、由大效应位点主导，避免因果效应正负相消时 Σest/Σtrue 不稳定
            b_est_v = per_snp["beta_est"].to_numpy(dtype=float)
            b_true_v = per_snp["beta_true"].to_numpy(dtype=float)
            denom = float(np.nansum(b_true_v ** 2))
            aggregate_ratio = (float(np.nansum(b_est_v * b_true_v)) / denom
                               if denom > 1e-12 else float("nan"))
            se_calib = float(np.nanmedian(
                per_snp["se_est"].to_numpy(dtype=float) /
                np.maximum(per_snp["beta_sd"].to_numpy(dtype=float), 1e-9)))
            effect_by_h2[f"{h2:g}"] = {
                "h2": h2,
                "mean_beta_ratio": aggregate_ratio,
                "median_beta_ratio": float(np.median(ratios)) if ratios.size else float("nan"),
                "sd_beta_ratio": float(np.std(ratios, ddof=1)) if ratios.size > 1 else 0.0,
                "se_calibration": se_calib,
                "n_causal_snps": int(len(per_snp)),
                "definition": ("mean_beta_ratio = 回归斜率 Σ(β_est·β_true)/Σ(β_true²)（符号感知、由主效 QTL 主导的汇总比值）；"
                               "median_beta_ratio = 逐位点 β_est/β_true 的中位数；"
                               "se_calibration = 各位点平均 SE / 其 β_est 跨重复 SD（理想≈1）"),
            }

    primary = f"{h2_primary:g}" if f"{h2_primary:g}" in power_by_h2 else list(power_by_h2)[-1]
    power_primary = dict(power_by_h2[primary])
    effect_primary = dict(effect_by_h2.get(primary) or {})
    conc = {
        "rank_auc": _summ([r["rank_auc"] for r in concordance_rows])["mean"],
        "top_causal_in_top100": _summ([r["top_causal_in_top100"]
                                       for r in concordance_rows])["mean"],
        "top1_is_causal_rate": _summ([r["top1_is_causal"]
                                      for r in concordance_rows])["mean"],
        "n_observations": len(concordance_rows),
        "definition": ("rank_auc = 因果 SNP 的 p 值小于随机 SNP 的概率（0.5=无信息）；"
                       "top_causal_in_top100 = 排名前 100 位点中因果 SNP 占全部因果 SNP 的比例"),
    }

    # --- 验收 gates ---------------------------------------------------------
    gate_rows: list[dict[str, Any]] = []

    def gate(name: str, value: Any, ok_range, detail: str) -> None:
        if value is None or not np.isfinite(value):
            status = "warn"
        else:
            lo, hi = ok_range
            status = "pass" if (lo is None or value >= lo) and (hi is None or value <= hi) \
                else "warn"
        gate_rows.append({"name": name, "status": status,
                          "value": value, "range": ok_range, "detail": detail})

    lam_mean = null["lambda_mean"]
    gate("null_lambda_mean", lam_mean, lambda_pass,
         f"零模型 λ 均值应在 {lambda_pass}")
    fpr05 = null["type1"].get("alpha_0.05")
    gate("null_fpr_alpha_0.05", fpr05, gates["fpr_alpha_0.05"],
         f"α=0.05 经验 FPR 应在 {gates['fpr_alpha_0.05']}")
    fpr5e8 = null["type1"].get("alpha_5e-8")
    gate("null_fpr_alpha_5e-08", fpr5e8, (0.0, 1e-3),
         "α=5e-8 经验 FPR 应接近 0")
    if np.isfinite(conc["rank_auc"]):
        gate("causal_rank_auc", conc["rank_auc"], (gates["rank_auc_min"], 1.0),
             f"因果 SNP 排序 AUC 应 > {gates['rank_auc_min']}")
    if np.isfinite(conc["top_causal_in_top100"]):
        gate("top_causal_in_top100", conc["top_causal_in_top100"],
             (gates["top100_min"], 1.0),
             f"top100 因果命中率应 > {gates['top100_min']}")
    any_hit = (power_primary.get("alpha_5e-8") or {}).get("any_hit_rate")
    if any_hit is not None:
        gate("qtl_any_hit_rate_5e-08", any_hit, (0.0, 1.0),
             "h²=主设定时至少检出 1 个因果 SNP 的重复比例（信息性，不作硬性门槛）")
    n_warn = sum(1 for g in gate_rows if g["status"] != "pass")
    null_lambda_status = "pass"
    if lam_mean is not None and np.isfinite(lam_mean):
        lo_w, hi_w = gates["lambda_warn"]
        null_lambda_status = ("pass" if lambda_pass[0] <= lam_mean <= lambda_pass[1]
                              else ("warn" if lo_w <= lam_mean <= hi_w else "fail"))

    # --- 表格与图 -----------------------------------------------------------
    tables: dict[str, str] = {}
    try:
        tables["calibration_null"] = _rel_to(t, t.tables_dir / "calibration_null.csv")
        utils.write_csv(pd.DataFrame(null_rows), t.tables_dir / "calibration_null.csv")
        tables["calibration_power"] = _rel_to(t, t.tables_dir / "calibration_power.csv")
        utils.write_csv(pd.DataFrame(power_rows), t.tables_dir / "calibration_power.csv")
        if effect_rows:
            tables["calibration_effect_bias"] = _rel_to(
                t, t.tables_dir / "calibration_effect_bias.csv")
            utils.write_csv(pd.DataFrame(effect_rows),
                            t.tables_dir / "calibration_effect_bias.csv")
        tables["calibration_concordance"] = _rel_to(
            t, t.tables_dir / "calibration_concordance.csv")
        utils.write_csv(pd.DataFrame(concordance_rows),
                        t.tables_dir / "calibration_concordance.csv")
        tables["calibration_gates"] = _rel_to(t, t.tables_dir / "calibration_gates.csv")
        utils.write_csv(pd.DataFrame(gate_rows), t.tables_dir / "calibration_gates.csv")
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"写校准表格失败: {type(exc).__name__}: {exc}")
    figures: list[str] = []
    try:
        p = _plot(t, "null", null_lambdas, t.cal_dir / "plots" / "calibration_null_lambda.png")
        if p:
            figures.append(_rel_to(t, p))
        alphas = sorted(ALPHAS)
        prim_rows = [r for r in power_rows if abs(r["h2"] - h2_primary) < 1e-9]
        pdp = pd.DataFrame({
            "alpha": alphas,
            "power": [float(np.nanmean([r.get(f"power_{a:g}", np.nan) for r in prim_rows]))
                      for a in alphas]})
        p = _plot(t, "power", pdp, t.cal_dir / "plots" / "calibration_power.png")
        if p:
            figures.append(_rel_to(t, p))
        if effect_rows:
            dfe = pd.DataFrame(effect_rows)
            p = _plot(t, "effect", dfe, t.cal_dir / "plots" / "calibration_effect_bias.png")
            if p:
                figures.append(_rel_to(t, p))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"绘图失败: {type(exc).__name__}: {exc}")

    # --- 汇总状态 -----------------------------------------------------------
    if null_lambda_status == "fail":
        status = "fail"
        errors.append(f"零模型 λ 均值 {lam_mean:.3f} 超出可接受范围 {gates['lambda_warn']}")
    elif null_lambda_status == "warn" or n_warn > 0:
        status = "warn"
    else:
        status = "pass"
    if fpr05 is not None and np.isfinite(fpr05) and not (
            gates["fpr_warn"][0] <= fpr05 <= gates["fpr_warn"][1]):
        errors.append(f"α=0.05 的 FPR={fpr05:.4f} 明显偏离 0.05"

                      f"（可接受 {gates['fpr_warn']}）")
        status = "fail"

    metrics: dict[str, Any] = {
        "stage": "calibration",
        "status": status,
        "model": model,
        "phenotype_type": pheno_type,
        "n_samples": int(n_samples),
        "n_snps_used": int(len(raw.snps)),
        "snp_source": snp_source,
        "n_replicates_requested": n_reps,
        "n_replicates_done": reps_eff,
        "h2_scenarios": h2_list,
        "seed": seed,
        "truth_source": eff.get("truth_path"),
        "effect_source": eff.get("source"),
        "n_causal_used": int(len(eff["causal"])),
        "n_background_snps": int(eff.get("n_bg") or 0),
        "null": null,
        "power": power_primary,
        "power_by_h2": power_by_h2,
        "effect_bias": effect_primary or None,
        "effect_bias_by_h2": effect_by_h2,
        "concordance": conc,
        "lambda_qtl": lambda_qtl,
        "gates": gate_rows,
        "thresholds": {"alpha_gw": ALPHA_GW, "alphas": list(ALPHAS),
                       "lambda_pass": list(lambda_pass),
                       "fpr_alpha_0.05": list(gates["fpr_alpha_0.05"])},
        "tables": tables,
        "figures": figures,
        "warnings": warnings,
        "errors": errors,
    }
    metrics["summary"] = _summary(metrics)
    return _finish(t, started, metrics)


def _rel_to(t, path: Path) -> str:
    try:
        return str(Path(path).relative_to(t.run_dir)).replace("\\", "/")
    except ValueError:
        return str(path)


def _summary(m: dict[str, Any]) -> str:
    lines: list[str] = []
    null = m.get("null") or {}
    if null:
        lines.append(f"零模型校准（{null.get('n_replicates')} 次重复）：λGC 均值 "
                     f"{null.get('lambda_mean'):.4f}（SD {null.get('lambda_sd'):.4f}）；"
                     f"经验假阳性率 α=0.05 为 {null.get('type1', {}).get('alpha_0.05', float('nan')):.4f}、"
                     f"α=1e-3 为 {null.get('type1', {}).get('alpha_1e-3', float('nan')):.5f}、"
                     f"α=5e-8 为 {null.get('type1', {}).get('alpha_5e-8', float('nan')):.2g}。")
    power = m.get("power") or {}
    if power:
        p5 = power.get("alpha_5e-8") or {}
        lines.append(f"已知 QTL（h²={power.get('h2')}，{power.get('n_causal')} 个因果位点）："
                     f"5e-8 阈值下平均检出 {p5.get('mean', float('nan')):.3f} 的因果位点，"
                     f"{p5.get('any_hit_rate', float('nan')):.2f} 的重复至少检出 1 个；"
                     f"1e-4 阈值平均检出 "
                     f"{(power.get('alpha_1e-4') or {}).get('mean', float('nan')):.3f}。")
    eb = m.get("effect_bias") or {}
    if eb:
        if eb.get("mean_beta_ratio") is not None:
            lines.append(f"效应量标定：β_est/β_true 回归斜率 {eb.get('mean_beta_ratio'):.3f}，"
                         f"SE 标定比 {eb.get('se_calibration'):.2f}（理想≈1）。")
        elif eb.get("sign_concordance") is not None:
            lines.append(f"效应方向标定：logistic 下估计 log-OR 与真实效应符号一致率 "
                         f"{eb.get('sign_concordance'):.3f}（理想≈1；log-OR 与 liability β "
                         "尺度不同，只比方向）。")
    conc = m.get("concordance") or {}
    if conc:
        lines.append(f"排序一致性：因果 SNP rank AUC={conc.get('rank_auc'):.3f}，"
                     f"top100 命中率={conc.get('top_causal_in_top100'):.3f}。")
    lines.append("结论：本数据集的 p 值校准与检出能力如上；λ 明显偏离 1 或 FPR 偏离 α 时，"
                 "GWAS 结果需谨慎解释（可结合更严格质控/协变量再分析）。")
    return "\n".join(lines)


def _finish(t, started: float, metrics: dict[str, Any]) -> dict[str, Any]:
    metrics["elapsed_s"] = round(time.time() - started, 2)
    try:
        t.cal_dir.mkdir(parents=True, exist_ok=True)
        utils.write_json(t.cal_dir / "metrics.json", metrics)
    except Exception as exc:  # noqa: BLE001
        t.logger.warning("写 calibration metrics.json 失败: %s", exc)
    for msg in metrics.get("warnings") or []:
        t.logger.warning("[calibration] %s", msg)
    for msg in metrics.get("errors") or []:
        t.logger.error("[calibration] %s", msg)
    t.logger.info("仿真校准完成：status=%s，耗时 %.1fs",
                  metrics.get("status"), metrics.get("elapsed_s") or 0.0)
    return metrics


if __name__ == "__main__":  # pragma: no cover - 手工对 run 目录执行
    import argparse
    import json

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description="对已有 run 目录执行仿真校准")
    ap.add_argument("run_dir")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    res = run_calibration(a.run_dir)
    print(f"status={res.get('status')}  elapsed={res.get('elapsed_s')}s")
    print(res.get("summary", ""))
    if res.get("null"):
        print("null:", json.dumps({k: res["null"][k] for k in
                                   ("n_replicates", "lambda_mean", "lambda_sd", "type1")},
                                  ensure_ascii=False))
    if res.get("power"):
        print("power:", json.dumps(res["power"], ensure_ascii=False))
    if a.json:
        print(json.dumps(res, ensure_ascii=False, indent=2, default=str))
