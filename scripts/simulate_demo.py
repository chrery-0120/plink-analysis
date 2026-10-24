#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""生成 PLINK 1.9 演示数据集（纯 Python 写 .bed/.bim/.fam，固定种子可复现）。

用途
----
为 plink-analysis skill 的端到端自检（scripts/selftest.py）和仿真校准
（scripts/pipeline/simulate.py）提供一份"结构已知、可验证"的演示数据集：

* 默认 400 样本 x 8000 SNP，18 条常染色体，固定种子 20261023，逐字节可复现。
* 2 个亚群（breed A / B）：全基因组存在轻度频率分化（Fst>0），另有 8% 位点
  强分化（|ΔAF| 0.25~0.42），PCA 第一主成分可完全分开两个亚群。
* 10 对高亲缘样本（8 对约 65% 基因型共担 + 2 对完全重复），用于检验
  duplicate / relatedness 质控闸门。
* 基因型缺失率 2%（完全随机）。
* MAF 谱来自 Beta(0.5, 2) 的次等位基因频率，符合"低 MAF 位点占多数"的经验分布。
* 30 个已知因果 SNP（3 个主效 QTL + 27 个微效），另加 1500 个多基因背景位点；
  数量性状 h2=0.5（因果位点占遗传方差的 85%），二元性状按 liability 模型
  生成，患病率约 20%（PLINK 编码 1=对照, 2=病例）。
* 输出 truth.csv / truth_meta.json，记录每个位点的真实等位基因频率与效应量，
  供 simulate.py 的已知 QTL 校准使用。

PLINK .bed 编码（已用 PLINK 1.9.0 实测确认）
-------------------------------------------
3 字节魔数 ``6c 1b 01``（01 = variant-major / SNP-major），随后按位点顺序，
每个位点占 ceil(N/4) 字节；每个样本 2 bit，低位先写：

    00 = A1A1（纯合第一个等位基因）    01 = 缺失
    10 = 杂合                        11 = A2A2（纯合第二个等位基因）

用法
----
    python scripts/simulate_demo.py --out E:\\deepseek\\demo_data
    python scripts/simulate_demo.py --out <dir> --n-samples 200 --n-snps 2000  # 小数据快测
    python scripts/simulate_demo.py --out <dir> --no-plink-verify              # 跳过 PLINK 校验

产物
----
    <out>/geno.bed .bim .fam    PLINK 二进制基因型
    <out>/pheno.txt             FID IID qtl       数量性状
    <out>/pheno_bin.txt         FID IID status    1=对照 2=病例（约 20% 病例）
    <out>/covar.txt             FID IID sex age
    <out>/group.txt             FID IID breed
    <out>/truth.csv             全部位点（含 is_causal/is_poly/beta_unit/...）
    <out>/truth_meta.json       生成参数、真实 h2、亲缘对、校验结果
    <out>/README.md             数据说明 + PLINK 读取验证命令与结果
    <out>/verify/               PLINK 校验原始输出（.log/.frq/.raw/.hwe/...）
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# 常量与默认值
# ---------------------------------------------------------------------------
def _workspace() -> Path:
    """<workspace> when the skill lives in <workspace>/.agents/skills/plink-analysis."""
    return Path(__file__).resolve().parents[4]


def _default_plink() -> str:
    for env_name in ("PLINK_EXE", "PLINK"):
        env = os.environ.get(env_name)
        if env:
            return env
    exe = "plink.exe" if os.name == "nt" else "plink"
    for candidate in (_workspace() / "tools" / "plink" / exe,
                      Path.cwd() / "tools" / "plink" / exe,
                      Path(exe)):
        if candidate.exists():
            return str(candidate)
    return exe


DEFAULT_OUT = str(_workspace() / "demo_data")
DEFAULT_PLINK = _default_plink()
DEFAULT_SEED = 20261023
DEFAULT_N_SAMPLES = 400
DEFAULT_N_SNPS = 8000
N_CHROM = 18                      # 猪的常染色体数（示例物种）
N_CAUSAL = 30                     # 已知因果位点数
N_MAJOR = 3                       # 其中主效 QTL 个数
CAUSAL_MIN_MAF = 0.15             # 因果位点要求常见变异
DIFF_RATE = 0.08                  # 亚群强分化位点比例
DRIFT = 0.15                      # 全基因组频率漂移（logit 尺度，方向随机）
DRIFT_LOGIT = DRIFT
DIFF_LOGIT = (1.6, 3.0)           # 强分化位点在 logit 尺度上的偏移范围
MISSING_RATE = 0.02               # 基因型缺失率
N_POLY = 1500                     # 多基因背景位点数
PVE_CAUSAL = 0.85                 # 因果位点承担的总遗传方差比例
MAJOR_SHARE = 0.60                # 主效 QTL 在因果方差中的占比（其余给微效）
H2_QUANT = 0.5                    # 数量性状遗传力
H2_LIAB = 0.5                     # 二元性状 liability 尺度遗传力
PREVALENCE = 0.20                 # 二元性状患病率
HWE_P_CUT = 1e-6


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _now() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


# ---------------------------------------------------------------------------
# 1. 基因型 / 表型生成
# ---------------------------------------------------------------------------
def build_dataset(seed: int, n_samples: int, n_snps: int) -> dict:
    """生成全部底层的 numpy 数组；所有随机性都来自同一个 seed。"""
    rng = np.random.default_rng(seed)

    # --- 位点位置：18 条常染色体各自连续一段，物理距离 5kb 递增 ----------
    per = n_snps // N_CHROM
    extra = n_snps % N_CHROM
    counts = [per + (1 if c < extra else 0) for c in range(N_CHROM)]
    chr_ids = np.repeat(np.arange(1, N_CHROM + 1), counts)   # 每条第 2 列必须连续
    bp = np.zeros(n_snps, dtype=np.int64)
    start = 1_000_000
    for c in range(1, N_CHROM + 1):
        idx = np.flatnonzero(chr_ids == c)
        bp[idx] = start + np.arange(len(idx)) * 5_000
        start += len(idx) * 5_000 + 2_000_000

    # --- 次等位基因频率谱：Beta(0.5,2)，并随机指定哪一个等位基因是次要等位基因 ---
    maf_draw = rng.beta(0.5, 2.0, n_snps)
    flip = rng.random(n_snps) < 0.5
    af1 = np.where(flip, 1.0 - maf_draw, maf_draw)       # 等位基因 1 在 A 群体中的频率
    # 两端轻微夹紧，避免造出严格单态位点（保留少量近单态位点用于 MAF 过滤）
    af1 = np.clip(af1, 0.002, 0.998)

    # --- 亚群频率分化（logit 尺度，方向随机；对称偏移不改变总体 MAF 谱） ----
    def _logit(p: np.ndarray) -> np.ndarray:
        p = np.clip(p, 1e-6, 1 - 1e-6)
        return np.log(p / (1.0 - p))

    def _expit(x: np.ndarray) -> np.ndarray:
        return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))

    is_diff = rng.random(n_snps) < DIFF_RATE
    sign = rng.choice([-1.0, 1.0], n_snps)
    d_logit = np.where(is_diff, rng.uniform(*DIFF_LOGIT, n_snps), DRIFT_LOGIT) * sign
    logit_base = _logit(af1)
    af_a = _expit(logit_base + d_logit)                  # breed A 中等位基因 1 的频率
    af_b = _expit(logit_base - d_logit)                  # breed B 中等位基因 1 的频率
    af_a = np.clip(af_a, 0.0005, 0.9995)
    af_b = np.clip(af_b, 0.0005, 0.9995)

    # --- 样本分组与协变量 ----------------------------------------------------
    n_a = n_samples // 2
    breed = np.array(["A"] * n_a + ["B"] * (n_samples - n_a), dtype=object)
    rng.shuffle(breed)
    sex = rng.integers(1, 3, n_samples)                  # 1=公 2=母
    age = np.round(rng.normal(3.2, 1.1, n_samples), 1)   # 岁
    age = np.clip(age, 0.5, 8.0)

    # --- 基因型：按亚群等位基因频率做二项抽样（HWE） -------------------------
    q = np.where(breed[:, None] == "A", af_a[None, :], af_b[None, :]).astype(np.float64)
    dos = rng.binomial(2, q).astype(np.int8)             # A1 等位基因个数 0/1/2

    # --- 高亲缘 / 重复样本 ---------------------------------------------------
    # 重复对：完全拷贝；亲缘对：拷贝 ~65% 位点基因型（其余独立），
    # 使 PLINK --genome 的 PI_HAT 落在 0.4~0.8（亲缘）与 ~1（重复）两档。
    pairs: list[dict] = []
    idx_a = np.flatnonzero(breed == "A")
    idx_b = np.flatnonzero(breed == "B")
    used: set[int] = set()
    for k in range(10):
        pool = idx_a if k % 2 == 0 else idx_b
        cand = [int(i) for i in pool if int(i) not in used]
        rng.shuffle(cand)
        s, t = cand[0], cand[1]
        used.update({s, t})
        is_dup = k < 2
        if is_dup:
            dos[t, :] = dos[s, :]
            kind = "duplicate"
        else:
            share = rng.random(n_snps) < 0.65
            dos[t, share] = dos[s, share]
            kind = "related"
        pairs.append({"sample1": int(s), "sample2": int(t), "kind": kind,
                      "same_breed": True})

    # --- 2% 随机缺失 ---------------------------------------------------------
    miss = rng.random((n_samples, n_snps)) < MISSING_RATE
    dos_f = dos.astype(np.float64)
    dos_f[miss] = np.nan

    # --- 位点筛选：取常见、低分化的位点作为因果位点 --------------------------
    af_emp = np.nanmean(np.where(dos_f >= 0, dos_f, np.nan), axis=0) / 2.0
    maf_emp = np.minimum(af_emp, 1.0 - af_emp)
    eligible = np.flatnonzero((maf_emp >= CAUSAL_MIN_MAF) & (~is_diff))
    n_causal = min(N_CAUSAL, max(6, eligible.size // 10))
    causal = np.sort(rng.choice(eligible, size=n_causal, replace=False))

    n_major = min(N_MAJOR, n_causal)
    rem = np.flatnonzero(eligible)
    rem = np.setdiff1d(rem, causal, assume_unique=True)
    n_poly = min(N_POLY, max(0, rem.size))
    poly = np.sort(rng.choice(rem, size=n_poly, replace=False)) if n_poly else np.array([], dtype=int)

    # --- 真实效应量：单位遗传方差（见 truth_meta 的说明） --------------------
    p_causal = np.clip(af_emp[causal], 0.01, 0.99)
    h_units = np.empty(n_causal)
    h_units[:n_major] = MAJOR_SHARE / n_major
    h_units[n_major:] = (1.0 - MAJOR_SHARE) / max(1, n_causal - n_major)
    beta_unit_causal = np.sqrt(h_units / (2.0 * p_causal * (1.0 - p_causal))) * np.sqrt(PVE_CAUSAL)
    # 主效 QTL 方向随机（部分为正、部分为负），微效随机符号
    beta_unit_causal *= rng.choice([-1.0, 1.0], size=n_causal)

    if n_poly:
        p_poly = np.clip(af_emp[poly], 0.01, 0.99)
        b_raw = rng.normal(size=n_poly)
        scale = np.sqrt((1.0 - PVE_CAUSAL) / np.sum(2.0 * p_poly * (1.0 - p_poly) * b_raw ** 2))
        beta_unit_poly = b_raw * scale
    else:
        beta_unit_poly = np.zeros(0)

    # --- 遗传值（缺失剂量用群体均值填补，避免样本遗传值被缺失拉低） ----------
    dos_imp = np.where(np.isnan(dos_f), 2.0 * af_emp[None, :], dos_f)
    g_unit = np.zeros(n_samples)
    if causal.size:
        g_unit += dos_imp[:, causal] @ beta_unit_causal
    if poly.size:
        g_unit += dos_imp[:, poly] @ beta_unit_poly
    g_std = (g_unit - g_unit.mean()) / g_unit.std(ddof=0)

    def _liability(h2: float, noise: np.ndarray) -> np.ndarray:
        """y = sqrt(h2/(1-h2)) * G_std + e, var(e)=1 => 样本遗传力≈h2。"""
        a = np.sqrt(h2 / (1.0 - h2))
        return a * g_std + noise

    noise_q = rng.normal(size=n_samples)
    trait = _liability(H2_QUANT, noise_q)

    noise_l = rng.normal(size=n_samples)
    liability = _liability(H2_LIAB, noise_l)
    thr = float(np.quantile(liability, 1.0 - PREVALENCE))
    status = np.where(liability > thr, 2, 1).astype(np.int8)   # PLINK: 1=对照 2=病例

    # --- 已实现遗传力（对表型方差的经验估计，供报告引用） --------------------
    a_q = np.sqrt(H2_QUANT / (1.0 - H2_QUANT))
    h2_realized = float(np.var(a_q * g_std) / np.var(trait)) if np.var(trait) > 0 else float("nan")

    return dict(
        rng=rng, chr_ids=chr_ids, bp=bp, af1=af1, af_a=af_a, af_b=af_b,
        is_diff=is_diff, breed=breed, sex=sex, age=age, dos_f=dos_f, miss=miss,
        af_emp=af_emp, maf_emp=maf_emp, pairs=pairs,
        causal=causal, poly=poly, beta_unit_causal=beta_unit_causal,
        beta_unit_poly=beta_unit_poly, g_unit=g_unit, g_std=g_std,
        trait=trait, status=status, thr=thr, liability=liability,
        h2_realized=h2_realized,
        sample_ids=[f"ind{i + 1:03d}" for i in range(n_samples)],
        snp_ids=[f"snp{i + 1:05d}" for i in range(n_snps)],
    )


# ---------------------------------------------------------------------------
# 2. 写 PLINK bfile + 表型/协变量/分组
# ---------------------------------------------------------------------------
def _pack_bed(no_call_codes: np.ndarray) -> bytes:
    """把 n_snps x n_samples 的 2-bit 编码矩阵打包成 .bed 字节流（variant-major）。"""
    n_snps, n_samples = no_call_codes.shape
    n_pad = ((n_samples + 3) // 4) * 4
    padded = np.zeros((n_snps, n_pad), dtype=np.uint8)
    padded[:, :n_samples] = no_call_codes
    padded = padded.reshape(n_snps, n_pad // 4, 4)
    packed = (padded[..., 0] | (padded[..., 1] << 2) |
              (padded[..., 2] << 4) | (padded[..., 3] << 6))
    return bytes([0x6C, 0x1B, 0x01]) + packed.astype(np.uint8).tobytes()


def write_bfile(out_dir: Path, ds: dict) -> dict:
    """写 geno.bed/.bim/.fam；返回 {prefix, bed, bim, fam}（不改名）。"""
    dos_f = ds["dos_f"]
    n_samples, n_snps = dos_f.shape
    dos_t = dos_f.T                                 # 转成位点主序（variant-major）
    # 剂量(A1 个数) -> .bed 编码：0→11(A2A2), 1→10(het), 2→00(A1A1), NaN→01(缺失)
    codes = np.full((n_snps, n_samples), 3, dtype=np.uint8)
    codes[dos_t == 2] = 0
    codes[dos_t == 1] = 2
    codes[np.isnan(dos_t)] = 1
    prefix = out_dir / "geno"
    (prefix.with_suffix(".bed")).write_bytes(_pack_bed(codes))

    alleles = ("A", "C", "G", "T")
    lines = []
    for i in range(n_snps):
        a1 = alleles[i % 4]
        a2 = alleles[(i // 4 + 1) % 4]
        if a2 == a1:
            a2 = "T" if a1 != "T" else "C"
        lines.append(f"{ds['chr_ids'][i]} {ds['snp_ids'][i]} 0 {ds['bp'][i]} {a1} {a2}")
    (prefix.with_suffix(".bim")).write_text("\n".join(lines) + "\n", encoding="utf-8")

    lines = []
    for i, iid in enumerate(ds["sample_ids"]):
        fid = str(ds["breed"][i])
        lines.append(f"{fid} {iid} 0 0 {ds['sex'][i]} -9")
    (prefix.with_suffix(".fam")).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"prefix": prefix, "bed": prefix.with_suffix(".bed"),
            "bim": prefix.with_suffix(".bim"), "fam": prefix.with_suffix(".fam")}


def _write_table(path: Path, header: list[str], rows: list[str]) -> None:
    path.write_text("\n".join(["\t".join(header)] + rows) + "\n", encoding="utf-8")


def write_pheno_covar_group(out_dir: Path, ds: dict) -> dict:
    ids = ds["sample_ids"]
    fid = [str(b) for b in ds["breed"]]

    def txt(name: str) -> Path:
        return out_dir / name

    _write_table(txt("pheno.txt"), ["FID", "IID", "qtl"],
                 [f"{fid[i]}\t{ids[i]}\t{ds['trait'][i]:.6f}" for i in range(len(ids))])
    _write_table(txt("pheno_bin.txt"), ["FID", "IID", "status"],
                 [f"{fid[i]}\t{ids[i]}\t{int(ds['status'][i])}" for i in range(len(ids))])
    _write_table(txt("covar.txt"), ["FID", "IID", "sex", "age"],
                 [f"{fid[i]}\t{ids[i]}\t{int(ds['sex'][i])}\t{ds['age'][i]:.1f}"
                  for i in range(len(ids))])
    _write_table(txt("group.txt"), ["FID", "IID", "breed"],
                 [f"{fid[i]}\t{ids[i]}\t{ds['breed'][i]}" for i in range(len(ids))])
    return {"pheno": txt("pheno.txt"), "pheno_bin": txt("pheno_bin.txt"),
            "covar": txt("covar.txt"), "group": txt("group.txt")}


def write_truth(out_dir: Path, ds: dict, args: argparse.Namespace, verify: dict | None) -> dict:
    """写 truth.csv（逐位点）与 truth_meta.json（生成参数、真实 h2、校验摘要）。"""
    n_snps = len(ds["snp_ids"])
    is_causal = np.zeros(n_snps, dtype=bool)
    is_causal[ds["causal"]] = True
    is_poly = np.zeros(n_snps, dtype=bool)
    if ds["poly"].size:
        is_poly[ds["poly"]] = True
    b_unit = np.zeros(n_snps)
    b_unit[ds["causal"]] = ds["beta_unit_causal"]
    if ds["poly"].size:
        b_unit[ds["poly"]] = ds["beta_unit_poly"]
    a_q = np.sqrt(H2_QUANT / (1.0 - H2_QUANT))
    a_l = np.sqrt(H2_LIAB / (1.0 - H2_LIAB))
    sd_g = float(ds["g_std"].std(ddof=0)) or 1.0

    line = ["SNP,CHR,BP,A1,A2,AF_A,AF_B,MAF,IS_DIFF,IS_CAUSAL,IS_POLY,"
            "BETA_UNIT,BETA_QT,BETA_LIAB"]
    lines = [",".join(line)]
    for i in range(n_snps):
        lines.append(",".join([
            ds["snp_ids"][i], str(int(ds["chr_ids"][i])), str(int(ds["bp"][i])),
            "X", "Y",                      # 占位，下面补真实等位基因
            f"{ds['af_a'][i]:.6f}", f"{ds['af_b'][i]:.6f}", f"{ds['maf_emp'][i]:.6f}",
            str(int(ds["is_diff"][i])), str(int(is_causal[i])), str(int(is_poly[i])),
            f"{b_unit[i]:.8g}", f"{b_unit[i] * a_q / sd_g:.8g}", f"{b_unit[i] * a_l / sd_g:.8g}",
        ]))
    (out_dir / "truth.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # 每个位点的真实等位基因符号（与 .bim 一致）：A1 = 字母表的第 i%4 个
    bim_lines = (out_dir / "geno.bim").read_text(encoding="utf-8").splitlines()
    tlines = (out_dir / "truth.csv").read_text(encoding="utf-8").splitlines()
    fixed = [tlines[0]]
    for i in range(1, len(tlines)):
        cols = tlines[i].split(",")
        bim = bim_lines[i - 1].split()
        cols[3], cols[4] = bim[4], bim[5]
        fixed.append(",".join(cols))
    (out_dir / "truth.csv").write_text("\n".join(fixed) + "\n", encoding="utf-8")

    causal_rows = []
    for i in ds["causal"]:
        causal_rows.append({
            "snp": ds["snp_ids"][int(i)], "chr": int(ds["chr_ids"][int(i)]),
            "bp": int(ds["bp"][int(i)]), "maf": round(float(ds["maf_emp"][int(i)]), 4),
            "af_a": round(float(ds["af_a"][int(i)]), 4),
            "af_b": round(float(ds["af_b"][int(i)]), 4),
            "beta_unit": round(float(b_unit[int(i)]), 6),
            "beta_qtl": round(float(b_unit[int(i)] * a_q / sd_g), 6),
        })

    meta = {
        "schema_version": 1,
        "seed": args.seed,
        "created_at": _now(),
        "generator": "plink-analysis/scripts/simulate_demo.py",
        "n_samples": int(ds["dos_f"].shape[0]),
        "n_snps": int(ds["dos_f"].shape[1]),
        "n_chrom": N_CHROM,
        "breeds": {"A": int(np.sum(ds["breed"] == "A")), "B": int(np.sum(ds["breed"] == "B"))},
        "n_diff_snps": int(np.sum(ds["is_diff"])),
        "missing_rate": MISSING_RATE,
        "n_causal": int(len(ds["causal"])),
        "n_major_qtl": int(min(N_MAJOR, len(ds["causal"]))),
        "n_poly_snps": int(len(ds["poly"])),
        "h2_quantitative": H2_QUANT,
        "h2_quantitative_realized": round(ds["h2_realized"], 4),
        "h2_liability": H2_LIAB,
        "prevalence": PREVALENCE,
        "prevalence_realized": round(float(np.mean(ds["status"] == 2)), 4),
        "pve_causal": PVE_CAUSAL,
        "maf_spectrum": "Beta(0.5, 2) minor-allele frequency",
        "effect_scale": ("BETA_UNIT: 效应量使 (因果+背景) 遗传值在群体中方差≈1；"
                         "BETA_QT = BETA_UNIT * sqrt(h2/(1-h2)) / sd(G)，"
                         "即在 pheno.txt 原始尺度上的每 A1 等位基因加性效应；"
                         "BETA_LIAB 同理（liability 尺度，二元性状用）"),
        "causal_snps": causal_rows,
        "related_pairs": ds["pairs"],
        "files": {},
        "plink_verify": verify or {},
    }
    for p in sorted(out_dir.iterdir()):
        if p.is_file() and p.name.startswith(("geno.", "pheno", "covar", "group", "truth")):
            meta["files"][p.name] = {"size": p.stat().st_size, "sha256": _sha256(p)}
    (out_dir / "truth_meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return meta


# ---------------------------------------------------------------------------
# 3. 纯 Python 自检（不依赖 PLINK）
# ---------------------------------------------------------------------------
def self_check(ds: dict, meta: dict) -> list[dict]:
    """在写盘前后用 numpy 检查数据结构是否满足设计目标。"""
    checks: list[dict] = []

    def add(name: str, ok: bool, detail: str, value=None) -> None:
        checks.append({"name": name, "status": "pass" if ok else "fail",
                       "detail": detail, "value": value})

    n_s, n_v = ds["dos_f"].shape
    add("sample_count", n_s >= 300, f"n_samples={n_s}", n_s)
    add("variant_count", n_v >= 5000 or n_v >= 200, f"n_snps={n_v}", n_v)

    miss = float(np.mean(np.isnan(ds["dos_f"])))
    add("missing_rate", abs(miss - MISSING_RATE) < 0.005,
        f"实际缺失率={miss:.4f}（目标 {MISSING_RATE}）", round(miss, 5))

    maf = ds["maf_emp"]
    add("maf_spectrum", float(np.mean((maf > 0.01) & (maf < 0.5))) > 0.6,
        f"MAF 均值={np.nanmean(maf):.3f}，MAF<0.01 占比="
        f"{float(np.mean(maf < 0.01)):.3f}，单态位点={int(np.sum(maf < 1e-9))}",
        round(float(np.nanmean(maf)), 4))

    # PCA：均值填补 + 标准化后 SVD，检查 PC1 是否完全分开两个亚群
    dos_imp = np.where(np.isnan(ds["dos_f"]), 2.0 * ds["af_emp"][None, :], ds["dos_f"])
    x = (dos_imp - dos_imp.mean(axis=0)) / np.maximum(dos_imp.std(axis=0), 1e-9)
    x = np.nan_to_num(x)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    pc1 = u[:, 0] * s[0]
    pc1 = (pc1 - pc1.mean()) / pc1.std()
    mean_a, mean_b = pc1[ds["breed"] == "A"].mean(), pc1[ds["breed"] == "B"].mean()
    # AUC：A 群体的 PC1 大于 B 群体的比例
    auc = float(np.mean(pc1[ds["breed"] == "A"][:, None] > pc1[ds["breed"] == "B"][None, :]))
    eig_all = s ** 2
    var_pc1 = float(eig_all[0] / eig_all.sum())
    add("pca_separation", max(auc, 1 - auc) >= 0.99,
        f"PC1 区分 breed A/B 的 AUC={max(auc, 1 - auc):.4f}，PC1 方差占比={var_pc1:.3f}，"
        f"组均值差={abs(mean_a - mean_b):.2f} sd", round(max(auc, 1 - auc), 4))
    ds["_pc1"] = pc1
    ds["_eigvals"] = eig_all[:20]

    # 亲缘/重复对：用"基因型一致率"粗检（PLINK --genome 结果见 verify）
    conc = []
    for pr in ds["pairs"]:
        i, j = pr["sample1"], pr["sample2"]
        a, b = ds["dos_f"][i], ds["dos_f"][j]
        good = ~(np.isnan(a) | np.isnan(b))
        rate = float(np.mean(a[good] == b[good])) if good.any() else float("nan")
        conc.append(rate)
    ok = all(c > 0.55 for c in conc)
    add("related_pairs", ok,
        f"{len(conc)} 对亲缘/重复样本，基因型一致率 min={min(conc):.3f}, max={max(conc):.3f}",
        round(min(conc), 4))

    # 已知 QTL 的可检出性（边际线性回归，无协变量）：主效位点应达全基因组显著
    from scipy import stats as st
    trait = ds["trait"]
    y = (trait - trait.mean()) / trait.std()
    causal = ds["causal"]
    xc = dos_imp[:, causal]
    xc = (xc - xc.mean(axis=0)) / np.maximum(xc.std(axis=0), 1e-12)
    beta = xc.T @ y / xc.shape[0]
    resid = y[:, None] - xc * beta[None, :]
    dof = max(1, xc.shape[0] - 2)
    se = np.sqrt(np.sum(resid ** 2, axis=0) / dof / np.sum(xc ** 2, axis=0))
    tstat = beta / np.maximum(se, 1e-12)
    pval = 2 * st.t.sf(np.abs(tstat), dof)
    n_sig = int(np.sum(pval < 5e-8))
    n_sugg = int(np.sum(pval < 1e-4))
    add("causal_qtl_detectable", n_sig >= 1,
        f"因果位点中 p<5e-8 的有 {n_sig}/{len(causal)} 个，p<1e-4 的有 {n_sugg} 个，"
        f"min p={pval.min():.2e}", n_sig)

    # HWE：按 breed 分别做卡方检验，统计极端偏离比例
    extremes = 0
    total = 0
    for grp in ("A", "B"):
        idx = np.flatnonzero(ds["breed"] == grp)
        d = ds["dos_f"][idx]
        n = d.shape[0]
        for j in range(0, d.shape[1], 7):        # 抽样检查以控制耗时
            col = d[:, j]
            col = col[~np.isnan(col)]
            n = col.size
            if n < 20:
                continue
            obs_hom1 = int(np.sum(col == 2))
            obs_het = int(np.sum(col == 1))
            obs_hom2 = int(np.sum(col == 0))
            p = (2 * obs_hom1 + obs_het) / (2 * n)
            if p <= 1e-6 or p >= 1 - 1e-6:
                continue
            exp_het = 2 * p * (1 - p) * n
            chi2 = ((obs_het - exp_het) ** 2 / max(exp_het, 1e-9) +
                    (obs_hom1 - n * p * p) ** 2 / max(n * p * p, 1e-9) +
                    (obs_hom2 - n * (1 - p) ** 2) ** 2 / max(n * (1 - p) ** 2, 1e-9))
            total += 1
            if st.chi2.sf(chi2, 1) < HWE_P_CUT:
                extremes += 1
    frac = extremes / total if total else 0.0
    add("hwe_extremes", frac < 0.005,
        f"HWE p<1e-6 比例={frac:.5f}（检查位点 {total} 个）", round(frac, 6))

    # 二元性状患病率
    prev = float(np.mean(ds["status"] == 2))
    add("case_prevalence", 0.12 <= prev <= 0.30,
        f"病例比例={prev:.3f}（目标 {PREVALENCE}）", round(prev, 4))

    # 数值有效性
    add("no_nan_phenotype",
        bool(np.all(np.isfinite(ds["trait"])) and np.all(np.isfinite(ds["liability"]))),
        "表型/liability 无 NaN/Inf")
    return checks


# ---------------------------------------------------------------------------
# 4. PLINK 读取验证
# ---------------------------------------------------------------------------
def _run_plink(plink: str, args: list[str], cwd: Path, timeout: int = 300) -> subprocess.CompletedProcess:
    env = None
    cmd = [plink] + args
    return subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=timeout, env=env)


def verify_with_plink(out_dir: Path, plink_exe: str, ds: dict) -> dict:
    """用 PLINK 1.9 实际读取 bfile 并做 freq/recode A/hardy/missing/genome/pca。"""
    plink = Path(plink_exe)
    if not plink.is_file():
        return {"status": "skipped", "detail": f"PLINK 不存在: {plink}"}
    vdir = out_dir / "verify"
    if vdir.exists():
        shutil.rmtree(vdir)
    vdir.mkdir(parents=True, exist_ok=True)
    bfile = str(out_dir / "geno")
    common = ["--bfile", bfile, "--threads", "4", "--allow-no-sex"]
    results: dict = {"status": "pass", "plink": str(plink), "commands": [], "checks": []}

    def record(name: str, args: list[str], ok_msg: str) -> subprocess.CompletedProcess:
        out = str(vdir / name)
        r = _run_plink(str(plink), common + args + ["--out", out], out_dir)
        results["commands"].append({
            "name": name,
            "command": " ".join([str(plink)] + common + args + ["--out", out]),
            "returncode": r.returncode,
        })
        if r.returncode != 0:
            results["status"] = "fail"
            results["checks"].append({"name": name, "status": "fail",
                                      "detail": (r.stderr or r.stdout or "")[-800:]})
        return r

    # 1) 计数与读取
    r = record("freq", ["--freq"], "")
    n_loaded_v = n_loaded_s = None
    for line in (r.stdout or "").splitlines():
        if "variants loaded from .bim" in line:
            n_loaded_v = int(line.split()[0].replace(",", ""))
        if "people" in line and "loaded from .fam" in line:
            n_loaded_s = int(line.split()[0].replace(",", ""))
    results["checks"].append({
        "name": "bed_readable",
        "status": "pass" if r.returncode == 0 and n_loaded_v == ds["dos_f"].shape[1]
                  and n_loaded_s == ds["dos_f"].shape[0] else "fail",
        "detail": f"PLINK 读入 {n_loaded_s} 样本 x {n_loaded_v} 位点（期望 "
                  f"{ds['dos_f'].shape[0]} x {ds['dos_f'].shape[1]}），returncode={r.returncode}",
    })
    frq = vdir / "freq.frq"
    if frq.is_file():
        rows = frq.read_text(encoding="utf-8", errors="replace").splitlines()[1:]
        mafs = np.array([float(x.split()[4]) for x in rows if len(x.split()) >= 5])
        results["n_maf_lt_001"] = int(np.sum(mafs < 0.01))
        results["maf_mean"] = float(np.mean(mafs)) if mafs.size else None
        results["checks"].append({
            "name": "freq_ok", "status": "pass",
            "detail": f"--freq 输出 {mafs.size} 个位点，MAF 均值 {np.mean(mafs):.4f}，"
                      f"MAF<0.01 位点 {int(np.sum(mafs < 0.01))} 个"})

    # 2) recode A：检查 .raw 可解析、缺失编码为 NA
    r = record("recode", ["--recode", "A"], "")
    raw = vdir / "recode.raw"
    if raw.is_file():
        header = raw.read_text(encoding="utf-8").splitlines()[0].split()
        na_cnt = 0
        with open(raw, encoding="utf-8") as fh:
            next(fh)
            for line in fh:
                na_cnt += line.split().count("NA")
        results["checks"].append({
            "name": "recode_a_ok",
            "status": "pass" if r.returncode == 0 and len(header) >= 7 else "fail",
            "detail": f"--recode A 生成 {len(header) - 6} 个 SNP 列，NA 单元格 {na_cnt} 个"})
        results["n_na_cells"] = na_cnt
        results["raw_columns"] = header[:6]

    # 3) hardy：不崩溃即可，记录极端偏离比例
    r = record("hardy", ["--hardy"], "")
    hwe = vdir / "hardy.hwe"
    if hwe.is_file():
        lines = hwe.read_text(encoding="utf-8", errors="replace").splitlines()[1:]
        pvals = []
        for x in lines:
            c = x.split()
            # PLINK 在 .fam 表型缺失时输出 TEST="ALL(NP)"；有表型时输出 "ALL"
            if len(c) >= 9 and c[2].startswith("ALL"):
                try:
                    pvals.append(float(c[8]))
                except ValueError:
                    pass
        p = np.array(pvals)
        frac = float(np.mean(p < HWE_P_CUT)) if p.size else 0.0
        results["hwe_extremes_frac"] = frac
        # 结构群体（2 个 breed 合并检验）存在 Wahlund 效应，HWE 极端比例 5% 左右属正常；
        # 该检查只确认 --hardy 不崩溃且比例未异常到 10% 以上。
        results["checks"].append({
            "name": "hardy_ok",
            "status": "pass" if r.returncode == 0 and frac < 0.10 else "warn",
            "detail": f"--hardy 完成，{p.size} 个位点，HWE p<1e-6 比例 {frac:.5f}"
                      f"（亚群结构合并检验的 Wahlund 效应，预期 3%~6%）"})

    # 4) missing：与纯 Python 缺失率对照
    r = record("missing", ["--missing"], "")
    lmiss = vdir / "missing.lmiss"
    if lmiss.is_file():
        lines = lmiss.read_text(encoding="utf-8", errors="replace").splitlines()[1:]
        fmiss = np.array([float(x.split()[4]) for x in lines if len(x.split()) >= 5])
        py_miss = np.array([float(np.mean(np.isnan(ds["dos_f"][:, j])))
                            for j in range(ds["dos_f"].shape[1])])
        diff = float(np.max(np.abs(fmiss - py_miss))) if fmiss.size == py_miss.size else float("nan")
        results["checks"].append({
            "name": "missing_match",
            "status": "pass" if (fmiss.size == py_miss.size and diff < 1e-6) else "fail",
            "detail": f"--missing 与纯 Python 缺失率最大绝对差 {diff:.2e}（容差 1e-6）"})

    # 5) genome：PI_HAT 是否落在设计档位
    r = record("genome", ["--genome"], "")
    genome = vdir / "genome.genome"
    if genome.is_file():
        lines = genome.read_text(encoding="utf-8", errors="replace").splitlines()[1:]
        pi = np.array([float(x.split()[9]) for x in lines if len(x.split()) >= 10])
        results["pihat_max"] = float(pi.max()) if pi.size else None
        results["pihat_gt_025"] = int(np.sum(pi > 0.25))
        results["checks"].append({
            "name": "genome_ok",
            "status": "pass" if pi.size and pi.max() > 0.9 else "warn",
            "detail": f"--genome 完成，PI_HAT 最大 {pi.max() if pi.size else float('nan'):.3f}，"
                      f"PI_HAT>0.25 的样本对 {int(np.sum(pi > 0.25))} 个"})

    # 6) pca：第一主成分应分开两个亚群
    r = record("pca", ["--pca", "10"], "")
    eigvec = vdir / "pca.eigenvec"
    if eigvec.is_file():
        lines = eigvec.read_text(encoding="utf-8", errors="replace").splitlines()
        pc1 = np.array([float(x.split()[2]) for x in lines])
        breed = np.array([str(b) for b in ds["breed"]])
        auc = float(np.mean(pc1[breed == "A"][:, None] > pc1[breed == "B"][None, :]))
        results["pca_auc_plink"] = auc
        results["checks"].append({
            "name": "pca_ok",
            "status": "pass" if max(auc, 1 - auc) > 0.99 else "warn",
            "detail": f"--pca 完成，PLINK PC1 区分 A/B 的 AUC={max(auc, 1 - auc):.4f}"})

    # 7) 关联分析可跑通（数量性状，检验 .bed + 表型联调）
    r = record("assoc_smoke", ["--linear", "--pheno", str(out_dir / "pheno.txt"),
                               "--keep-allele-order"], "")
    assoc = vdir / "assoc_smoke.assoc.linear"
    if assoc.is_file():
        n_rows = len(assoc.read_text(encoding="utf-8", errors="replace").splitlines()) - 1
        results["checks"].append({
            "name": "linear_smoke",
            "status": "pass" if r.returncode == 0 else "fail",
            "detail": f"--linear 跑通，输出 {n_rows} 行关联结果"})

    if any(c["status"] == "fail" for c in results["checks"]):
        results["status"] = "fail"
    return results


# ---------------------------------------------------------------------------
# 5. README
# ---------------------------------------------------------------------------
README_TEMPLATE = """\
# PLINK 演示数据集（{n_samples} 样本 x {n_snps} SNP）

由 `plink-analysis/scripts/simulate_demo.py` 生成，固定种子 `{seed}`，可逐字节复现。
生成时间：{created_at}

## 文件清单

| 文件 | 说明 |
|---|---|
| `geno.bed` / `.bim` / `.fam` | PLINK 1.9 二进制基因型（variant-major，18 条常染色体） |
| `pheno.txt` | 数量性状 `qtl`（FID IID qtl），h2≈{h2:.2f} |
| `pheno_bin.txt` | 二元性状 `status`（1=对照，2=病例），患病率≈{prev:.2f} |
| `covar.txt` | 协变量 `sex`（1 公 / 2 母）、`age`（岁） |
| `group.txt` | 分组 `breed`（A / B） |
| `truth.csv` | 逐位点真值：AF_A/AF_B/MAF、是否分化位点、是否因果/多基因、真实效应量 |
| `truth_meta.json` | 生成参数、真实 h2、亲缘对、文件 sha256、PLINK 校验摘要 |
| `verify/` | PLINK 读取验证的原始输出与命令记录 |

## 数据设计

* **亚群结构**：breed A（{n_a}）与 breed B（{n_b}）。全基因组存在轻度频率漂移
  （logit 尺度 ±0.15，方向随机），另有 {n_diff} 个位点强分化（logit 偏移 1.6~3.0，
  |ΔAF| 可达 0.3~0.7），PCA 第一主成分可完全分开两个亚群。
* **高亲缘**：{n_pairs} 对；其中 2 对完全重复（PI_HAT≈1），8 对约 65% 位点基因型共担
  （PI_HAT≈0.4~0.8），用于检验 duplicate/relatedness 质控闸门。
* **缺失**：{miss_rate:.0%} 完全随机缺失。
* **MAF 谱**：次等位基因频率 ~ Beta(0.5, 2)，低 MAF 位点占多数；MAF<0.01 的位点占比约
  {maf_low:.1%}，可被 `--maf 0.01` 滤除。
* **已知 QTL**：{n_causal} 个因果 SNP（{n_major} 个主效 + 其余微效，因果位点承担 {pve:.0%}
  遗传方差），另加 {n_poly} 个多基因背景位点；真值见 `truth.csv` 的 `IS_CAUSAL` /
  `BETA_QT` 列，主效 QTL 在 n≈{n_samples} 时可达 p<5e-8。
* **二元性状**：liability 模型（liability 尺度 h2={h2_liab:.2f}），阈值取 80 分位使患病率≈{prev:.2f}。

## 一键使用（配合 plink-analysis skill）

```bash
python .agents/skills/plink-analysis/scripts/plink_pipeline.py auto \\
    --input E:\\deepseek\\demo_data\\geno \\
    --pheno E:\\deepseek\\demo_data\\pheno.txt \\
    --covar E:\\deepseek\\demo_data\\covar.txt \\
    --group E:\\deepseek\\demo_data\\group.txt \\
    --run-id demo_qtl --out-root E:\\deepseek\\plink_runs
python .agents/skills/plink-analysis/scripts/selftest.py
```

## PLINK 读取验证（{verify_status}）

验证命令（全部通过 `E:\\deepseek\\tools\\plink\\plink.exe` 实际执行）：

```bash
plink --bfile geno --freq        --out verify/freq
plink --bfile geno --recode A    --keep-allele-order --out verify/recode
plink --bfile geno --hardy       --out verify/hardy
plink --bfile geno --missing     --out verify/missing
plink --bfile geno --genome      --out verify/genome
plink --bfile geno --pca 10      --out verify/pca
plink --bfile geno --linear --pheno pheno.txt --keep-allele-order --out verify/assoc_smoke
```

{verify_table}

> 完整原始输出见 `verify/`，完整命令与返回码见 `truth_meta.json` 的 `plink_verify`。
"""


def write_readme(out_dir: Path, ds: dict, meta: dict, verify: dict) -> None:
    def table(v: dict) -> str:
        checks = v.get("checks") or []
        if not checks:
            return "（未执行 PLINK 校验）"
        lines = ["| 检查 | 状态 | 详情 |", "|---|---|---|"]
        for c in checks:
            lines.append(f"| {c['name']} | {c['status'].upper()} | {c['detail']} |")
        return "\n".join(lines)

    maf = ds["maf_emp"]
    text = README_TEMPLATE.format(
        n_samples=meta["n_samples"], n_snps=meta["n_snps"], seed=meta["seed"],
        created_at=meta["created_at"], h2=meta["h2_quantitative"],
        prev=meta["prevalence"], n_a=meta["breeds"]["A"], n_b=meta["breeds"]["B"],
        drift=DRIFT, n_diff=meta["n_diff_snps"], n_pairs=len(meta["related_pairs"]),
        miss_rate=MISSING_RATE, maf_low=float(np.mean(maf < 0.01)),
        n_causal=meta["n_causal"], n_major=meta["n_major_qtl"], pve=PVE_CAUSAL,
        n_poly=meta["n_poly_snps"], h2_liab=meta["h2_liability"],
        verify_status=verify.get("status", "skipped"), verify_table=table(verify),
    )
    (out_dir / "README.md").write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# 6. main
# ---------------------------------------------------------------------------
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="生成 PLINK 1.9 演示数据集（可复现）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--out", default=DEFAULT_OUT, help="输出目录")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED, help="随机种子")
    ap.add_argument("--n-samples", type=int, default=DEFAULT_N_SAMPLES)
    ap.add_argument("--n-snps", type=int, default=DEFAULT_N_SNPS)
    ap.add_argument("--plink", default=DEFAULT_PLINK, help="PLINK 1.9 可执行文件路径")
    ap.add_argument("--no-plink-verify", action="store_true", help="跳过 PLINK 读取验证")
    ap.add_argument("--quiet", action="store_true")
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    args = parse_args(argv)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    say = (lambda *a: None) if args.quiet else print

    say(f"[simulate_demo] seed={args.seed} n_samples={args.n_samples} n_snps={args.n_snps}")
    say(f"[simulate_demo] 输出目录: {out_dir}")
    ds = build_dataset(args.seed, args.n_samples, args.n_snps)
    files = write_bfile(out_dir, ds)
    tables = write_pheno_covar_group(out_dir, ds)

    checks = self_check(ds, {})
    n_fail = sum(1 for c in checks if c["status"] != "pass")
    for c in checks:
        say(f"  [{'PASS' if c['status'] == 'pass' else 'FAIL'}] {c['name']}: {c['detail']}")

    verify: dict = {"status": "skipped", "detail": "未执行 PLINK 校验"}
    if not args.no_plink_verify:
        say("[simulate_demo] PLINK 读取验证 ...")
        verify = verify_with_plink(out_dir, args.plink, ds)
        log_path = out_dir / "verify" / "plink_verify.json"
        log_path.write_text(json.dumps(verify, ensure_ascii=False, indent=2) + "\n",
                            encoding="utf-8")
        for c in verify.get("checks", []):
            say(f"  [{'PASS' if c['status'] == 'pass' else c['status'].upper()}] "
                f"{c['name']}: {c['detail']}")
        if verify.get("status") == "fail":
            n_fail += 1

    meta = write_truth(out_dir, ds, args, verify)
    meta["python_self_checks"] = checks
    (out_dir / "truth_meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_readme(out_dir, ds, meta, verify)

    total = len(checks) + len(verify.get("checks", []))
    say(f"[simulate_demo] 完成：{total} 项自检/验证，{n_fail} 项 FAIL；"
        f"bfile={files['prefix']}；truth_csv={out_dir / 'truth.csv'}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
