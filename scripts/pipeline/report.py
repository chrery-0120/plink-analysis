"""Report generation for the ``plink-analysis`` skill (DESIGN.md section 9).

Public API
----------
``stage_report(run_dir, stage, metrics, config, summary_override=None) -> Path``
    Write ``<run_dir>/<stage_dir>/report.md`` (Chinese, fixed section structure).
``final_report(run_dir, config, state) -> dict``
    Read every stage ``metrics.json``/``report.md`` and write
    ``05_reports/final_report.md``, ``05_reports/final_report.html``
    (self-contained: inline CSS + base64 images) and ``05_reports/summary.json``.
    Returns ``{"md": <str path>, "html": <str path>, "summary": <str path>}``.
``render_markdown_to_html(md_text, out_path, *, title=..., base_dir=None) -> Path``
    Small Markdown -> HTML renderer used by the final report.

Everything is defensive: missing metrics fields, empty tables, missing figures and
skipped stages never raise; the report is still produced and clearly marks the
stage as not executed.  All text is UTF-8; numbers use scientific notation where
appropriate; the HTML is intended to be opened directly in a browser or printed
to PDF from the browser.
"""

from __future__ import annotations

import base64
import csv
import html as _html
import json
import math
import mimetypes
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

try:  # pandas is part of the documented environment; keep report importable without it
    import pandas as pd
except Exception:  # pragma: no cover
    pd = None


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
STAGE_DIRS: dict[str, str] = {
    "intake": "00_intake",
    "qc": "01_qc",
    "popgen": "02_popgen",
    "gwas": "03_gwas",
    "calibration": "04_calibration",
}
STAGE_ORDER: list[str] = ["intake", "qc", "popgen", "gwas", "calibration"]
REPORTS_DIR = "05_reports"

_STATUS_BADGES = {
    "pass": "✅ PASS",
    "warn": "⚠️ WARN",
    "fail": "❌ FAIL",
    "skipped": "⏭️ SKIPPED",
    "running": "⏳ RUNNING",
    "unknown": "❔ UNKNOWN",
}
_STATUS_RANK = {"pass": 0, "skipped": 0, "unknown": 2, "running": 2, "warn": 1, "fail": 3}

_SPECIES_ZH = {
    "auto": "自动推断",
    "pig": "猪",
    "sus_scrofa": "猪（Sus scrofa）",
    "sus_scrofa_domesticus": "家猪（Sus scrofa domesticus）",
    "cattle": "牛",
    "bos_taurus": "牛（Bos taurus）",
    "sheep": "绵羊",
    "ovis_aries": "绵羊（Ovis aries）",
    "goat": "山羊",
    "capra_hircus": "山羊（Capra hircus）",
    "chicken": "鸡",
    "gallus_gallus": "鸡（Gallus gallus）",
    "horse": "马",
    "equus_caballus": "马（Equus caballus）",
    "human": "人",
    "homo_sapiens": "人（Homo sapiens）",
    "other": "其他",
}

# ---------------------------------------------------------------------------
# Stage metadata
# ---------------------------------------------------------------------------
_STAGE_META: dict[str, dict[str, Any]] = {
    "intake": {
        "name": "数据接入与格式识别",
        "dir": "00_intake",
        "purpose": (
            "将用户提供的基因型、表型、协变量与分组数据接入分析流程，识别文件格式"
            "（bed/ped/vcf/gen/23andMe）、染色体命名体系与物种，核对样本/位点计数与文件完整性，"
            "并在进入质量控制之前暴露编码、ID 与表型缺失值问题。数据接入是保证后续所有统计结果可解释的第一道关口："
            "若格式或 ID 在此阶段静默出错，后面的质控与关联分析都会建立在不完整的数据之上。"
        ),
        "methods": [
            "使用 PLINK 1.9 读取并校验基因型文件，识别 bed/ped/vcf/gen/23andMe 等格式；"
            "对 VCF 等多等位格式按需拆分或仅保留双等位位点。",
            "统计样本数、位点数与染色体分布，核对 .fam/.bim/.bed 三者 ID 一致性与文件大小。",
            "解析表型编码（-9/NA/0 等缺失码），判断数量性状或二元性状；检查协变量与分组文件的列名、"
            "缺失与类别水平。",
            "依据染色体集合与标签推断物种（参考：猪 18 对、牛 29 对、绵羊 26 对、山羊 29 对、鸡 38 对常染色体），"
            "必要时以 --chr-set 覆盖推断结果。",
            "预检记录文件存在性、可读性、SHA-256 与格式警告；任何无法解析的输入必须在进入质控前暴露。",
        ],
        "param_sections": [("输入数据", "input"), ("物种设置", "species")],
        "handling": [
            "确认输入路径与文件权限；Windows 中文路径建议加引号并检查控制台编码。",
            "ID 不一致时按 references/05-error-playbook.md 使用 --double-id / --const-fid 或修正表型文件。",
        ],
    },
    "qc": {
        "name": "质量控制",
        "dir": "01_qc",
        "purpose": (
            "去除低质量样本与位点，降低基因型错误、缺失、偏分离与群体分层对下游 PCA、亲缘推断和 GWAS 的干扰。"
            "质控标准必须在分析开始前固定并记录（见 run.json 的 qc 段），避免根据结果反复调整阈值造成"
            "结果导向的过滤（p-hacking）。本阶段同时输出样本/位点保留率与各项诊断分布，为后续解释提供依据。"
        ),
        "methods": [
            "样本层面：--mind 过滤高缺失样本；--het 计算杂合度并剔除 ±k 倍标准差异常样本；"
            "--genome 估计 PI_HAT 识别重复与近亲样本；--check-sex 校验性别一致性；基于 PCA 识别离群样本。",
            "位点层面：--geno 缺失率过滤；--maf 最小等位基因频率过滤；--hwe 哈迪-温伯格平衡检验过滤；"
            "保留双等位位点并按需去除单态位点。",
            "使用 --indep-pairwise 进行 LD 剪枝，生成独立位点集供 PCA、亲缘与校准复用。",
            "所有阈值来自 run.json 且记录在报告中；删除数量与删除原因逐项统计，保证可追溯。",
        ],
        "param_sections": [("质控过滤参数", "qc")],
        "handling": [
            "如阈值过严导致样本/位点大量丢失，可在 run.json 放宽 mind/geno/maf 后使用 --force 重跑。",
            "杂合度整体偏移可能提示样本污染、近交或芯片批次问题，建议结合 PCA 与亲缘结果复核。",
        ],
    },
    "popgen": {
        "name": "群体结构与亲缘分析",
        "dir": "02_popgen",
        "purpose": (
            "评估样本的群体分层、亲缘关系、连锁不平衡与近交水平，为 GWAS 协变量（PCA 主成分）选择、"
            "独立性假设检查以及保种/育种决策提供依据。群体结构未被校正会同时造成假阳性和功效损失，"
            "因此本阶段的结果直接决定 GWAS 模型是否可信。"
        ),
        "methods": [
            "基于 LD 剪枝位点运行 --pca，报告各主成分方差解释比例与离群样本。",
            "使用 --genome 估计样本对 PI_HAT，统计重复样本（≥ duplicate_pihat）与近亲样本（≥ related_pihat）。",
            "使用 --r2 估计连锁不平衡，按物理距离分箱汇总平均 r²，评估标记密度与 LD 剪枝参数是否合理。",
            "使用 --homozyg 检测纯合片段（ROH），计算每个样本的 FROH 作为近交程度的基因组代理指标。",
            "使用 --fst 估计群体间遗传分化系数，识别分化位点与品种间差异。",
        ],
        "param_sections": [("群体结构参数", "popgen")],
        "handling": [
            "若 PCA 显示明显分层，需在 GWAS 中纳入足够数量的主成分或分群体分析。",
            "若近亲样本对比例高，建议剔除或改用可处理亲缘结构的混合模型工具。",
        ],
    },
    "gwas": {
        "name": "全基因组关联分析",
        "dir": "03_gwas",
        "purpose": (
            "在控制群体结构与固定效应的前提下，检验全基因组标记与表型的统计关联，给出效应量、标准误、"
            "显著性阈值与 lead SNP，为候选基因定位、标记辅助选择与后续功能验证提供统计证据。"
        ),
        "methods": [
            "依据表型类型选择 --linear（数量性状）或 --logistic（二元性状）；协变量包含 PCA 主成分"
            "以及用户提供的性别、年龄、批次等固定效应。",
            "输出 beta、SE、OR（二元性状）与 p 值；采用 Bonferroni、FDR（BH）与建议性阈值三档判定，"
            "并记录阈值来源。",
            "以 λGC 与 QQ 图评估 p 值校准；λGC 偏离 1 提示群体分层残余、亲缘结构或表型转换问题。",
            "使用 --clump 对显著位点做 LD 聚类并报告 lead SNP；通过 LD 与物理位置解释信号区间。",
            "PLINK 1.9 不提供混合模型（LMM），在强亲缘/强分层群体中应结合外部工具复核。",
        ],
        "param_sections": [("关联分析参数", "gwas")],
        "handling": [
            "λGC 明显 > 1.1 时，增加 PCA 主成分、剔除近亲样本或改用混合模型后复核。",
            "显著位点稀少时，先排查表型编码、协变量缺失与样本量，再考虑放宽阈值（须在报告中声明）。",
        ],
    },
    "calibration": {
        "name": "仿真校准与独立校验",
        "dir": "04_calibration",
        "purpose": (
            "用仿真数据与独立重算检验流程本身的统计校准性：零模型下的 λGC 与假阳性率是否受控、"
            "已知效应下的功效是否合理、效应量估计是否无偏，以及关键统计量能否从原始数据独立复现。"
            "本阶段为“结果可信”提供证据，而不只是给出显著性列表。"
        ),
        "methods": [
            "零模型仿真：使用与真实数据无关的表型重复 R 次关联分析，统计 λGC 分布与"
            "α = 0.05 / 1e-3 / 5e-8 的经验假阳性率，期望分别为 0.05 / 0.001 / 5e-8。",
            "已知 QTL 仿真：按设定遗传力 h² 与因果位点数生成表型，估计不同阈值下的统计功效、"
            "效应量偏差与因果位点在结果中的排名（AUC）。",
            "独立重算：从原始基因型/表型重新计算 MAF、缺失率、λGC、beta/SE，并与主流程结果比对"
            "（容差见 validation.json），以发现实现或文件读写错误。",
        ],
        "param_sections": [("仿真校准参数", "calibration"), ("运行时参数", "runtime")],
        "handling": [
            "若零模型 λ 均值偏离 [0.97, 1.03]，检查表型构建、协变量泄漏与主成分数量。",
            "若功效显著低于预期，确认因果位点 MAF、样本量与阈值设置，再评估研究设计是否足够。",
        ],
    },
}

_PARAM_DESC: dict[str, str] = {
    "mind": "样本缺失率上限（PLINK --mind）",
    "geno": "位点缺失率上限（PLINK --geno）",
    "maf": "最小等位基因频率下限（PLINK --maf）",
    "hwe": "HWE 检验 p 值下限（PLINK --hwe）",
    "het_sd": "杂合度离群判定标准差倍数（--het）",
    "duplicate_pihat": "重复样本 PI_HAT 阈值",
    "related_pihat": "近亲样本 PI_HAT 阈值",
    "sex_check": "性别一致性检查模式（--check-sex）",
    "pca_outlier_sd": "PCA 离群判定标准差倍数",
    "remove_related": "是否剔除近亲样本",
    "remove_pca_outliers": "是否剔除 PCA 离群样本",
    "ld_prune.window": "LD 剪枝窗口（SNP 数，--indep-pairwise）",
    "ld_prune.step": "LD 剪枝步长（SNP 数）",
    "ld_prune.r2": "LD 剪枝 r² 阈值",
    "pca": "输出主成分数（--pca）",
    "ld_window_kb": "LD 计算窗口（kb，--r2）",
    "ld_bin_kb": "LD 距离分箱宽度（kb）",
    "roh.snp": "ROH 检测最少 SNP 数（--homozyg）",
    "roh.kb": "ROH 检测最小长度（kb）",
    "roh.density": "ROH 区域 SNP 密度阈值（每 kb）",
    "roh.gap": "ROH 允许的最大间隔（kb）",
    "roh.het_missing": "ROH 内允许的杂合/缺失比例",
    "model": "关联模型：auto/linear/logistic",
    "covariates": "协变量设置（auto 表示自动选择）",
    "n_pcs": "纳入 GWAS 的 PCA 主成分数",
    "ci": "效应量置信区间水平",
    "clump.p1": "clump 显著阈值 p1",
    "clump.p2": "clump 次级阈值 p2",
    "clump.r2": "clump 的 r² 聚类阈值",
    "clump.kb": "clump 的物理距离窗口（kb）",
    "genomic_control": "基因组控制（none / beta / 1）",
    "enabled": "是否启用仿真校准",
    "n_replicates": "仿真重复次数",
    "n_causal": "因果位点数",
    "h2": "设定遗传力 h²",
    "n_snps_subset": "仿真使用的 SNP 子集大小",
    "alpha": "显著性阈值 α",
    "threads": "PLINK 线程数",
    "memory_mb": "PLINK 内存上限（MB）",
    "timeout_s": "单步命令超时（秒）",
    "keep_intermediate": "是否保留中间文件",
    "genotype.type": "基因型数据格式",
    "genotype.prefix": "PLINK 前缀（不含扩展名）",
    "genotype.path": "基因型文件/目录路径",
    "phenotype.path": "表型文件路径",
    "phenotype.type": "表型类型（quantitative/binary/auto）",
    "phenotype.missing": "表型缺失值编码",
    "covariates.path": "协变量文件路径",
    "covariates.columns": "使用的协变量列",
    "covariates.categorical": "按类别处理的协变量",
    "groups.path": "分组文件路径",
    "groups.name": "分组变量名（如品种）",
    "pedigree.path": "系谱文件路径",
    "name": "物种名称",
    "chr_set": "染色体组设置（--chr-set）",
    "autosome_count": "常染色体数目",
    "sex_system": "性别决定系统",
    "nonautosomal": "非常染色体标签",
}

_METRIC_LABELS: dict[str, str] = {
    "genotype_type": "基因型格式",
    "n_samples": "样本数",
    "n_variants": "位点数",
    "chromosomes": "染色体集合",
    "species_inferred": "物种推断结果",
    "files": "输入文件",
    "preflight": "预检",
    "phenotype": "表型",
    "covariates": "协变量",
    "groups": "分组",
    "n_samples_in": "质控前样本数",
    "n_samples_out": "质控后样本数",
    "n_variants_in": "质控前位点数",
    "n_variants_out": "质控后位点数",
    "removed": "剔除数量",
    "missing_rate": "缺失率",
    "maf": "MAF",
    "hwe": "HWE",
    "het": "杂合度",
    "sex": "性别检查",
    "relatedness": "亲缘关系",
    "pca": "PCA",
    "pruned": "LD 剪枝",
    "gates": "质控判定",
    "overall": "总体判定",
    "n_pcs": "主成分数",
    "var_explained": "方差解释比例",
    "n_outliers": "离群样本数",
    "outliers": "离群样本",
    "n_pairs_gt_025": "PI_HAT > 0.25 的样本对数",
    "max": "最大值",
    "min": "最小值",
    "mean": "均值",
    "sd": "标准差",
    "n_roh": "ROH 数量",
    "mean_roh_per_indiv": "平均每样本 ROH 数",
    "mean_froh": "平均 FROH",
    "top_samples": "FROH 最高的样本",
    "fst": "FST",
    "global": "全局值",
    "top": "极值位点",
    "diversity": "遗传多样性",
    "mean_het": "平均杂合度",
    "n_groups": "分组数",
    "model": "关联模型",
    "phenotype_type": "表型类型",
    "n_tested": "检验位点数",
    "lambda_gc": "λGC",
    "lambda_source": "λGC 计算来源",
    "thresholds": "显著性阈值",
    "n_significant": "Bonferroni 显著位点数",
    "n_fdr05": "FDR < 0.05 位点数",
    "n_suggestive": "建议性阈值位点数",
    "lead_snps": "Lead SNP 列表",
    "inflation": "p 值膨胀",
    "null": "零模型",
    "n_replicates": "重复次数",
    "lambda_mean": "λGC 均值",
    "lambda_sd": "λGC 标准差",
    "type1": "I 类错误率",
    "power": "统计功效",
    "effect_bias": "效应量偏差",
    "concordance": "一致性",
    "validation": "独立校验",
    "n_checks": "校验项数",
    "n_pass": "通过项数",
    "n_warn": "警告项数",
    "n_fail": "失败项数",
    "checks": "校验明细",
    "sample_mean": "样本平均缺失率",
    "variant_mean": "位点平均缺失率",
    "n_below_001": "MAF < 0.01 的位点数",
    "n_p_lt_1e-6": "HWE p < 1e-6 的位点数",
    "lambda": "λ",
    "checked": "是否检查",
    "discordant": "性别不一致样本",
    "n_duplicates": "重复样本对数",
    "n_related": "近亲样本对数",
    "max_pihat": "最大 PI_HAT",
    "mean_pihat": "平均 PI_HAT",
    "eigenvalues": "特征值",
    "n_variants_out": "质控后位点数",
    "n_pairs": "位点对数",
    "bins": "LD 距离分箱",
    "distance_bp": "距离（bp）",
    "mean_r2": "平均 r²",
    "n_roh": "ROH 数量",
    "froh": "FROH",
    "total_kb": "ROH 总长度（kb）",
    "snp": "位点",
    "iid": "样本 IID",
    "chr": "染色体",
    "bp": "物理位置",
    "p": "P 值",
    "p_distribution.min": "最小 P 值",
    "p_distribution.median": "P 值中位数",
    "p_distribution.n_p_lt_5e8": "P < 5e-8 的位点数",
    "beta": "效应量 beta",
    "se": "标准误",
    "af": "等位基因频率",
    "a1": "A1",
    "a2": "A2",
    "alpha_5e-8": "α = 5e-8",
    "alpha_1e-4": "α = 1e-4",
    "alpha_0.05": "α = 0.05",
    "mean_beta_ratio": "平均 beta 比值",
    "se_calibration": "SE 校准系数",
    "top_causal_in_top100": "因果位点进入前 100 比例",
    "rank_auc": "排名 AUC",
    # --- nested keys that would otherwise lose context --------------------
    "removed.samples.missingness": "剔除样本 · 缺失率",
    "removed.samples.het_outlier": "剔除样本 · 杂合度离群",
    "removed.samples.sex_discordant": "剔除样本 · 性别不一致",
    "removed.samples.duplicate": "剔除样本 · 重复样本",
    "removed.samples.related": "剔除样本 · 近亲样本",
    "removed.samples.pca_outlier": "剔除样本 · PCA 离群",
    "removed.variants.missingness": "剔除位点 · 缺失率",
    "removed.variants.maf": "剔除位点 · MAF 过低",
    "removed.variants.hwe": "剔除位点 · HWE 偏分离",
    "maf.mean": "平均 MAF",
    "maf.n_below_001": "MAF < 0.01 的位点数",
    "hwe.n_p_lt_1e-6": "HWE p < 1e-6 的位点数",
    "hwe.lambda": "HWE λ",
    "het.mean": "杂合度均值 F",
    "het.sd": "杂合度标准差",
    "sex.checked": "是否检查性别",
    "sex.thresholds.F": "性别判定阈值 F（雌性）",
    "sex.thresholds.M": "性别判定阈值 M（雄性）",
    "pca.eigenvalues": "PCA 特征值",
    "pca.var_explained": "PCA 方差解释比例",
    "pca.n_outliers": "PCA 离群样本数",
    "pruned.n_variants": "LD 剪枝后位点数",
    "preflight.status": "输入预检总体状态",
    "species_inferred.name": "推断物种名称",
    "species_inferred.chr_set": "染色体组设置",
    "species_inferred.sex_system": "性别决定系统",
    "species_inferred.autosomes": "常染色体数",
    "phenotype.present": "是否提供表型",
    "phenotype.type": "表型类型",
    "phenotype.n": "表型样本数",
    "phenotype.missing": "表型缺失数",
    "covariates.present": "是否提供协变量",
    "covariates.columns": "协变量列",
    "groups.present": "是否提供分组",
    "groups.levels": "分组水平",
    "inflation.mean_chi2": "平均 χ²",
    "inflation.lambda": "膨胀因子 λ",
    "thresholds.bonferroni": "Bonferroni 阈值",
    "thresholds.fdr_0.05": "FDR < 0.05 阈值",
    "thresholds.suggestive": "建议性阈值",
    "null.n_replicates": "零模型重复次数",
    "null.lambda_mean": "零模型 λGC 均值",
    "null.lambda_sd": "零模型 λGC 标准差",
    "null.type1.alpha_0.05": "零模型 α=0.05 假阳性率",
    "null.type1.alpha_1e-3": "零模型 α=1e-3 假阳性率",
    "null.type1.alpha_5e-8": "零模型 α=5e-8 假阳性率",
    "power.h2": "设定遗传力 h²",
    "power.n_causal": "因果位点数",
    "concordance.top_causal_in_top100": "因果位点进入前 100 比例",
    "concordance.rank_auc": "因果位点排名 AUC",
}

_TABLE_LABELS: dict[str, str] = {
    # exact table stems (highest priority)
    "sample_qc": "样本水平质控表",
    "qc_summary": "质控汇总表",
    "lead_snps": "Lead SNP 列表",
    "ld_bins": "LD 衰减分箱",
    "validation_checks": "独立校验明细",
    "preflight": "输入预检表",
    "pca": "PCA 结果表",
    # keyword fallbacks for arbitrary stage table names
    "lead": "Lead SNP 列表",
    "qc": "质控汇总表",
    "missing": "缺失率汇总",
    "maf": "MAF 汇总",
    "hwe": "HWE 汇总",
    "het": "杂合度汇总",
    "related": "亲缘关系汇总",
    "kinship": "亲缘矩阵摘要",
    "ld": "LD 衰减分箱",
    "roh": "ROH 汇总",
    "fst": "FST 汇总",
    "validation": "独立校验明细",
    "calibration": "校准结果",
    "power": "功效结果",
    "effect": "效应量偏差",
    "null": "零模型结果",
    "summary": "汇总",
    "top_samples": "FROH 最高的样本",
    "outliers": "离群样本",
    "top": "极值位点",
    "bins": "LD 分箱",
    "checks": "校验明细",
    "levels": "分组水平",
    "files": "输入文件",
}

_FIGURE_CAPTIONS: dict[str, str] = {
    "manhattan": "全基因组曼哈顿图（-log10 P，含显著性阈值）",
    "fst_manhattan": "全基因组 FST 曼哈顿图",
    "hwe_qq": "HWE p 值 QQ 图（质控诊断）",
    "qq": "GWAS p 值 QQ 图（含 λGC 标注）",
    "pca_scatter": "PCA 散点图（PC1 × PC2）",
    "kin": "亲缘关系热图",
    "related": "亲缘关系（PI_HAT）分布",
    "ld": "LD 衰减曲线",
    "roh": "样本 FROH 分布",
    "fst": "FST 分布",
    "missing": "样本/位点缺失率分布",
    "maf": "MAF 频谱",
    "het": "杂合度（F）分布",
    "calibration": "零模型 λGC 校准",
    "power": "统计功效曲线",
    "effect": "效应量校准",
    "pheno": "表型分布",
    "scree": "PCA 碎石图",
    "kinship": "亲缘关系热图",
}


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------
def _now_str() -> str:
    return datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z")


def _read_json(path: Path | str, default: Any = None) -> Any:
    try:
        with Path(path).open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return default


def _write_text(path: Path | str, text: str) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    return out


def _write_json(path: Path | str, obj: Any) -> Path:
    return _write_text(path, json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def _rel_path(path: Path | str, base: Path | str) -> str:
    try:
        return Path(os.path.relpath(str(path), str(base))).as_posix()
    except Exception:  # pragma: no cover
        return str(path)


def _human_size(size: int | float) -> str:
    try:
        value = float(size)
    except Exception:
        return "—"
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024.0 or unit == "GB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024.0
    return f"{value:.1f} GB"


def _fmt(value: Any, *, digits: int = 4) -> str:
    """Format any scalar/list/dict value for a Chinese report table."""
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "是" if value else "否"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "—"
        if value == 0:
            return "0"
        av = abs(value)
        if av < 1e-3 or av >= 1e6:
            return f"{value:.3e}"
        if av >= 100:
            return f"{value:,.2f}".rstrip("0").rstrip(".")
        return f"{value:.{digits}g}"
    if isinstance(value, Mapping):
        try:
            return json.dumps(value, ensure_ascii=False)
        except Exception:
            return str(value)
    if isinstance(value, (list, tuple, set)):
        items = list(value)
        if len(items) <= 10:
            return ", ".join(_fmt(v, digits=digits) for v in items)
        return ", ".join(_fmt(v, digits=digits) for v in items[:10]) + f", …（共 {len(items)} 项）"
    text = str(value).strip()
    return text if text else "—"


def _clean_cell(value: Any, *, max_len: int = 80) -> str:
    text = _fmt(value)
    if text in {"", "nan", "None", "NaN", "NA"}:
        return "—"
    text = text.replace("|", "\\|").replace("\r", " ").replace("\n", " ")
    if len(text) > max_len:
        text = text[: max_len - 1] + "…"
    return text


def _md_table(headers: Sequence[Any], rows: Sequence[Sequence[Any]], *, aligns: Sequence[str] | None = None) -> list[str]:
    """Build a GitHub-flavoured Markdown table from already cleaned rows."""
    head = [str(h) for h in headers]
    body = [list(map(str, row)) for row in rows]
    if aligns is None:
        aligns = []
        for index in range(len(head)):
            column = [row[index] for row in body if index < len(row) and row[index] not in ("", "—")]
            numeric = bool(column) and all(
                re.fullmatch(r"[-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?", cell.replace(" ", "")) for cell in column
            )
            aligns.append("right" if numeric else "left")
    separators = {"left": ":---", "center": ":---:", "right": "---:"}
    lines = [
        "| " + " | ".join(head) + " |",
        "| " + " | ".join(separators.get(a, ":---") for a in aligns) + " |",
    ]
    for row in body:
        cells = [str(c).replace("|", "\\|").replace("\n", " ") for c in row]
        if len(cells) < len(head):
            cells += ["—"] * (len(head) - len(cells))
        lines.append("| " + " | ".join(cells[: len(head)]) + " |")
    return lines


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value] if value.strip() else []
    if isinstance(value, Mapping):
        return [f"{k}: {_fmt(v)}" for k, v in value.items()]
    if isinstance(value, (list, tuple, set)):
        out: list[str] = []
        for item in value:
            if isinstance(item, Mapping):
                out.append(_fmt(item))
            elif item is not None:
                out.append(str(item))
        return out
    return [str(value)]


def _norm_status(value: Any) -> str:
    if value is None:
        return "unknown"
    text = str(value).strip().lower()
    alias = {
        "pass": "pass",
        "passed": "pass",
        "ok": "pass",
        "success": "pass",
        "warn": "warn",
        "warning": "warn",
        "fail": "fail",
        "failed": "fail",
        "error": "fail",
        "skipped": "skipped",
        "skip": "skipped",
        "running": "running",
        "unknown": "unknown",
        "": "unknown",
    }
    return alias.get(text, text if text in _STATUS_BADGES else "unknown")


def _badge(status: Any) -> str:
    return _STATUS_BADGES.get(_norm_status(status), "❔ UNKNOWN")


def _status_rank(status: Any) -> int:
    return _STATUS_RANK.get(_norm_status(status), 2)


def _overall_status(stages: Mapping[str, Mapping[str, Any]]) -> str:
    statuses = [_norm_status(item.get("status")) for item in stages.values()]
    if not statuses:
        return "unknown"
    active = [s for s in statuses if s not in ("skipped", "unknown")]
    if not active:
        return "skipped" if all(s == "skipped" for s in statuses) else "unknown"
    if any(s == "fail" for s in active):
        return "fail"
    if any(s in ("warn", "running") for s in active):
        return "warn"
    return "pass"


def _get_path(node: Any, dotted: str) -> Any:
    current = node
    for part in str(dotted).split("."):
        if not isinstance(current, Mapping) or part not in current:
            return None
        current = current[part]
    return current


def _iter_kv(node: Any, prefix: str = "") -> Iterable[tuple[str, Any]]:
    if not isinstance(node, Mapping):
        return
    for key, value in node.items():
        path = f"{prefix}{key}"
        if isinstance(value, Mapping):
            yield from _iter_kv(value, path + ".")
        elif isinstance(value, (list, tuple)):
            if not value:
                yield path, "—"
            elif all(isinstance(v, (str, int, float, bool)) or v is None for v in value):
                yield path, ", ".join(str(v) for v in value)
            else:
                yield path, f"{len(value)} 项"
        else:
            yield path, value


def _label_for(key: str, labels: Mapping[str, str]) -> str:
    text = str(key)
    if text in labels:
        return labels[text]
    parts = text.split(".")
    for index in range(1, len(parts)):
        candidate = ".".join(parts[index:])
        if candidate in labels:
            return labels[candidate]
    return labels.get(parts[-1], text)


def _table_label(name: str) -> str:
    lowered = str(name).lower()
    if lowered in _TABLE_LABELS:
        return _TABLE_LABELS[lowered]
    for key in sorted(_TABLE_LABELS, key=len, reverse=True):
        if key in lowered and len(key) >= 3:
            return _TABLE_LABELS[key]
    return name


def _param_desc(path: str) -> str:
    text = str(path)
    if text in _PARAM_DESC:
        return _PARAM_DESC[text]
    parts = text.split(".")
    for index in range(1, len(parts)):
        candidate = ".".join(parts[index:])
        if candidate in _PARAM_DESC:
            return _PARAM_DESC[candidate]
    return _PARAM_DESC.get(parts[-1], "")


def _fmt_lambda(value: Any) -> str:
    try:
        return f"{float(value):.3f}"
    except Exception:
        return _fmt(value)


# ---------------------------------------------------------------------------
# P-value display helpers
# ---------------------------------------------------------------------------
_P_TEXT_FIELDS = ("p_text", "pvalue_text", "p_value_text", "p_text_value")
_TEXT_FIELD_SUFFIX = "_text"
_DETAIL_TEXT_COLUMNS = {"p_text", "p_value_text", "pvalue_text", "min_text", "pmin_text"}


def _p_value_text(value: Any) -> str:
    """Format a numeric p-value for display.

    PLINK writes ``0`` for p-values that underflow double precision; those are
    shown as ``<1e-300`` instead of a misleading ``0`` (or ``p = 0``).
    """
    if value is None:
        return "—"
    try:
        number = float(value)
    except Exception:
        text = str(value).strip()
        return text if text else "—"
    if not math.isfinite(number):
        return "—"
    if number <= 0:
        return "<1e-300"
    if number < 1e-3 or number >= 1e6:
        return f"{number:.3e}"
    return _fmt(number)


def _p_display(item: Any) -> str:
    """Display a p-value for a SNP record, preferring its pre-formatted text."""
    if not isinstance(item, Mapping):
        return _p_value_text(item)
    for key in _P_TEXT_FIELDS:
        text = item.get(key)
        if text is not None and str(text).strip() not in {"", "nan", "None"}:
            return str(text).strip()
    for key in ("p", "P", "pval", "p_value", "pvalue"):
        if key in item:
            return _p_value_text(item.get(key))
    return "—"


def _min_p_display(mapping: Any) -> str | None:
    """Display the minimum p-value of a distribution, preferring ``min_text``."""
    if not isinstance(mapping, Mapping):
        return None
    for key in ("min_text", "p_text", "min_p_text", "pmin_text"):
        text = mapping.get(key)
        if text is not None and str(text).strip() not in {"", "nan", "None"}:
            return str(text).strip()
    for key in ("min", "min_p", "minimum"):
        if key in mapping:
            return _p_value_text(mapping.get(key))
    return None


def _calibration_evidence(item: Mapping[str, Any]) -> tuple[bool, bool, int, dict]:
    """Summarise whether a calibration stage produced usable evidence."""
    metrics = item.get("metrics") if isinstance(item.get("metrics"), Mapping) else {}
    null = metrics.get("null") if isinstance(metrics.get("null"), Mapping) else {}
    power = metrics.get("power") if isinstance(metrics.get("power"), Mapping) else {}
    validation = _merge_validation(item)
    has_null = bool(
        null.get("n_replicates")
        or null.get("lambda_mean") is not None
        or (isinstance(null.get("type1"), Mapping) and null.get("type1"))
    )
    has_power = bool(power)
    try:
        n_checks = int(validation.get("n_checks") or 0)
    except Exception:
        n_checks = 0
    return has_null, has_power, n_checks, validation


def _skip_reason(item: Mapping[str, Any], default: str) -> str:
    """Extract the human-readable reason a stage was skipped/not completed."""
    metrics = item.get("metrics") if isinstance(item.get("metrics"), Mapping) else {}
    state = item.get("state") if isinstance(item.get("state"), Mapping) else {}
    for source in (state, metrics):
        for key in ("skip_reason", "reason", "error", "summary", "note", "message", "detail"):
            value = source.get(key)
            if value is not None and str(value).strip():
                return " ".join(str(value).split())
    return default


def _short_reason(reason: Any, max_chars: int = 160) -> str:
    """Keep a reason concise (first sentence) for inline report notes."""
    text = " ".join(str(reason or "").split())
    if not text:
        return "未提供具体原因"
    for separator in ("。", "；", ";", "！", "?"):
        head = text.split(separator, 1)[0].strip()
        if head and len(head) <= max_chars:
            suffix = separator if separator in {"。", "；", "！"} else ""
            return head + suffix
    if len(text) > max_chars:
        return text[: max_chars - 1] + "…"
    return text


def _collapse_p_text_columns(
    headers: Sequence[str], rows: Sequence[Sequence[str]]
) -> tuple[list[str], list[list[str]]]:
    """Prefer ``P_TEXT``/``MIN_TEXT`` columns over numeric ``P``/``MIN`` in CSV tables."""
    new_headers = [str(h) for h in headers]
    new_rows = [[("" if cell is None else str(cell)) for cell in row] for row in rows]
    lower = [h.strip().lower() for h in new_headers]
    text_idx = None
    for index, name in enumerate(lower):
        if name in {"p_text", "p_value_text", "pvalue_text", "min_text", "pmin_text"}:
            text_idx = index
            break
    numeric_idx = None
    for index, name in enumerate(lower):
        if index == text_idx:
            continue
        if name in {"p", "pval", "pvalue", "p_value", "min", "min_p"}:
            numeric_idx = index
            break
    if text_idx is not None:
        keep = [index for index in range(len(new_headers)) if index != numeric_idx]
        collapsed_headers: list[str] = []
        for index in keep:
            if index == text_idx:
                collapsed_headers.append("最小 P 值" if lower[index].startswith("min") else "P 值")
            else:
                collapsed_headers.append(new_headers[index])
        collapsed_rows: list[list[str]] = []
        for row in new_rows:
            values = []
            for index in keep:
                value = row[index] if index < len(row) else ""
                if (
                    index == text_idx
                    and str(value).strip() in {"", "—", "nan", "None", "NA"}
                    and numeric_idx is not None
                ):
                    value = row[numeric_idx] if numeric_idx < len(row) else ""
                values.append(value)
            collapsed_rows.append(values)
        new_headers, new_rows = collapsed_headers, collapsed_rows
        lower = [h.strip().lower() for h in new_headers]
    # Underflowed numeric p-values (PLINK prints 0) must not be displayed as 0.
    p_columns = [i for i, name in enumerate(lower) if name in {"p", "pval", "pvalue", "p_value"}]
    for row in new_rows:
        for index in p_columns:
            if index < len(row):
                try:
                    if float(row[index]) <= 0:
                        row[index] = "<1e-300"
                except Exception:
                    pass
    return new_headers, new_rows


def _detail_table_columns(items: Sequence[Mapping[str, Any]], max_cols: int = 9) -> list[str]:
    columns: list[str] = []
    for item in items:
        for key in item.keys():
            name = str(key)
            if name.lower() in _DETAIL_TEXT_COLUMNS:
                continue
            if name not in columns:
                columns.append(name)
    has_p_text = any(any(str(k).lower() in _DETAIL_TEXT_COLUMNS for k in item.keys()) for item in items)
    if has_p_text and not any(c in columns for c in ("p", "P", "pvalue", "p_value")):
        columns.append("p")
    return columns[:max_cols]


def _detail_cell(item: Mapping[str, Any], column: str) -> str:
    if column == "p":
        return _p_display(item)
    if column == "q":
        return _p_value_text(item.get("q"))
    if column == "min":
        min_text = _min_p_display(item)
        if min_text is not None:
            return min_text
    return _clean_cell(item.get(column))


def _renumber_figures(markdown: str) -> str:
    """Renumber every ``**图 N**`` label sequentially (per report)."""
    counter = [0]

    def repl(_match: re.Match[str]) -> str:
        counter[0] += 1
        return f"**图 {counter[0]}**"

    return re.sub(r"\*\*图 \d+\*\*", repl, str(markdown))


def _format_phenotype(value: Any) -> str:
    if isinstance(value, Mapping):
        if value.get("present") is False:
            return "未提供"
        parts: list[str] = []
        kind = value.get("type")
        if kind:
            labels = {"quantitative": "数量性状", "binary": "二元性状", "case_control": "病例-对照", "auto": "自动判断"}
            parts.append(labels.get(str(kind).lower(), str(kind)))
        if value.get("n") is not None:
            parts.append(f"n = {_fmt(value.get('n'))}")
        if value.get("missing") is not None:
            parts.append(f"缺失 {_fmt(value.get('missing'))}")
        return "，".join(parts) if parts else "—"
    return _fmt(value)


def _format_covariates(value: Any) -> str:
    if isinstance(value, Mapping):
        if value.get("present") is False:
            return "未提供"
        columns = value.get("columns") or value.get("names")
        if isinstance(columns, (list, tuple)) and columns:
            return "、".join(str(c) for c in columns)
        return "—"
    return _fmt(value)


def _format_groups(value: Any) -> str:
    if isinstance(value, Mapping):
        levels = value.get("levels")
        if isinstance(levels, list) and levels:
            parts = []
            for level in levels:
                if isinstance(level, Mapping):
                    parts.append(f"{_fmt(level.get('name'))}（{_fmt(level.get('n'))}）")
                else:
                    parts.append(_fmt(level))
            return "、".join(parts)
        if value.get("present") is False:
            return "未提供"
    return _fmt(value)


# ---------------------------------------------------------------------------
# Metrics normalisation
# ---------------------------------------------------------------------------
_STAGE_KEYS = {"stage", "status", "summary", "metrics", "tables", "figures", "warnings", "errors", "elapsed_s"}


def _extract_metrics(obj: Any) -> dict:
    """Accept either the stage return dict or the bare metrics dict."""
    if not isinstance(obj, Mapping):
        return {}
    if "metrics" in obj and isinstance(obj.get("metrics"), Mapping) and any(k in obj for k in _STAGE_KEYS):
        inner = dict(obj["metrics"])
        for key in ("status", "summary", "tables", "figures", "warnings", "errors", "stage", "overall"):
            if key in obj and key not in inner:
                inner[key] = obj[key]
        return inner
    return dict(obj)


def _parse_tables(raw: Any, stage_dir: Path) -> dict[str, Path]:
    tables: dict[str, Path] = {}
    if isinstance(raw, Mapping):
        for name, value in raw.items():
            if isinstance(value, Mapping):
                continue
            candidate = Path(str(value))
            if not candidate.is_absolute():
                candidate = stage_dir / candidate
            tables[str(name)] = candidate
    elif isinstance(raw, (list, tuple)):
        for value in raw:
            candidate = Path(str(value))
            if not candidate.is_absolute():
                candidate = stage_dir / candidate
            tables[candidate.stem] = candidate
    if not tables:
        table_dir = stage_dir / "tables"
        if table_dir.is_dir():
            for path in sorted(table_dir.glob("*.csv")):
                tables[path.stem] = path
    return {name: path for name, path in tables.items() if path.exists()}


def _parse_figures(raw: Any, stage_dir: Path) -> list[Path]:
    figures: list[Path] = []
    if isinstance(raw, (list, tuple)):
        for value in raw:
            candidate = Path(str(value))
            if not candidate.is_absolute():
                candidate = stage_dir / candidate
            if candidate.exists() and candidate.suffix.lower() == ".png":
                figures.append(candidate)
    if not figures:
        plot_dir = stage_dir / "plots"
        if plot_dir.is_dir():
            figures = sorted(plot_dir.glob("*.png"))
    seen: set[str] = set()
    unique: list[Path] = []
    for path in figures:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def _figure_caption(path: Path) -> str:
    stem = path.stem.lower()
    specific = {
        "calibration_effect_bias": "效应量校准（β_est / β_true）",
        "calibration_power": "已知 QTL 检出力（power）曲线",
        "calibration_null_lambda": "零模型 λGC 分布（重复模拟）",
        "calibration_null_qq": "零模型 p 值 QQ 图",
    }
    for key in sorted(specific, key=len, reverse=True):
        if key in stem:
            return specific[key]
    for key in sorted(_FIGURE_CAPTIONS, key=len, reverse=True):
        if key in stem:
            return _FIGURE_CAPTIONS[key]
    return path.stem


def _read_csv_rows(path: Path, *, max_rows: int = 20, max_cols: int = 15) -> tuple[list[str], list[list[str]], list[str]]:
    """Read a CSV table and truncate to ``max_rows`` x ``max_cols``."""
    headers: list[str] = []
    body: list[list[str]] = []
    total_rows = 0
    total_cols = 0
    try:
        if pd is not None:
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
            headers = [str(c) for c in frame.columns]
            total_rows, total_cols = frame.shape
            work = frame.iloc[:max_rows, :max_cols]
            for _, row in work.iterrows():
                body.append([_csv_cell(row[col]) for col in work.columns])
        else:  # pragma: no cover - pandas is present in the target environment
            with path.open("r", encoding="utf-8-sig", newline="") as fh:
                reader = list(csv.reader(fh))
            if reader:
                headers = [str(c).strip() for c in reader[0]]
                rows = reader[1:]
            else:
                rows = []
            total_rows, total_cols = len(rows), len(headers)
            for row in rows[:max_rows]:
                body.append([_csv_cell(c) for c in row[:max_cols]])
    except Exception:
        return [], [], [f"无法读取数据表 `{path.name}`（文件损坏或编码异常）。"]
    notes: list[str] = []
    if total_rows > max_rows:
        notes.append(f"原表共 {total_rows:,} 行，此处仅展示前 {max_rows} 行；完整数据见 `{path.name}`。")
    if total_cols > max_cols:
        notes.append(f"原表共 {total_cols:,} 列，此处仅展示前 {max_cols} 列。")
    headers, body = _collapse_p_text_columns(headers, body)
    return headers, body, notes


def _csv_cell(value: Any) -> str:
    if value is None:
        return "—"
    if pd is not None:
        try:
            if pd.isna(value):
                return "—"
        except Exception:
            pass
    text = str(value).strip()
    if text in {"", "nan", "NaN", "NA", "None"}:
        return "—"
    text = text.replace("|", "\\|").replace("\r", " ").replace("\n", " ")
    if len(text) > 80:
        text = text[:79] + "…"
    return text


def _gate_rows(metrics: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    for key in ("gates", "checks"):
        raw = metrics.get(key)
        if isinstance(raw, list) and raw:
            return [item for item in raw if isinstance(item, Mapping)]
    preflight = metrics.get("preflight")
    if isinstance(preflight, Mapping) and isinstance(preflight.get("checks"), list):
        return [item for item in preflight["checks"] if isinstance(item, Mapping)]
    return []


def _render_gates(metrics: Mapping[str, Any]) -> list[str]:
    rows = _gate_rows(metrics)
    lines: list[str] = []
    if not rows:
        status = _norm_status(metrics.get("overall") or metrics.get("status"))
        if status != "unknown":
            lines.append(
                f"本阶段未提供逐项判定表，总体状态为 {_badge(status)}；"
                "详细阈值与实际值见上文“关键指标”。"
            )
        else:
            lines.append("本阶段未提供质控判定信息。")
        return lines
    counts = {"pass": 0, "warn": 0, "fail": 0, "unknown": 0}
    table_rows: list[list[str]] = []
    for item in rows:
        status = _norm_status(item.get("status") or item.get("result"))
        counts[status if status in counts else "unknown"] += 1
        table_rows.append(
            [
                _clean_cell(item.get("name") or item.get("check") or item.get("id") or "—"),
                _badge(status),
                _clean_cell(item.get("value", item.get("observed", "—"))),
                _clean_cell(item.get("threshold", "—")),
                _clean_cell(item.get("detail", item.get("note", "—")), max_len=100),
            ]
        )
    lines.append(
        f"共 {len(rows)} 项判定：✅ PASS {counts['pass']} 项，⚠️ WARN {counts['warn']} 项，"
        f"❌ FAIL {counts['fail']} 项" + (f"，❔ 未知 {counts['unknown']} 项" if counts["unknown"] else "") + "。"
    )
    lines.append("")
    lines.extend(_md_table(["判定项", "状态", "实测值", "阈值", "说明"], table_rows))
    return lines


# ---------------------------------------------------------------------------
# Stage report helpers
# ---------------------------------------------------------------------------
def _discover_extra_metric_tables(metrics: Mapping[str, Any]) -> list[tuple[str, list[Mapping[str, Any]]]]:
    found: list[tuple[str, list[Mapping[str, Any]]]] = []
    skip = {"gates", "checks", "preflight"}

    def walk(node: Any, prefix: str = "") -> None:
        if isinstance(node, Mapping):
            for key, value in node.items():
                if str(key) in skip:
                    continue
                path = f"{prefix}{key}"
                if isinstance(value, list) and value and all(isinstance(v, Mapping) for v in value):
                    found.append((path, value))  # type: ignore[arg-type]
                elif isinstance(value, (Mapping, list)):
                    walk(value, path + ".")
        elif isinstance(node, list):
            for index, value in enumerate(node):
                if isinstance(value, Mapping):
                    walk(value, f"{prefix}[{index}].")

    walk(metrics)
    return found


def _ensure_stage_dir(run_dir: Path, stage: str) -> tuple[Path, dict[str, Any]]:
    meta = _STAGE_META.get(stage)
    if meta is None:
        stage_dir = run_dir / STAGE_DIRS.get(stage, stage)
        meta = {
            "name": stage or "未知阶段",
            "dir": stage_dir.name,
            "purpose": "（未定义的阶段类型）",
            "methods": [],
            "param_sections": [],
            "handling": [],
        }
    else:
        stage_dir = run_dir / meta["dir"]
    stage_dir.mkdir(parents=True, exist_ok=True)
    return stage_dir, meta


def _project_title(config: Mapping[str, Any], run_dir: Path | None = None) -> str:
    try:
        title = _get_path(config, "report.title")
    except Exception:
        title = None
    if title:
        return str(title)
    run_id = config.get("run_id") if isinstance(config, Mapping) else None
    if run_id:
        return str(run_id)
    if run_dir is not None:
        return run_dir.name
    return "PLINK 全流程分析"


def _species_label(config: Mapping[str, Any], metrics: Mapping[str, Any] | None = None) -> str:
    name = None
    if isinstance(config, Mapping):
        name = _get_path(config, "species.name")
    if metrics:
        inferred = metrics.get("species_inferred")
        if isinstance(inferred, Mapping) and inferred.get("name"):
            name = inferred.get("name")
    if not name:
        return "未指定"
    text = str(name)
    return _SPECIES_ZH.get(text.lower(), text)


def _stage_param_lines(stage: str, config: Mapping[str, Any], metrics: Mapping[str, Any]) -> list[str]:
    meta = _STAGE_META.get(stage, {})
    lines: list[str] = []
    for section_title, dotted in meta.get("param_sections", []):
        node = _get_path(config, dotted)
        if not isinstance(node, Mapping) or not node:
            continue
        rows: list[list[str]] = []
        for key, value in _iter_kv(node):
            full_key = f"{dotted}.{key}" if not key.startswith(dotted) else key
            rows.append(
                [
                    _clean_cell(full_key),
                    _clean_cell(value),
                    _clean_cell(_param_desc(full_key)),
                ]
            )
        if not rows:
            continue
        lines.append(f"#### {section_title}")
        lines.append("")
        lines.extend(_md_table(["参数", "取值", "说明"], rows))
        lines.append("")
    # metrics-derived actual settings (more authoritative than config for some stages)
    extra: list[list[str]] = []
    if stage == "gwas":
        for key, label in (
            ("model", "实际模型"),
            ("phenotype_type", "表型类型"),
            ("covariates", "实际协变量"),
        ):
            value = metrics.get(key)
            if value is not None:
                extra.append([label, _clean_cell(value)])
        if metrics.get("lambda_source"):
            extra.append(["λGC 计算来源", _clean_cell(metrics.get("lambda_source"))])
    if stage == "intake":
        species = metrics.get("species_inferred")
        if isinstance(species, Mapping) and species:
            for key, value in _iter_kv(species):
                extra.append([f"物种推断 · {_clean_cell(key)}", _clean_cell(value)])
    if extra:
        lines.append("#### 实际运行设置")
        lines.append("")
        lines.extend(_md_table(["项目", "取值"], extra))
        lines.append("")
    return lines


def _stage_key_metrics(stage: str, metrics: Mapping[str, Any]) -> str:
    """Short one-line metric digest used by the final report stage table."""
    try:
        if stage == "intake":
            return f"样本 {_fmt(metrics.get('n_samples'))}，位点 {_fmt(metrics.get('n_variants'))}，格式 {_fmt(metrics.get('genotype_type'))}"
        if stage == "qc":
            return (
                f"样本 {_fmt(metrics.get('n_samples_out'))}/{_fmt(metrics.get('n_samples_in'))}，"
                f"位点 {_fmt(metrics.get('n_variants_out'))}/{_fmt(metrics.get('n_variants_in'))}"
            )
        if stage == "popgen":
            pca = metrics.get("pca") if isinstance(metrics.get("pca"), Mapping) else {}
            rel = metrics.get("relatedness") if isinstance(metrics.get("relatedness"), Mapping) else {}
            var = pca.get("var_explained")
            pc1 = "—"
            if isinstance(var, list) and var:
                try:
                    pc1 = f"{float(var[0]) * 100:.2f}%"
                except Exception:
                    pc1 = _fmt(var[0])
            return f"PC1 方差 {pc1}，PI_HAT≥0.25 对 {_fmt(rel.get('n_pairs_gt_025'))}"
        if stage == "gwas":
            summary = (
                f"λGC {_fmt_lambda(metrics.get('lambda_gc'))}，显著 {_fmt(metrics.get('n_significant'))}，"
                f"检验位点 {_fmt(metrics.get('n_tested'))}"
            )
            min_p = _min_p_display(metrics.get("p_distribution"))
            if min_p:
                summary += f"，最小 P {min_p}"
            return summary
        if stage == "calibration":
            if _norm_status(metrics.get("status") or metrics.get("overall")) == "skipped":
                return "本次未执行"
            null = metrics.get("null") if isinstance(metrics.get("null"), Mapping) else {}
            val = metrics.get("validation") if isinstance(metrics.get("validation"), Mapping) else {}
            if (
                null.get("lambda_mean") is None
                and not null.get("n_replicates")
                and not (isinstance(null.get("type1"), Mapping) and null.get("type1"))
            ):
                summary = "零模型校准未获得结果"
            else:
                summary = f"零模型 λ 均值 {_fmt_lambda(null.get('lambda_mean'))}"
            n_checks = val.get("n_checks")
            try:
                if n_checks and int(n_checks) > 0:
                    summary += f"，校验通过 {_fmt(val.get('n_pass'))}/{_fmt(n_checks)}"
                else:
                    summary += "，独立校验结果缺失"
            except Exception:
                summary += "，独立校验结果缺失"
            return summary
    except Exception:  # pragma: no cover
        pass
    return "—"


def _stage_impact_lines(stage: str, status: str, metrics: Mapping[str, Any]) -> list[str]:
    badge = _badge(status)
    if status == "skipped":
        return [
            f"本阶段状态为 {badge}（未执行），其输出不会进入最终结果。",
            "最终报告中相应部分标注为“本阶段未执行”，结论仅基于已完成的阶段。",
        ]
    if stage == "intake":
        lines = [
            "数据接入结果决定后续所有阶段可用的样本与位点集合；若样本数或表型缺失比例异常，"
            "应在质控前修正输入文件。",
        ]
        if status == "fail":
            lines.append("❌ 数据接入失败，后续阶段不应继续；请先修复文件/编码问题。")
        elif status == "warn":
            lines.append("⚠️ 存在格式或编码警告，建议在解释结果时纳入局限性说明。")
        return lines
    if stage == "qc":
        lines = [
            f"质控后样本与位点分别保留 {_fmt(metrics.get('n_samples_out'))} 与 {_fmt(metrics.get('n_variants_out'))}；"
            "下游 PCA、亲缘与 GWAS 均基于该数据集。",
        ]
        if status == "fail":
            lines.append("❌ 质控未通过：下游关联结果可能存在系统性偏倚，不建议直接作为育种/论文结论依据。")
        elif status == "warn":
            lines.append("⚠️ 部分阈值触发警告：结果可用于探索性分析，正式报告建议补充敏感性分析。")
        else:
            lines.append("✅ 质控通过：样本与位点质量满足预设阈值。")
        return lines
    if stage == "popgen":
        lines = [
            "群体结构与亲缘结果用于 GWAS 协变量选择；PCA 主成分数量与亲缘样本处理方式会直接影响关联结果。",
        ]
        rel = metrics.get("relatedness")
        n_rel = rel.get("n_pairs_gt_025") if isinstance(rel, Mapping) else None
        if n_rel:
            lines.append(f"⚠️ 检测到 {_fmt(n_rel)} 对近亲样本（PI_HAT ≥ 0.25），建议剔除或使用可处理亲缘结构的模型。")
        if status == "warn":
            lines.append("⚠️ 本阶段存在警告，需在 GWAS 中确认结构校正是否充分。")
        elif status == "fail":
            lines.append("❌ 群体结构分析失败，GWAS 的协变量校正缺乏依据。")
        return lines
    if stage == "gwas":
        lines = []
        lam = metrics.get("lambda_gc")
        try:
            lam_value = float(lam)
        except Exception:
            lam_value = None
        if lam_value is not None and lam_value > 1.1:
            lines.append(f"⚠️ λGC = {lam_value:.3f} 明显大于 1，提示 p 值膨胀，建议增加主成分/剔除近亲或改用混合模型。")
        elif lam_value is not None and lam_value < 0.9:
            lines.append(f"⚠️ λGC = {lam_value:.3f} 偏小，可能提示过度校正或表型转换问题。")
        n_sig = metrics.get("n_significant")
        if n_sig in (0, "0", 0.0):
            lines.append("本次分析未检出 Bonferroni 显著位点；可能受样本量、遗传力或标记密度限制。")
        lead = metrics.get("lead_snps")
        if isinstance(lead, list) and lead:
            top = lead[0] if isinstance(lead[0], Mapping) else {}
            lines.append(
                f"排名第一的 lead SNP 为 { _fmt(top.get('snp')) }（{ _fmt(top.get('chr')) }:{_fmt(top.get('bp'))}，"
                f"P = {_p_display(top)}），可作为候选区域定位起点。"
            )
        if status == "fail":
            lines.insert(0, "❌ GWAS 未通过：结果不可用于下游解释，请检查报错与日志。")
        return lines or ["GWAS 结果用于候选位点筛选，需通过独立群体或功能实验验证。"]
    if stage == "calibration":
        lines = [
            "校准与独立校验为结果可信度提供证据；若零模型 λ 偏离或独立重算失败，"
            "GWAS 显著性应谨慎解读。",
        ]
        if status == "fail":
            lines.append("❌ 校准/校验未通过：流程实现或数据读写可能存在错误，主结果不可信。")
        elif status == "warn":
            lines.append("⚠️ 校准/校验存在警告：建议人工复核对应检查项。")
        return lines
    return []


def _render_results_section(stage: str, stage_dir: Path, metrics: Mapping[str, Any], tables: Mapping[str, Path], figures: Sequence[Path], status: str) -> list[str]:
    lines: list[str] = []
    if status == "skipped":
        lines.append("本阶段未执行，无结果数据。")
        lines.append("")
        return lines
    if not metrics:
        lines.append("本阶段未产生可解析的 metrics 数据（metrics.json 缺失或为空）。")
        lines.append("")
        return lines

    # 3.x scalar metrics
    scalar_rows: list[list[str]] = []
    seen: set[str] = set()
    sub = [0]

    def next_heading(text: str) -> str:
        sub[0] += 1
        return f"#### 3.{sub[0]} {text}"

    _SCALAR_SKIP = {
        "tables",
        "figures",
        "gates",
        "warnings",
        "errors",
        "status",
        "summary",
        "stage",
        "overall",
        "checks",
    }

    def add_scalar(path: str, value: Any) -> None:
        if path in seen:
            return
        seen.add(path)
        scalar_rows.append([_label_for(path, _METRIC_LABELS), _clean_cell(value, max_len=120)])

    def walk(node: Any, prefix: str = "") -> None:
        if not isinstance(node, Mapping):
            return
        for key, value in node.items():
            name = str(key)
            if name in _SCALAR_SKIP or name.lower().endswith(_TEXT_FIELD_SUFFIX):
                continue
            path = f"{prefix}{name}"
            if isinstance(value, Mapping):
                # Use pre-formatted "*_text" siblings (e.g. p_distribution.min_text)
                texts = {
                    str(k)[: -len(_TEXT_FIELD_SUFFIX)].lower(): str(v)
                    for k, v in value.items()
                    if str(k).lower().endswith(_TEXT_FIELD_SUFFIX)
                }
                for child_key, child_value in value.items():
                    child_name = str(child_key)
                    if child_name.lower().endswith(_TEXT_FIELD_SUFFIX):
                        continue
                    child_path = f"{path}.{child_name}"
                    if child_name.lower() in texts:
                        add_scalar(child_path, texts[child_name.lower()])
                        continue
                    if isinstance(child_value, Mapping):
                        walk(child_value, child_path + ".")
                    elif isinstance(child_value, list):
                        if not child_value or all(isinstance(v, Mapping) for v in child_value):
                            continue
                        add_scalar(child_path, _fmt(child_value))
                    else:
                        add_scalar(child_path, child_value)
            elif isinstance(value, list):
                if not value or all(isinstance(v, Mapping) for v in value):
                    continue
                add_scalar(path, _fmt(value))
            else:
                add_scalar(path, value)

    walk(metrics)
    if scalar_rows:
        if len(scalar_rows) > 40:
            note = f"（共 {len(scalar_rows)} 项，展示前 40 项）"
            scalar_rows = scalar_rows[:40]
        else:
            note = ""
        lines.append(next_heading("关键指标"))
        lines.append("")
        if note:
            lines.append(note)
            lines.append("")
        lines.extend(_md_table(["指标", "取值"], scalar_rows))
        lines.append("")

    # 3.x detail tables (lists of dicts)
    detail = _discover_extra_metric_tables(metrics)
    if detail:
        lines.append(next_heading("明细表"))
        lines.append("")
        for key, items in detail[:5]:
            title = _table_label(key)
            columns = _detail_table_columns(items)
            rows: list[list[str]] = []
            for item in items[:10]:
                rows.append([_detail_cell(item, column) for column in columns])
            lines.append(f"**{title}**（{len(items)} 条" + (f"，展示前 {len(rows)} 条" if len(items) > len(rows) else "") + "）")
            lines.append("")
            lines.extend(_md_table([_label_for(c, _METRIC_LABELS) for c in columns], rows))
            lines.append("")

    # 3.x CSV tables
    if tables:
        lines.append(next_heading("数据表"))
        lines.append("")
        for name, path in list(tables.items())[:6]:
            headers, rows, notes = _read_csv_rows(path)
            lines.append(f"**{_table_label(name)}**（`{_rel_path(path, stage_dir)}`）")
            lines.append("")
            if headers:
                lines.extend(_md_table(headers, rows))
                if notes:
                    lines.append("")
                    lines.append("> " + " ".join(notes))
            else:
                lines.append("（无法读取该表）")
            lines.append("")
        if len(tables) > 6:
            lines.append(f"> 另有 {len(tables) - 6} 个数据表未在正文展开，详见阶段目录 `tables/`。")
            lines.append("")

    # 3.x figures
    if figures:
        lines.append(next_heading("图表"))
        lines.append("")
        for index, figure in enumerate(figures, start=1):
            try:
                rel = figure.relative_to(stage_dir).as_posix()
            except Exception:
                rel = figure.name
            caption = _figure_caption(figure)
            lines.append(f"**图 {index}**　{caption}")
            lines.append("")
            lines.append(f"![{caption}]({rel})")
            lines.append("")
    return lines


def stage_report(
    run_dir: Path,
    stage: str,
    metrics: dict,
    config: dict,
    summary_override: str | None = None,
) -> Path:
    """Write ``<run_dir>/<stage dir>/report.md`` and return its path."""
    run_path = Path(run_dir)
    stage_key = str(stage or "").strip().lower()
    stage_dir, meta = _ensure_stage_dir(run_path, stage_key)
    cfg: Mapping[str, Any] = config if isinstance(config, Mapping) else {}
    wrapper: Mapping[str, Any] = metrics if isinstance(metrics, Mapping) else {}
    inner = _extract_metrics(wrapper)

    status = _norm_status(
        wrapper.get("status")
        or inner.get("status")
        or inner.get("overall")
        or ("unknown" if inner else "skipped")
    )
    run_id = _project_title(cfg, run_path)
    title = f"{run_id} — {meta['name']}报告"
    summary = summary_override or wrapper.get("summary") or inner.get("summary") or ""
    warnings = _as_list(wrapper.get("warnings") or inner.get("warnings"))
    errors = _as_list(wrapper.get("errors") or inner.get("errors"))
    tables = _parse_tables(wrapper.get("tables") if wrapper.get("tables") else inner.get("tables"), stage_dir)
    figures = _parse_figures(wrapper.get("figures") if wrapper.get("figures") else inner.get("figures"), stage_dir)

    lines: list[str] = []
    lines.append(f"# {title}")
    lines.append("")
    lines.append(
        f"> **运行编号**：`{cfg.get('run_id', run_path.name)}` ｜ **阶段**：`{stage_key}`（`{meta['dir']}/`）"
        f" ｜ **状态**：{_badge(status)} ｜ **生成时间**：{_now_str()}"
    )
    if cfg.get("created_at"):
        lines.append(f"> **运行创建时间**：{_fmt(cfg.get('created_at'))} ｜ **物种**：{_species_label(cfg, inner)}")
    lines.append("")

    # ---- 阶段概要 -------------------------------------------------------
    lines.append("## 阶段概要")
    lines.append("")
    if summary:
        lines.append(str(summary).strip())
    else:
        auto = f"本阶段状态为 {_badge(status)}。"
        if stage_key == "qc":
            auto += (
                f"质控前 {_fmt(inner.get('n_samples_in'))} 个样本、{_fmt(inner.get('n_variants_in'))} 个位点，"
                f"质控后 {_fmt(inner.get('n_samples_out'))} 个样本、{_fmt(inner.get('n_variants_out'))} 个位点。"
            )
        lines.append(auto)
    lines.append("")
    if status == "skipped":
        lines.append("> ⏭️ **本阶段未执行**：以下章节仅保留报告结构，无实际结果；最终结论不依赖本阶段输出。")
        lines.append("")

    # ---- 一、阶段目的 ---------------------------------------------------
    lines.append("## 一、阶段目的")
    lines.append("")
    lines.append(meta["purpose"])
    lines.append("")

    # ---- 二、方法与参数 -------------------------------------------------
    lines.append("## 二、方法与参数")
    lines.append("")
    if meta["methods"]:
        lines.append("### 2.1 分析方法")
        lines.append("")
        for item in meta["methods"]:
            lines.append(f"- {item}")
        lines.append("")
    else:
        lines.append("本阶段为通用阶段，未定义专门的分析方法说明。")
        lines.append("")
    lines.append("### 2.2 参数设置")
    lines.append("")
    if status == "skipped":
        lines.append("本阶段未执行，未读取或应用参数。")
        lines.append("")
    else:
        param_lines = _stage_param_lines(stage_key, cfg, inner)
        if param_lines:
            lines.extend(param_lines)
        else:
            lines.append("未找到本阶段对应的配置参数（run.json 中相应段落缺失），使用流程默认值。")
            lines.append("")
        runtime = cfg.get("runtime") if isinstance(cfg.get("runtime"), Mapping) else {}
        if runtime:
            lines.append(
                f"> 运行时：线程 {_fmt(runtime.get('threads'))}，内存上限 {_fmt(runtime.get('memory_mb'))} MB，"
                f"单步超时 {_fmt(runtime.get('timeout_s'))} 秒。"
            )
            lines.append("")

    # ---- 三、结果 -------------------------------------------------------
    lines.append("## 三、结果")
    lines.append("")
    lines.extend(_render_results_section(stage_key, stage_dir, inner, tables, figures, status))

    # ---- 四、质控判定 ---------------------------------------------------
    lines.append("## 四、质控判定")
    lines.append("")
    if status == "skipped":
        lines.append("本阶段未执行，无判定结果。")
    else:
        lines.extend(_render_gates(inner))
    lines.append("")

    # ---- 五、异常与处理 -------------------------------------------------
    lines.append("## 五、异常与处理")
    lines.append("")
    if warnings:
        lines.append("**警告（WARN）**")
        lines.append("")
        for item in warnings:
            lines.append(f"- ⚠️ {item}")
        lines.append("")
    if errors:
        lines.append("**错误（ERROR）**")
        lines.append("")
        for item in errors:
            lines.append(f"- ❌ {item}")
        lines.append("")
    if not warnings and not errors:
        lines.append("本阶段未记录警告或错误。")
        lines.append("")
    else:
        lines.append("**处理建议**")
        lines.append("")
        for hint in meta.get("handling", []):
            lines.append(f"- {hint}")
        lines.append("")

    # ---- 六、对后续分析的影响 -------------------------------------------
    lines.append("## 六、对后续分析的影响")
    lines.append("")
    for item in _stage_impact_lines(stage_key, status, inner):
        lines.append(f"- {item}")
    lines.append("")

    # ---- 七、产物清单 ---------------------------------------------------
    lines.append("## 七、产物清单")
    lines.append("")
    artifact_rows: list[list[str]] = []
    for path in sorted(stage_dir.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(stage_dir).as_posix()
        if rel == "report.md":
            continue
        artifact_rows.append([rel, _human_size(path.stat().st_size), _artifact_description(rel)])
    if artifact_rows:
        lines.extend(_md_table(["文件", "大小", "说明"], artifact_rows[:60]))
        if len(artifact_rows) > 60:
            lines.append(f"> 另有 {len(artifact_rows) - 60} 个文件未列出。")
    else:
        lines.append("本阶段目录暂无文件产物。")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("> 本报告由 plink-analysis skill 自动生成；最终汇总见 `../05_reports/final_report.md`（如已生成）。")

    report_path = stage_dir / "report.md"
    _write_text(report_path, _renumber_figures("\n".join(lines) + "\n"))
    return report_path


def _artifact_description(rel: str) -> str:
    name = Path(rel).name
    if name == "metrics.json":
        return "阶段指标（机器可读）"
    if name == "validation.json":
        return "独立校验结果"
    if name == "report.md":
        return "阶段报告（Markdown）"
    if "/tables/" in rel or rel.startswith("tables/"):
        return "数据表（CSV）"
    if "/plots/" in rel or rel.startswith("plots/"):
        return "图（PNG）"
    if "/data/" in rel or rel.startswith("data/"):
        return "中间数据"
    if name.endswith(".log"):
        return "日志"
    return "文件"


def _merge_validation(item: Mapping[str, Any]) -> dict:
    """Merge ``validation.json`` (authoritative checks) with metrics.validation counts."""
    metrics = item.get("metrics") if isinstance(item.get("metrics"), Mapping) else {}
    metric_validation = metrics.get("validation") if isinstance(metrics.get("validation"), Mapping) else {}
    file_validation = item.get("validation") if isinstance(item.get("validation"), Mapping) else {}
    merged: dict[str, Any] = {}
    merged.update(metric_validation)
    merged.update(file_validation)
    checks = merged.get("checks")
    if isinstance(checks, list) and checks:
        statuses = [_norm_status(check.get("status")) for check in checks if isinstance(check, Mapping)]
        merged.setdefault("n_checks", len(statuses))
        merged.setdefault("n_pass", statuses.count("pass"))
        merged.setdefault("n_warn", statuses.count("warn"))
        merged.setdefault("n_fail", statuses.count("fail"))
    if merged.get("overall") is None and merged.get("n_checks") is not None:
        if merged.get("n_fail"):
            merged["overall"] = "fail"
        elif merged.get("n_warn"):
            merged["overall"] = "warn"
        else:
            merged["overall"] = "pass"
    return merged


# ---------------------------------------------------------------------------
# Final report helpers
# ---------------------------------------------------------------------------
def _collect_stages(run_dir: Path, state: Mapping[str, Any] | None, config: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    stages: dict[str, dict[str, Any]] = {}
    state_stages: Mapping[str, Any] = {}
    if isinstance(state, Mapping) and isinstance(state.get("stages"), Mapping):
        state_stages = state["stages"]
    for stage in STAGE_ORDER:
        stage_dir = run_dir / STAGE_DIRS[stage]
        raw = _read_json(stage_dir / "metrics.json", None)
        inner = _extract_metrics(raw) if raw is not None else {}
        state_item = state_stages.get(stage) if isinstance(state_stages, Mapping) else None
        state_item = state_item if isinstance(state_item, Mapping) else {}
        report_path = stage_dir / "report.md"
        if state_item.get("status"):
            status = _norm_status(state_item.get("status"))
        elif raw is not None:
            status = _norm_status(inner.get("status") or inner.get("overall") or "unknown")
        else:
            status = "skipped"
        stages[stage] = {
            "stage": stage,
            "dir": stage_dir,
            "name": _STAGE_META[stage]["name"],
            "raw": raw if isinstance(raw, Mapping) else {},
            "metrics": inner,
            "state": dict(state_item),
            "status": status,
            "report": report_path,
            "report_exists": report_path.exists(),
            "validation": _read_json(stage_dir / "validation.json", None),
        }
    return stages


def _summary_of_stage(item: Mapping[str, Any]) -> str:
    metrics = item.get("metrics") if isinstance(item.get("metrics"), Mapping) else {}
    summary = metrics.get("summary")
    if summary:
        return str(summary).strip()
    report_path = item.get("report")
    if isinstance(report_path, Path) and report_path.exists():
        try:
            text = report_path.read_text(encoding="utf-8", errors="replace")
            match = re.search(r"##\s*阶段概要\s*\n+(.+?)(?:\n\s*\n|\n##)", text, re.S)
            if match:
                return match.group(1).strip()
        except Exception:
            pass
    return ""


def _stage_summary_block(item: Mapping[str, Any], *, max_chars: int = 420) -> list[str]:
    """Quote the stage report conclusion (read from report.md when available)."""
    summary = " ".join(_summary_of_stage(item).split())
    if not summary:
        return []
    if len(summary) > max_chars:
        summary = summary[: max_chars - 1] + "…"
    return [f"> **阶段结论**：{summary}", ""]


def _collect_figures_for_final(stages: Mapping[str, Mapping[str, Any]], figures_dir: Path) -> dict[str, list[tuple[str, str]]]:
    """Copy every stage figure into 05_reports/figures and return caption/path pairs."""
    result: dict[str, list[tuple[str, str]]] = {stage: [] for stage in STAGE_ORDER}
    figures_dir.mkdir(parents=True, exist_ok=True)
    for stage in STAGE_ORDER:
        stage_dir = stages[stage]["dir"]
        plot_dir = stage_dir / "plots" if isinstance(stage_dir, Path) else None
        if not plot_dir or not plot_dir.is_dir():
            continue
        figures = sorted(plot_dir.glob("*.png"))
        for index, source in enumerate(figures[:10]):
            filename = source.name if source.name.startswith(f"{stage}_") else f"{stage}_{source.name}"
            destination = figures_dir / filename
            try:
                if destination.exists():
                    destination.unlink()
                destination.write_bytes(source.read_bytes())
            except Exception:
                continue
            result[stage].append((_figure_caption(source), f"figures/{destination.name}"))
    return result


def _render_figure_block(captions: Sequence[tuple[str, str]], start: int = 1) -> list[str]:
    lines: list[str] = []
    for offset, (caption, rel) in enumerate(captions):
        lines.append(f"**图 {start + offset}**　{caption}")
        lines.append("")
        lines.append(f"![{caption}]({rel})")
        lines.append("")
    return lines


def _percent(numerator: Any, denominator: Any) -> str:
    try:
        num = float(numerator)
        den = float(denominator)
        if den == 0:
            return "—"
        return f"{num / den * 100:.1f}%"
    except Exception:
        return "—"


def _stage_status_table(stages: Mapping[str, Mapping[str, Any]], run_dir: Path) -> list[str]:
    rows: list[list[str]] = []
    for stage in STAGE_ORDER:
        item = stages[stage]
        report_rel = _rel_path(item["report"], run_dir) if item["report_exists"] else "—"
        elapsed = item.get("state", {}).get("elapsed_s")
        rows.append(
            [
                _STAGE_META[stage]["name"],
                f"`{STAGE_DIRS[stage]}/`",
                _badge(item["status"]),
                _clean_cell(elapsed),
                _clean_cell(
                    "本次未执行"
                    if _norm_status(item["status"]) == "skipped"
                    else _stage_key_metrics(stage, item["metrics"])
                ),
                f"[{report_rel}]({ '..' if run_dir.name else '.' }/{report_rel})" if item["report_exists"] else "—",
            ]
        )
    lines = ["| 阶段 | 目录 | 状态 | 耗时(s) | 关键结果 | 阶段报告 |", "| :--- | :--- | :---: | ---: | :--- | :--- |"]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return lines


def _section_executive_summary(
    stages: Mapping[str, Mapping[str, Any]], config: Mapping[str, Any], overall: str
) -> list[str]:
    lines = ["## 一、执行摘要", ""]
    intake = stages["intake"]["metrics"]
    qc = stages["qc"]["metrics"]
    popgen = stages["popgen"]["metrics"]
    gwas = stages["gwas"]["metrics"]
    calibration = stages["calibration"]["metrics"]

    bullets: list[str] = []
    if intake:
        bullets.append(
            f"数据概况：{_fmt(intake.get('n_samples'))} 个样本 × {_fmt(intake.get('n_variants'))} 个位点，"
            f"格式 {_fmt(intake.get('genotype_type'))}，物种 {_species_label(config, intake)}。"
        )
    if qc:
        bullets.append(
            f"质量控制：样本保留率 {_percent(qc.get('n_samples_out'), qc.get('n_samples_in'))}"
            f"（{_fmt(qc.get('n_samples_out'))}/{_fmt(qc.get('n_samples_in'))}），"
            f"位点保留率 {_percent(qc.get('n_variants_out'), qc.get('n_variants_in'))}"
            f"（{_fmt(qc.get('n_variants_out'))}/{_fmt(qc.get('n_variants_in'))}），"
            f"总体判定 {_badge(qc.get('overall') or qc.get('status'))}。"
        )
    if popgen:
        pca = popgen.get("pca") if isinstance(popgen.get("pca"), Mapping) else {}
        rel = popgen.get("relatedness") if isinstance(popgen.get("relatedness"), Mapping) else {}
        var = pca.get("var_explained")
        pc1 = "—"
        if isinstance(var, list) and var:
            try:
                pc1 = f"{float(var[0]) * 100:.2f}%"
            except Exception:
                pc1 = _fmt(var[0])
        bullets.append(
            f"群体结构：PC1 解释 {pc1} 方差，PCA 离群样本 {_fmt(pca.get('n_outliers'))} 个，"
            f"PI_HAT ≥ 0.25 的样本对 {_fmt(rel.get('n_pairs_gt_025'))} 对。"
        )
    if gwas:
        top = ""
        lead = gwas.get("lead_snps")
        if isinstance(lead, list) and lead and isinstance(lead[0], Mapping):
            first = lead[0]
            top = f"，最显著位点 {_fmt(first.get('snp'))}（P = {_p_display(first)}）"
        else:
            min_p = _min_p_display(gwas.get("p_distribution"))
            if min_p:
                top = f"，最小 P = {min_p}"
        bullets.append(
            f"GWAS：{_fmt(gwas.get('model'))} 模型，检验 {_fmt(gwas.get('n_tested'))} 个位点，"
            f"λGC = {_fmt_lambda(gwas.get('lambda_gc'))}，Bonferroni 显著 {_fmt(gwas.get('n_significant'))} 个{top}。"
        )
    cal_item = stages["calibration"]
    cal_status = _norm_status(cal_item.get("status"))
    has_null, has_power, n_checks, validation = _calibration_evidence(cal_item)
    if cal_status == "skipped":
        reason = _short_reason(_skip_reason(cal_item, "仿真校准未执行（可能未启用或依赖未就绪）")).rstrip("。；; ")
        reason = reason or "未提供具体原因"
        bullets.append(
            f"校准与独立校验：**本次未执行**（原因：{reason}）；p 值校准与统计效力未经独立验证。"
        )
    elif not has_null and n_checks == 0 and not has_power:
        reason = _short_reason(_skip_reason(cal_item, "未获得有效的仿真校准结果")).rstrip("。；; ")
        reason = reason or "未提供具体原因"
        bullets.append(
            f"校准与独立校验：未获得有效校准结果（原因：{reason}）；p 值校准与统计效力未经独立验证。"
        )
    elif calibration or n_checks:
        parts: list[str] = []
        if has_null:
            null = calibration.get("null") if isinstance(calibration.get("null"), Mapping) else {}
            parts.append(f"零模型 λGC 均值 {_fmt_lambda(null.get('lambda_mean'))}")
        else:
            parts.append("零模型校准未获得结果")
        if n_checks > 0:
            parts.append(f"独立校验通过 {_fmt(validation.get('n_pass'))}/{_fmt(n_checks)} 项")
        else:
            parts.append("独立校验结果缺失")
        parts.append(f"总体判定 {_badge(calibration.get('overall') or cal_status)}")
        bullets.append("校准与校验：" + "，".join(parts) + "。")
    if not bullets:
        bullets.append("各阶段尚未产生可用的 metrics.json，本报告仅包含结构与复现信息。")

    lines.extend(f"- {item}" for item in bullets)
    lines.append("")
    lines.append(f"**综合判定：{_badge(overall)}**")
    lines.append("")
    verdict = {
        "pass": "所有已执行阶段均通过预设阈值，结果可用于候选位点筛选与育种决策参考；仍需注意观察性关联不等于因果。",
        "warn": "流程整体可用，但存在警告项；建议按各阶段报告中的处理建议复核后，再用于正式结论。",
        "fail": "存在失败阶段，结果不可直接用于结论；请先依据阶段报告与 logs 修复问题并重跑。",
        "skipped": "未执行任何分析阶段，本报告仅包含复现命令与结构占位。",
        "unknown": "阶段状态不完整，请检查 state.json 与各阶段 metrics.json。",
    }.get(_norm_status(overall), "请检查各阶段报告。")
    lines.append(verdict)
    lines.append("")
    recommendations = _recommendations(stages)
    if recommendations:
        lines.append("**建议关注**")
        lines.append("")
        lines.extend(f"- {item}" for item in recommendations)
        lines.append("")
    return lines


def _recommendations(stages: Mapping[str, Mapping[str, Any]]) -> list[str]:
    tips: list[str] = []
    qc = stages["qc"]["metrics"]
    gwas = stages["gwas"]["metrics"]
    cal = stages["calibration"]["metrics"]
    try:
        lam = float(gwas.get("lambda_gc"))
        if lam > 1.1:
            tips.append(
                "λGC 明显大于 1，建议增加 GWAS 的 PCA 主成分数、剔除近亲样本，或使用可处理亲缘结构的混合模型复核。"
            )
    except Exception:
        pass
    try:
        n_samples = int(qc.get("n_samples_out") or 0)
        if 0 < n_samples < 200:
            tips.append(f"质控后样本量仅 {n_samples}，检验功效有限；显著位点需谨慎解释，建议扩大样本或独立群体验证。")
    except Exception:
        pass
    cal_item = stages["calibration"]
    cal_status = _norm_status(cal_item.get("status"))
    has_null, _has_power, _n_checks, _validation = _calibration_evidence(cal_item)
    if cal_status in ("skipped", "unknown") or not has_null:
        tips.append("仿真校准未获得有效结果，无法对 p 值校准与统计效力提供独立证据。")
    if _norm_status(stages["popgen"]["status"]) == "skipped":
        tips.append("群体结构与亲缘分析未执行，GWAS 协变量选择缺少数据支持。")
    if _norm_status(qc.get("overall")) == "fail":
        tips.append("质控未通过，后续结果存在系统性偏倚风险。")
    if not tips:
        tips.append("结果解读时应结合表型测量误差、芯片覆盖度与群体历史，避免将统计关联直接解释为因果。")
    return tips


def _section_data_overview(stages: Mapping[str, Mapping[str, Any]], config: Mapping[str, Any], figures: Mapping[str, list[tuple[str, str]]]) -> list[str]:
    lines = ["## 三、数据接入与数据概况", ""]
    intake = stages["intake"]
    metrics = intake["metrics"]
    if intake["status"] == "skipped" or not metrics:
        lines.append("> ⏭️ 数据接入阶段未执行或未产生结果，无法呈现数据概况。")
        lines.append("")
        return lines
    lines.extend(_stage_summary_block(intake))
    rows = [
        ["基因型格式", _fmt(metrics.get("genotype_type"))],
        ["样本数", _fmt(metrics.get("n_samples"))],
        ["位点数", _fmt(metrics.get("n_variants"))],
        ["染色体", _fmt(metrics.get("chromosomes"))],
        ["推断物种", _species_label(config, metrics)],
        ["表型", _format_phenotype(metrics.get("phenotype"))],
        ["协变量", _format_covariates(metrics.get("covariates"))],
        ["分组", _format_groups(metrics.get("groups"))],
        ["预检状态", _badge((metrics.get("preflight") or {}).get("status") if isinstance(metrics.get("preflight"), Mapping) else "unknown")],
    ]
    lines.extend(_md_table(["项目", "取值"], rows))
    lines.append("")
    species = metrics.get("species_inferred")
    if isinstance(species, Mapping) and species:
        lines.append("**物种推断明细**")
        lines.append("")
        lines.extend(_md_table(["项目", "取值"], [[_label_for(f"species_inferred.{k}", _METRIC_LABELS), _clean_cell(v)] for k, v in _iter_kv(species)]))
        lines.append("")
    gates = _gate_rows(metrics)
    if gates:
        lines.append("**输入文件预检**")
        lines.append("")
        gate_rows: list[list[str]] = []
        for item in gates:
            gate_rows.append(
                [
                    _clean_cell(item.get("name") or item.get("check") or "—"),
                    _badge(item.get("status")),
                    _clean_cell(item.get("detail", "—"), max_len=120),
                ]
            )
        lines.extend(_md_table(["检查项", "状态", "说明"], gate_rows))
        lines.append("")
    files = metrics.get("files")
    if isinstance(files, list) and files:
        lines.append("**输入文件清单**")
        lines.append("")
        file_rows: list[list[str]] = []
        for item in files[:12]:
            if isinstance(item, Mapping):
                file_rows.append(
                    [
                        _clean_cell(item.get("path", "—"), max_len=100),
                        _human_size(item.get("size", 0)),
                        _clean_cell(item.get("sha256", "—"), max_len=24),
                    ]
                )
        if file_rows:
            lines.extend(_md_table(["文件", "大小", "SHA-256（前 24 位）"], file_rows))
            lines.append("")
    lines.extend(_render_figure_block(figures.get("intake", [])))
    return lines


def _section_qc(stages: Mapping[str, Mapping[str, Any]], figures: Mapping[str, list[tuple[str, str]]]) -> list[str]:
    lines = ["## 四、质量控制（QC）", ""]
    item = stages["qc"]
    metrics = item["metrics"]
    if item["status"] == "skipped" or not metrics:
        lines.append("> ⏭️ 质控阶段未执行或未产生结果，后续关联分析的样本/位点质量无法确认。")
        lines.append("")
        return lines
    lines.extend(_stage_summary_block(item))
    lines.append("### 4.1 质控前后对比")
    lines.append("")
    rows = [
        ["样本数", _fmt(metrics.get("n_samples_in")), _fmt(metrics.get("n_samples_out")), _percent(metrics.get("n_samples_out"), metrics.get("n_samples_in"))],
        ["位点数", _fmt(metrics.get("n_variants_in")), _fmt(metrics.get("n_variants_out")), _percent(metrics.get("n_variants_out"), metrics.get("n_variants_in"))],
    ]
    lines.extend(_md_table(["项目", "质控前", "质控后", "保留率"], rows))
    lines.append("")
    removed = metrics.get("removed")
    if isinstance(removed, Mapping):
        lines.append("### 4.2 剔除明细")
        lines.append("")
        detail_rows: list[list[str]] = []
        samples = removed.get("samples")
        variants = removed.get("variants")
        if isinstance(samples, Mapping):
            for key, value in samples.items():
                detail_rows.append(["样本", _label_for(str(key), {"missingness": "缺失率", "het_outlier": "杂合度离群", "sex_discordant": "性别不一致", "duplicate": "重复样本", "related": "近亲样本", "pca_outlier": "PCA 离群"}), _fmt(value)])
        if isinstance(variants, Mapping):
            for key, value in variants.items():
                detail_rows.append(["位点", _label_for(str(key), {"missingness": "缺失率", "maf": "MAF 过低", "hwe": "HWE 偏分离", "monomorphic": "单态位点"}), _fmt(value)])
        if detail_rows:
            lines.extend(_md_table(["层面", "剔除原因", "数量"], detail_rows))
            lines.append("")
    lines.append("### 4.3 质量指标")
    lines.append("")
    metric_rows: list[list[str]] = []

    def add(label: str, value: Any) -> None:
        if value is not None:
            metric_rows.append([label, _clean_cell(value)])

    missing = metrics.get("missing_rate") if isinstance(metrics.get("missing_rate"), Mapping) else {}
    add("样本平均缺失率", missing.get("sample_mean"))
    add("位点平均缺失率", missing.get("variant_mean"))
    maf = metrics.get("maf") if isinstance(metrics.get("maf"), Mapping) else {}
    add("平均 MAF", maf.get("mean"))
    add("MAF < 0.01 位点数", maf.get("n_below_001"))
    hwe = metrics.get("hwe") if isinstance(metrics.get("hwe"), Mapping) else {}
    add("HWE p < 1e-6 位点数", hwe.get("n_p_lt_1e-6"))
    add("HWE λ", hwe.get("lambda"))
    het = metrics.get("het") if isinstance(metrics.get("het"), Mapping) else {}
    add("杂合度均值 F", het.get("mean"))
    add("杂合度标准差", het.get("sd"))
    sex = metrics.get("sex") if isinstance(metrics.get("sex"), Mapping) else {}
    add("性别检查", "已执行" if sex.get("checked") else ("未执行" if sex else None))
    add("性别不一致样本数", len(sex.get("discordant")) if isinstance(sex.get("discordant"), list) else sex.get("discordant"))
    rel = metrics.get("relatedness") if isinstance(metrics.get("relatedness"), Mapping) else {}
    add("重复样本对数", rel.get("n_duplicates"))
    add("近亲样本对数", rel.get("n_related"))
    add("最大 PI_HAT", rel.get("max_pihat"))
    add("平均 PI_HAT", rel.get("mean_pihat"))
    pca = metrics.get("pca") if isinstance(metrics.get("pca"), Mapping) else {}
    add("PCA 离群样本数", pca.get("n_outliers"))
    pruned = metrics.get("pruned") if isinstance(metrics.get("pruned"), Mapping) else {}
    add("LD 剪枝后位点数", pruned.get("n_variants"))
    if metric_rows:
        lines.extend(_md_table(["指标", "取值"], metric_rows))
        lines.append("")
    gates = _gate_rows(metrics)
    if gates:
        lines.append("### 4.4 质控判定")
        lines.append("")
        lines.extend(_render_gates(metrics))
        lines.append("")
    lines.extend(_render_figure_block(figures.get("qc", [])))
    return lines


def _section_popgen(stages: Mapping[str, Mapping[str, Any]], figures: Mapping[str, list[tuple[str, str]]]) -> list[str]:
    lines = ["## 五、群体结构与亲缘分析", ""]
    item = stages["popgen"]
    metrics = item["metrics"]
    if item["status"] == "skipped" or not metrics:
        lines.append("> ⏭️ 群体结构与亲缘分析未执行或未产生结果，GWAS 的群体分层校正缺少依据。")
        lines.append("")
        return lines
    lines.extend(_stage_summary_block(item))
    pca = metrics.get("pca") if isinstance(metrics.get("pca"), Mapping) else {}
    if pca:
        lines.append("### 5.1 主成分分析（PCA）")
        lines.append("")
        var = pca.get("var_explained")
        pca_rows: list[list[str]] = [["主成分数", _fmt(pca.get("n_pcs"))], ["离群样本数", _fmt(pca.get("n_outliers"))]]
        if isinstance(var, list) and var:
            # NOTE: stage_popgen/stage_qc report var_explained as a FRACTION OF
            # ALL N SAMPLES (PLINK 1.9 eigenvalues sum to ~N because --pca works
            # on the variance-standardised GRM).  Do not renormalise by the sum
            # of the first k eigenvalues here -- that would inflate PC1 ~N/k.
            for index, value in enumerate(var[:10], start=1):
                try:
                    pct = f"{float(value) * 100:.3f}%"
                except Exception:
                    pct = _fmt(value)
                pca_rows.append([f"PC{index} 方差解释", pct])
            try:
                cumulative = sum(float(v) for v in var[:5]) * 100
                pca_rows.append(["PC1–PC5 累计方差（占总方差）", f"{cumulative:.2f}%"])
            except Exception:
                pass
        lines.extend(_md_table(["项目", "取值"], pca_rows))
        lines.append("")
        outliers = pca.get("outliers")
        if isinstance(outliers, list) and outliers:
            lines.append(f"**离群样本**：{_fmt(outliers)}")
            lines.append("")
    rel = metrics.get("relatedness") if isinstance(metrics.get("relatedness"), Mapping) else {}
    if rel:
        lines.append("### 5.2 亲缘关系")
        lines.append("")
        lines.extend(
            _md_table(
                ["项目", "取值"],
                [
                    ["PI_HAT ≥ 0.25 的样本对数", _fmt(rel.get("n_pairs_gt_025"))],
                    ["最大 PI_HAT", _fmt(rel.get("max"))],
                ],
            )
        )
        lines.append("")
    ld = metrics.get("ld") if isinstance(metrics.get("ld"), Mapping) else {}
    if ld:
        lines.append("### 5.3 连锁不平衡（LD）")
        lines.append("")
        lines.append(f"共统计 {_fmt(ld.get('n_pairs'))} 个位点对；按物理距离分箱的平均 r² 如下表。")
        lines.append("")
        bins = ld.get("bins")
        if isinstance(bins, list) and bins and all(isinstance(b, Mapping) for b in bins):
            columns = ["distance_bp", "mean_r2", "n_pairs"]
            rows = [[_clean_cell(b.get(c)) for c in columns] for b in bins[:10]]
            lines.extend(_md_table(["距离（bp）", "平均 r²", "位点对数"], rows))
            lines.append("")
    roh = metrics.get("roh") if isinstance(metrics.get("roh"), Mapping) else {}
    if roh:
        lines.append("### 5.4 纯合片段（ROH）与近交")
        lines.append("")
        lines.extend(
            _md_table(
                ["指标", "取值"],
                [
                    ["ROH 数量", _fmt(roh.get("n_roh"))],
                    ["平均每样本 ROH 数", _fmt(roh.get("mean_roh_per_indiv"))],
                    ["平均 FROH", _fmt(roh.get("mean_froh"))],
                ],
            )
        )
        lines.append("")
        top = roh.get("top_samples")
        if isinstance(top, list) and top and isinstance(top[0], Mapping):
            columns = ["iid", "froh", "n_roh", "total_kb"]
            lines.extend(_md_table([_label_for(c, _METRIC_LABELS) for c in columns], [[_clean_cell(t.get(c)) for c in columns] for t in top[:10]]))
            lines.append("")
    fst = metrics.get("fst") if isinstance(metrics.get("fst"), Mapping) else {}
    if fst:
        lines.append("### 5.5 群体分化（FST）")
        lines.append("")
        lines.extend(
            _md_table(
                ["指标", "取值"],
                [["全局 FST", _fmt(fst.get("global"))], ["参与位点数", _fmt(fst.get("n_variants"))]],
            )
        )
        lines.append("")
        top = fst.get("top")
        if isinstance(top, list) and top and isinstance(top[0], Mapping):
            columns = ["snp", "chr", "bp", "fst"]
            lines.extend(_md_table([_label_for(c, _METRIC_LABELS) for c in columns], [[_clean_cell(t.get(c)) for c in columns] for t in top[:10]]))
            lines.append("")
    diversity = metrics.get("diversity") if isinstance(metrics.get("diversity"), Mapping) else {}
    if diversity:
        lines.append("### 5.6 遗传多样性")
        lines.append("")
        lines.extend(
            _md_table(
                ["指标", "取值"],
                [["平均杂合度", _fmt(diversity.get("mean_het"))], ["分组数", _fmt(diversity.get("n_groups"))]],
            )
        )
        lines.append("")
    lines.extend(_render_figure_block(figures.get("popgen", [])))
    return lines


def _section_gwas(stages: Mapping[str, Mapping[str, Any]], figures: Mapping[str, list[tuple[str, str]]]) -> list[str]:
    lines = ["## 六、全基因组关联分析（GWAS）", ""]
    item = stages["gwas"]
    metrics = item["metrics"]
    if item["status"] == "skipped" or not metrics:
        lines.append("> ⏭️ GWAS 阶段未执行或未产生结果，本报告无法给出关联结论。")
        lines.append("")
        return lines
    lines.extend(_stage_summary_block(item))
    sub = [0]

    def next_heading(text: str) -> str:
        sub[0] += 1
        return f"### 6.{sub[0]} {text}"

    lines.append(next_heading("模型与参数"))
    lines.append("")
    covariates = metrics.get("covariates")
    rows = [
        ["关联模型", _fmt(metrics.get("model"))],
        ["表型类型", _fmt(metrics.get("phenotype_type"))],
        ["参与样本数", _fmt(metrics.get("n_samples"))],
        ["检验位点数", _fmt(metrics.get("n_tested"))],
        ["协变量", _fmt(covariates)],
        ["λGC", _fmt_lambda(metrics.get("lambda_gc"))],
        ["λGC 来源", _fmt(metrics.get("lambda_source"))],
        ["总体判定", _badge(metrics.get("overall") or metrics.get("status"))],
    ]
    lines.extend(_md_table(["项目", "取值"], rows))
    lines.append("")
    thresholds = metrics.get("thresholds")
    if isinstance(thresholds, Mapping) and thresholds:
        lines.append(next_heading("显著性阈值与结果计数"))
        lines.append("")
        threshold_rows = [[_clean_cell(k), _clean_cell(v)] for k, v in thresholds.items()]
        lines.extend(_md_table(["阈值类型", "数值"], threshold_rows))
        lines.append("")
    lines.extend(
        _md_table(
            ["统计口径", "数量"],
            [
                ["Bonferroni 显著", _fmt(metrics.get("n_significant"))],
                ["FDR < 0.05", _fmt(metrics.get("n_fdr05"))],
                ["建议性阈值", _fmt(metrics.get("n_suggestive"))],
            ],
        )
    )
    lines.append("")
    p_distribution = metrics.get("p_distribution") if isinstance(metrics.get("p_distribution"), Mapping) else {}
    if p_distribution:
        p_rows: list[list[str]] = []
        min_p_text = _min_p_display(p_distribution)
        if min_p_text:
            p_rows.append(["最小 P 值", min_p_text])
        if p_distribution.get("median") is not None:
            p_rows.append(["P 值中位数", _p_value_text(p_distribution.get("median"))])
        if p_distribution.get("n_p_lt_5e8") is not None:
            p_rows.append(["P < 5e-8 的位点数", _fmt(p_distribution.get("n_p_lt_5e8"))])
        if p_rows:
            lines.append("**P 值分布**")
            lines.append("")
            lines.extend(_md_table(["项目", "取值"], p_rows))
            lines.append("")
    lead = metrics.get("lead_snps")
    if isinstance(lead, list) and lead and all(isinstance(v, Mapping) for v in lead):
        lines.append(next_heading("Lead SNP"))
        lines.append("")
        available = ["snp", "chr", "bp", "a1", "a2", "af", "beta", "se", "p", "n"]
        columns = [c for c in available if any(c in entry for entry in lead)]
        if any("p_text" in entry for entry in lead) and "p" not in columns:
            columns.append("p")
        rows: list[list[str]] = []
        for entry in lead[:20]:
            rows.append([_p_display(entry) if column == "p" else _clean_cell(entry.get(column)) for column in columns])
        lines.extend(_md_table([_label_for(c, _METRIC_LABELS) for c in columns], rows))
        if len(lead) > 20:
            lines.append("")
            lines.append(f"> 共 {len(lead)} 个 lead SNP，此处展示前 20 个；完整结果见 `03_gwas/tables/`。")
        lines.append("")
    lambda_value = metrics.get("lambda_gc")
    try:
        lam = float(lambda_value)
    except Exception:
        lam = None
    if lam is not None:
        lines.append(next_heading("结果解读"))
        lines.append("")
        if lam > 1.1:
            lines.append(
                f"λGC = {lam:.3f} 提示 p 值整体膨胀，可能存在残余群体分层、亲缘结构或表型转换/离群值问题；"
                "在将显著位点视为候选之前，建议增加主成分、剔除近亲样本或使用混合模型复核。"
            )
        elif lam < 0.9:
            lines.append(
                f"λGC = {lam:.3f} 偏小，可能由协变量过度校正或表型分布偏态造成；需结合 QQ 图与校准阶段结果判断。"
            )
        else:
            lines.append(f"λGC = {lam:.3f}，p 值分布与零假设基本一致，未见明显系统性膨胀。")
        lines.append("")
    lines.extend(_render_figure_block(figures.get("gwas", [])))
    return lines


def _section_calibration(stages: Mapping[str, Mapping[str, Any]], figures: Mapping[str, list[tuple[str, str]]]) -> list[str]:
    lines = ["## 七、仿真校准与独立校验", ""]
    item = stages["calibration"]
    metrics = item["metrics"]
    status = _norm_status(item["status"])
    has_null, has_power, n_checks, validation = _calibration_evidence(item)
    default_reason = (
        "仿真校准模块未启用或依赖未就绪（例如 simulate.py 不可用）"
        if status == "skipped"
        else "零模型/已知 QTL 仿真未产生有效结果"
    )
    reason = _short_reason(_skip_reason(item, default_reason)).rstrip("。；; ")
    reason = reason or "未提供具体原因"
    if status == "skipped":
        lines.append(
            f"> ⏭️ **本次未执行**（原因：{reason}）。未提供零模型校准（λGC、假阳性率）、"
            "统计功效与独立重算证据；p 值校准与统计效力未经独立验证，主结果的可信度仅依赖方法与质控。"
        )
        lines.append("")
        return lines
    if not has_null and not has_power and n_checks == 0:
        lines.append(
            f"> ⚠️ **本次未获得有效的仿真校准结果**（原因：{reason}）。"
            "为避免误读，本节不展示 0/0 类占位统计；p 值校准与统计效力未经独立验证。"
        )
        lines.append("")
        return lines
    lines.extend(_stage_summary_block(item))
    if not has_null:
        lines.append(
            f"> ⚠️ 零模型校准未获得结果（{reason}）；λGC 校准与假阳性率无法评估，"
            "以下仅展示独立重算校验等可获得的证据。"
        )
        lines.append("")
    null = metrics.get("null") if isinstance(metrics.get("null"), Mapping) else {}
    sub = [0]

    def next_heading(text: str) -> str:
        sub[0] += 1
        return f"### 7.{sub[0]} {text}"

    if null:
        lines.append(next_heading("零模型校准"))
        lines.append("")
        type1 = null.get("type1") if isinstance(null.get("type1"), Mapping) else {}
        rows = [
            ["重复次数", _fmt(null.get("n_replicates"))],
            ["λGC 均值", _fmt_lambda(null.get("lambda_mean"))],
            ["λGC 标准差", _fmt_lambda(null.get("lambda_sd"))],
            ["α = 0.05 经验假阳性率", _fmt(type1.get("alpha_0.05"))],
            ["α = 1e-3 经验假阳性率", _fmt(type1.get("alpha_1e-3"))],
            ["α = 5e-8 经验假阳性率", _fmt(type1.get("alpha_5e-8"))],
        ]
        lines.extend(_md_table(["指标", "取值"], rows))
        lines.append("")
        try:
            lam_mean = float(null.get("lambda_mean"))
            if 0.97 <= lam_mean <= 1.03:
                lines.append(f"零模型 λGC 均值 {lam_mean:.3f} 落在可接受区间 [0.97, 1.03]，提示检验统计量校准良好。")
            else:
                lines.append(
                    f"⚠️ 零模型 λGC 均值 {lam_mean:.3f} 超出 [0.97, 1.03]，提示 p 值分布存在系统性偏差，"
                    "需检查表型构建、协变量与主成分处理。"
                )
            lines.append("")
        except Exception:
            pass
    power = metrics.get("power") if isinstance(metrics.get("power"), Mapping) else {}
    if power:
        lines.append(next_heading("已知 QTL 功效"))
        lines.append("")
        rows = [["设定遗传力 h²", _fmt(power.get("h2"))], ["因果位点数", _fmt(power.get("n_causal"))]]
        for key in ("alpha_5e-8", "alpha_1e-4", "alpha_5e-3"):
            value = power.get(key)
            if isinstance(value, Mapping):
                rows.append([f"{_label_for(key, _METRIC_LABELS)} 功效均值", _fmt(value.get("mean"))])
                rows.append([f"{_label_for(key, _METRIC_LABELS)} 功效标准差", _fmt(value.get("sd"))])
            elif value is not None:
                rows.append([f"{_label_for(key, _METRIC_LABELS)} 功效", _fmt(value)])
        lines.extend(_md_table(["指标", "取值"], rows))
        lines.append("")
    bias = metrics.get("effect_bias") if isinstance(metrics.get("effect_bias"), Mapping) else {}
    concordance = metrics.get("concordance") if isinstance(metrics.get("concordance"), Mapping) else {}
    if bias or concordance:
        lines.append(next_heading("效应量与排名一致性"))
        lines.append("")
        rows = []
        if bias:
            rows.append(["平均 beta 比值", _fmt(bias.get("mean_beta_ratio"))])
            rows.append(["SE 校准系数", _fmt(bias.get("se_calibration"))])
        if concordance:
            rows.append(["因果位点进入前 100 的比例", _fmt(concordance.get("top_causal_in_top100"))])
            rows.append(["因果位点排名 AUC", _fmt(concordance.get("rank_auc"))])
        lines.extend(_md_table(["指标", "取值"], rows))
        lines.append("")
    if isinstance(validation, Mapping) and validation:
        lines.append(next_heading("独立重算校验"))
        lines.append("")
        rows = [
            ["校验项数", _fmt(validation.get("n_checks"))],
            ["通过", _fmt(validation.get("n_pass"))],
            ["警告", _fmt(validation.get("n_warn"))],
            ["失败", _fmt(validation.get("n_fail"))],
            ["总体判定", _badge(validation.get("overall"))],
        ]
        lines.extend(_md_table(["指标", "取值"], rows))
        lines.append("")
        checks = validation.get("checks")
        if isinstance(checks, list) and checks and all(isinstance(c, Mapping) for c in checks):
            check_rows: list[list[str]] = []
            for check in checks:
                check_rows.append(
                    [
                        _clean_cell(check.get("name") or "—"),
                        _badge(check.get("status")),
                        _clean_cell(check.get("value", "—")),
                        _clean_cell(check.get("threshold", "—")),
                        _clean_cell(check.get("detail", "—"), max_len=120),
                    ]
                )
            lines.extend(_md_table(["校验项", "状态", "实测值", "阈值", "说明"], check_rows))
            lines.append("")
    elif item["status"] != "skipped":
        lines.append("> 未找到 validation.json；独立重算校验结果缺失。")
        lines.append("")
    lines.extend(_render_figure_block(figures.get("calibration", [])))
    return lines


def _section_methods(config: Mapping[str, Any], run_meta: Mapping[str, Any] | None = None) -> list[str]:
    lines = ["## 八、方法学与参数", ""]
    lines.append("### 8.1 方法学依据")
    lines.append("")
    lines.extend(
        [
            "- PLINK 1.9 是全基因组关联与群体遗传分析的标准工具链（Chang et al., 2015, *GigaScience* 4:7；"
            "https://www.cog-genomics.org/plink/1.9/）。",
            "- 质控阈值的推荐范围参考 Anderson et al., 2010, *PLoS One* 5:e11342（样本/位点过滤与 HWE 处理）。",
            "- PCA 作为群体分层校正的标准做法（Price et al., 2006, *Nature Genetics* 38:904；"
            "Patterson et al., 2006, *PLoS Genetics* 2:e190）。",
            "- 亲缘/IBD 估计与重复样本识别使用 PI_HAT（Purcell et al., 2007, *AJHG* 81:559）。",
            "- 多重检验采用 Bonferroni 与 Benjamini–Hochberg FDR；建议性阈值仅供探索。",
            "- ROH/FROH 作为近交代理指标，需注意芯片密度与群体历史对估计的影响（Browning & Browning, 2011）。",
            "- PLINK 1.9 不含线性混合模型（LMM）；强亲缘/强分层数据建议使用 GCTA、GEMMA、BOLT-LMM 等复核"
            "（Yang et al., 2011, *AJHG* 88:76）。",
            "- 仿真校准与独立重算用于验证流程实现，不替代独立队列的重复验证。",
        ]
    )
    lines.append("")
    lines.append("### 8.2 关键参数")
    lines.append("")
    for section_title, dotted in (
        ("质控", "qc"),
        ("群体结构", "popgen"),
        ("GWAS", "gwas"),
        ("校准", "calibration"),
        ("运行时", "runtime"),
    ):
        node = _get_path(config, dotted)
        if not isinstance(node, Mapping) or not node:
            continue
        rows = []
        for key, value in _iter_kv(node):
            full_key = f"{dotted}.{key}"
            rows.append(
                [
                    _clean_cell(full_key),
                    _clean_cell(value),
                    _clean_cell(_param_desc(full_key)),
                ]
            )
        if rows:
            lines.append(f"**{section_title}参数**")
            lines.append("")
            lines.extend(_md_table(["参数", "取值", "说明"], rows))
            lines.append("")
    environment: dict[str, Any] = {}
    if isinstance(run_meta, Mapping):
        for key in ("environment", "env", "versions", "version"):
            value = run_meta.get(key)
            if isinstance(value, Mapping):
                environment.update({str(k): v for k, v in value.items()})
        for key in ("plink_version", "python_version", "python", "os", "platform"):
            value = run_meta.get(key)
            if value is not None and not isinstance(value, Mapping):
                environment[key] = value
    if environment:
        lines.append("### 8.3 运行环境（来自 run.json）")
        lines.append("")
        lines.extend(_md_table(["项目", "取值"], [[_clean_cell(k), _clean_cell(v)] for k, v in environment.items()]))
        lines.append("")
    return lines


def _section_repro(config: Mapping[str, Any], run_dir: Path) -> list[str]:
    lines = ["## 九、复现命令", ""]
    genotype = _get_path(config, "input.genotype") or {}
    target = genotype.get("path") or genotype.get("prefix") if isinstance(genotype, Mapping) else None
    target = target or "<基因型前缀或路径>"
    argv = ["python scripts\\plink_pipeline.py auto", f'--input "{target}"']
    phenotype = _get_path(config, "input.phenotype") or {}
    if isinstance(phenotype, Mapping) and phenotype.get("path"):
        argv.append(f'--pheno "{phenotype["path"]}"')
    covariates = _get_path(config, "input.covariates") or {}
    if isinstance(covariates, Mapping) and covariates.get("path"):
        argv.append(f'--covar "{covariates["path"]}"')
    groups = _get_path(config, "input.groups") or {}
    if isinstance(groups, Mapping) and groups.get("path"):
        argv.append(f'--group "{groups["path"]}"')
    species = _get_path(config, "species.name")
    if species:
        argv.append(f"--species {species}")
    run_id = config.get("run_id") or run_dir.name
    argv.append(f"--run-id {run_id}")
    argv.append(f'--out-root "{run_dir.parent}"')
    cmd_auto = " ".join(argv)
    lines.append("在 skill 根目录（`plink-analysis/`）下执行：")
    lines.append("")
    lines.append("```powershell")
    skill_dir = Path(__file__).resolve().parents[2]
    lines.append(f"cd {skill_dir}")
    lines.append(cmd_auto)
    lines.append("```")
    lines.append("")
    lines.append("对已有 run 目录重跑或补跑阶段：")
    lines.append("")
    lines.append("```powershell")
    lines.append(f'python scripts\\plink_pipeline.py run --run-dir "{run_dir}" --stage all --force')
    lines.append(f'python scripts\\plink_pipeline.py report --run-dir "{run_dir}"')
    lines.append(f'python scripts\\plink_pipeline.py status --run-dir "{run_dir}"')
    lines.append("```")
    lines.append("")
    lines.append(
        "> 复现需保持相同版本的 PLINK（1.9）、Python 依赖与 run.json 中的参数；"
        "随机种子与仿真设置记录在 `run.json` / `04_calibration/metrics.json` 中。"
    )
    lines.append("")
    return lines


def _section_limitations(stages: Mapping[str, Mapping[str, Any]], config: Mapping[str, Any]) -> list[str]:
    lines = ["## 十、局限性", ""]
    items = [
        "统计关联不等于因果：GWAS 显著位点仅提示候选区域，仍需独立群体重复、功能注释与实验验证。",
        "PLINK 1.9 不提供混合模型（LMM）；PCA 协变量只能部分校正群体分层与亲缘结构，强结构群体中可能仍有残余膨胀。",
        "芯片标记的 ascertainment bias、未填补基因型与低 MAF 位点会限制稀有变异与因果位点的检出。",
        "多重检验 Bonferroni 较为保守；FDR 与建议性阈值提高了发现数但增加了假阳性风险，需在结论中区分。",
        "样本量决定检验功效；小样本下的效应量估计与显著性均不稳定。",
        "质控阈值、LD 剪枝参数与主成分数量会影响结果；本报告固定记录参数，但未做全参数敏感性分析。",
        "表型测量误差、缺失机制与协变量编码依赖上游数据质量，本流程无法纠正系统性表型偏差。",
        "仿真校准验证的是流程实现与统计性质，不能替代独立队列或前瞻性验证；仿真中的遗传结构是对真实群体的简化。",
        "FROH、PI_HAT 等指标受芯片密度与群体历史影响，跨研究比较时应谨慎。",
    ]
    qc = stages["qc"]["metrics"]
    try:
        n = int(qc.get("n_samples_out") or 0)
        if 0 < n < 200:
            items.insert(0, f"质控后样本量仅 {n}，统计功效有限，显著位点与效应量需谨慎解释。")
    except Exception:
        pass
    cal_item = stages.get("calibration")
    if cal_item:
        cal_status = _norm_status(cal_item.get("status"))
        has_null, _has_power, _n_checks, _validation = _calibration_evidence(cal_item)
        if cal_status in ("skipped", "unknown") or not has_null:
            items.insert(
                0,
                "本次运行未完成仿真校准：p 值校准（λGC、假阳性率）与统计效力未经独立验证，"
                "显著性结论需谨慎解读，建议补充仿真校准或独立队列验证。",
            )
    lines.extend(f"- {item}" for item in items)
    lines.append("")
    return lines


def _section_artifacts(stages: Mapping[str, Mapping[str, Any]], run_dir: Path, reports_dir: Path) -> list[str]:
    lines = ["## 十一、文件与产物清单", ""]
    rows: list[list[str]] = []
    for name in ("run.json", "state.json"):
        path = run_dir / name
        if path.exists():
            rows.append([name, _human_size(path.stat().st_size), "运行配置 / 状态"])
    for stage in STAGE_ORDER:
        item = stages[stage]
        stage_dir = item["dir"]
        if not isinstance(stage_dir, Path) or not stage_dir.exists():
            continue
        for path in sorted(stage_dir.rglob("*")):
            if path.is_file():
                rel = _rel_path(path, run_dir)
                rows.append([rel, _human_size(path.stat().st_size), _artifact_description(rel)])
    for name, description in (
        ("final_report.md", "最终报告（Markdown）"),
        ("final_report.html", "最终报告（自包含 HTML）"),
        ("summary.json", "机器可读摘要"),
    ):
        path = reports_dir / name
        if path.exists():
            rows.append([f"{REPORTS_DIR}/{name}", _human_size(path.stat().st_size), description])
    if rows:
        lines.extend(_md_table(["文件", "大小", "说明"], rows[:80]))
        if len(rows) > 80:
            lines.append(f"> 另有 {len(rows) - 80} 个文件未列出。")
    else:
        lines.append("运行目录中暂无文件。")
    lines.append("")
    lines.append("**阶段报告索引**")
    lines.append("")
    for stage in STAGE_ORDER:
        item = stages[stage]
        if item["report_exists"]:
            rel = _rel_path(item["report"], reports_dir)
            lines.append(f"- [{_STAGE_META[stage]['name']}报告]({rel})（{_badge(item['status'])}）")
    lines.append("")
    return lines


def _build_final_markdown(
    run_dir: Path,
    config: Mapping[str, Any],
    stages: Mapping[str, Mapping[str, Any]],
    run_meta: Mapping[str, Any],
    figures: Mapping[str, list[tuple[str, str]]],
) -> str:
    overall = _overall_status(stages)
    project = _project_title(config, run_dir)
    lines: list[str] = []
    lines.append(f"# {project} — PLINK 全流程分析最终报告")
    lines.append("")
    lines.append(
        f"> **运行编号**：`{config.get('run_id', run_dir.name)}` ｜ **生成时间**：{_now_str()} ｜ "
        f"**综合判定**：{_badge(overall)} ｜ **物种**：{_species_label(config, stages['intake']['metrics'])}"
    )
    if config.get("created_at"):
        lines.append(f"> **运行创建时间**：{_fmt(config.get('created_at'))} ｜ **报告版本**：1（plink-analysis skill）")
    lines.append("")
    lines.extend(_section_executive_summary(stages, config, overall))
    lines.append("## 二、阶段状态总览")
    lines.append("")
    lines.extend(_stage_status_table(stages, run_dir))
    lines.append("")
    lines.extend(_section_data_overview(stages, config, figures))
    lines.extend(_section_qc(stages, figures))
    lines.extend(_section_popgen(stages, figures))
    lines.extend(_section_gwas(stages, figures))
    lines.extend(_section_calibration(stages, figures))
    lines.extend(_section_methods(config, run_meta))
    lines.extend(_section_repro(config, run_dir))
    lines.extend(_section_limitations(stages, config))
    lines.extend(_section_artifacts(stages, run_dir, run_dir / REPORTS_DIR))
    lines.append("---")
    lines.append("")
    lines.append(
        "> 本报告由 [plink-analysis skill](https://www.cog-genomics.org/plink/1.9/) 自动生成，"
        "所有数值以各阶段 `metrics.json` 与 `tables/` 中的原始结果为准。"
        "阶段报告的详细诊断信息请查阅对应 `report.md`。"
    )
    return _renumber_figures("\n".join(lines) + "\n")


def _build_summary(
    run_dir: Path,
    config: Mapping[str, Any],
    stages: Mapping[str, Mapping[str, Any]],
    figures: Mapping[str, list[tuple[str, str]]],
    outputs: Mapping[str, Path],
) -> dict:
    intake = stages["intake"]["metrics"]
    qc = stages["qc"]["metrics"]
    popgen = stages["popgen"]["metrics"]
    gwas = stages["gwas"]["metrics"]
    calibration = stages["calibration"]["metrics"]
    validation = _merge_validation(stages["calibration"])
    pca = popgen.get("pca") if isinstance(popgen.get("pca"), Mapping) else {}
    rel = popgen.get("relatedness") if isinstance(popgen.get("relatedness"), Mapping) else {}
    null = calibration.get("null") if isinstance(calibration.get("null"), Mapping) else {}
    power = calibration.get("power") if isinstance(calibration.get("power"), Mapping) else {}
    cal_item = stages["calibration"]
    cal_status = _norm_status(cal_item.get("status"))
    has_null, has_power, n_checks_evidence, validation = _calibration_evidence(cal_item)
    calibration_reason = (
        _short_reason(_skip_reason(cal_item, "仿真校准未执行或未获得有效结果"))
        if (cal_status in ("skipped", "unknown") or not has_null or n_checks_evidence == 0)
        else None
    )
    calibration_skip_reason = calibration_reason if cal_status in ("skipped", "unknown") else None
    calibration_null_missing_reason = (
        calibration_reason
        if (cal_status not in ("skipped", "unknown") and not has_null)
        else None
    )
    figure_list = [rel_path for stage in STAGE_ORDER for _caption, rel_path in figures.get(stage, [])]

    def as_list(value: Any) -> list[Any]:
        return list(value) if isinstance(value, (list, tuple)) else []

    summary = {
        "schema_version": 1,
        "run_id": config.get("run_id") or run_dir.name,
        "generated_at": _now_str(),
        "overall_status": _overall_status(stages),
        "stages": {
            stage: {
                "status": stages[stage]["status"],
                "report": _rel_path(stages[stage]["report"], run_dir) if stages[stage]["report_exists"] else None,
                "metrics": _rel_path(stages[stage]["dir"] / "metrics.json", run_dir)
                if isinstance(stages[stage]["dir"], Path) and (stages[stage]["dir"] / "metrics.json").exists()
                else None,
                "elapsed_s": stages[stage].get("state", {}).get("elapsed_s"),
            }
            for stage in STAGE_ORDER
        },
        "data": {
            "genotype_type": intake.get("genotype_type"),
            "n_samples": intake.get("n_samples"),
            "n_variants": intake.get("n_variants"),
            "chromosomes": as_list(intake.get("chromosomes")),
            "species": intake.get("species_inferred"),
            "phenotype": intake.get("phenotype"),
            "covariates": intake.get("covariates"),
            "groups": intake.get("groups"),
        },
        "qc": {
            "n_samples_in": qc.get("n_samples_in"),
            "n_samples_out": qc.get("n_samples_out"),
            "n_variants_in": qc.get("n_variants_in"),
            "n_variants_out": qc.get("n_variants_out"),
            "overall": qc.get("overall"),
            "gates": as_list(qc.get("gates")),
            "removed": qc.get("removed"),
        },
        "popgen": {
            "pca": pca,
            "relatedness": rel,
            "roh": popgen.get("roh"),
            "fst": popgen.get("fst"),
            "ld": popgen.get("ld"),
            "overall": popgen.get("overall"),
        },
        "gwas": {
            "model": gwas.get("model"),
            "phenotype_type": gwas.get("phenotype_type"),
            "n_samples": gwas.get("n_samples"),
            "n_tested": gwas.get("n_tested"),
            "covariates": gwas.get("covariates"),
            "lambda_gc": gwas.get("lambda_gc"),
            "thresholds": gwas.get("thresholds"),
            "n_significant": gwas.get("n_significant"),
            "n_fdr05": gwas.get("n_fdr05"),
            "n_suggestive": gwas.get("n_suggestive"),
            "lead_snps": as_list(gwas.get("lead_snps")),
            "p_distribution": gwas.get("p_distribution"),
            "min_p_text": _min_p_display(gwas.get("p_distribution")),
            "overall": gwas.get("overall"),
        },
        "calibration": {
            "null": null,
            "power": power,
            "effect_bias": calibration.get("effect_bias"),
            "concordance": calibration.get("concordance"),
            "evidence": {"null": has_null, "power": has_power, "n_checks": n_checks_evidence},
            "skip_reason": calibration_skip_reason,
            "null_missing_reason": calibration_null_missing_reason,
            "overall": calibration.get("overall"),
        },
        "validation": {
            "n_checks": validation.get("n_checks"),
            "n_pass": validation.get("n_pass"),
            "n_warn": validation.get("n_warn"),
            "n_fail": validation.get("n_fail"),
            "overall": validation.get("overall"),
        },
        "warnings": {
            stage: _as_list(stages[stage]["metrics"].get("warnings"))
            for stage in STAGE_ORDER
            if _as_list(stages[stage]["metrics"].get("warnings"))
        },
        "errors": {
            stage: _as_list(stages[stage]["metrics"].get("errors"))
            for stage in STAGE_ORDER
            if _as_list(stages[stage]["metrics"].get("errors"))
        },
        "figures": figure_list,
        "outputs": {key: _rel_path(path, run_dir) for key, path in outputs.items()},
    }
    return summary


def final_report(run_dir: Path, config: dict, state: dict | None = None) -> dict:
    """Write final_report.md/.html and summary.json under ``05_reports``.

    Returns ``{"md": str, "html": str, "summary": str}`` (paths are strings so the
    dict can be printed or JSON-encoded safely by the CLI).
    """
    run_path = Path(run_dir)
    cfg: Mapping[str, Any] = config if isinstance(config, Mapping) else {}
    stages = _collect_stages(run_path, state, cfg)
    run_meta = _read_json(run_path / "run.json", {})
    run_meta = run_meta if isinstance(run_meta, Mapping) else {}
    reports_dir = run_path / REPORTS_DIR
    reports_dir.mkdir(parents=True, exist_ok=True)
    figures = _collect_figures_for_final(stages, reports_dir / "figures")

    markdown = _build_final_markdown(run_path, cfg, stages, run_meta, figures)
    md_path = reports_dir / "final_report.md"
    _write_text(md_path, markdown)

    html_path = reports_dir / "final_report.html"
    title = f"{_project_title(cfg, run_path)} — PLINK 全流程分析报告"
    render_markdown_to_html(markdown, html_path, title=title, base_dir=reports_dir)

    summary_path = reports_dir / "summary.json"
    outputs = {"md": md_path, "html": html_path, "summary": summary_path}
    summary = _build_summary(run_path, cfg, stages, figures, outputs)
    _write_json(summary_path, summary)

    return {"md": str(md_path), "html": str(html_path), "summary": str(summary_path)}


# ---------------------------------------------------------------------------
# Markdown -> self-contained HTML
# ---------------------------------------------------------------------------
_ASSETS_CSS = Path(__file__).resolve().parents[2] / "assets" / "report.css"

_FALLBACK_CSS = """
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Microsoft YaHei",sans-serif;
margin:0;padding:32px;background:#fff;color:#1f2328;line-height:1.7}
main.report{max-width:1100px;margin:0 auto}
h1,h2,h3{line-height:1.35;color:#123}
table{border-collapse:collapse;width:100%;margin:12px 0}
th,td{border:1px solid #d0d7de;padding:6px 10px;text-align:left;vertical-align:top}
th{background:#f2f4f7}
img{max-width:100%;height:auto}
code,pre{font-family:ui-monospace,Consolas,monospace}
.st-pass{background:#e6f4ea;color:#1e7b34;padding:1px 6px;border-radius:4px}
.st-warn{background:#fff4e5;color:#92400e;padding:1px 6px;border-radius:4px}
.st-fail{background:#fdecea;color:#a50e0e;padding:1px 6px;border-radius:4px}
.st-skip{background:#eef2f6;color:#4b5563;padding:1px 6px;border-radius:4px}
"""


def _load_css() -> str:
    try:
        if _ASSETS_CSS.exists():
            text = _ASSETS_CSS.read_text(encoding="utf-8")
            if text.strip():
                return text
    except Exception:
        pass
    return _FALLBACK_CSS


def _decorate_status(text: str) -> str:
    pattern = re.compile(r"((?:✅|⚠️|❌|⏭️|⏳)\s*)?\b(PASS|WARN|FAIL|SKIPPED|ERROR|RUNNING)\b")

    def repl(match: re.Match[str]) -> str:
        prefix = match.group(1) or ""
        word = match.group(2)
        css = {
            "PASS": "st-pass",
            "WARN": "st-warn",
            "FAIL": "st-fail",
            "ERROR": "st-fail",
            "SKIPPED": "st-skip",
            "RUNNING": "st-running",
        }.get(word, "st-unknown")
        return f'<span class="st {css}">{prefix}{word}</span>'

    return pattern.sub(repl, text)


def _resolve_image(path_text: str, base_dir: Path) -> Path | None:
    cleaned = path_text.strip().strip("<>").strip()
    if not cleaned or cleaned.startswith(("http://", "https://", "data:")):
        return None
    cleaned = _html.unescape(cleaned)
    candidate = Path(cleaned)
    options: list[Path] = []
    if candidate.is_absolute():
        options.append(candidate)
    else:
        options.append(base_dir / candidate)
        options.append(Path(cleaned))
    for option in options:
        try:
            if option.is_file():
                return option
        except OSError:
            continue
    # last resort: search below the base directory for the file name
    try:
        for found in base_dir.rglob(candidate.name):
            if found.is_file():
                return found
    except Exception:
        pass
    return None


def _image_html(path_text: str, alt: str, base_dir: Path) -> str:
    alt_escaped = _html.escape(alt or "figure", quote=True)
    stripped = path_text.strip()
    if stripped.startswith("data:"):
        return f'<img src="{_html.escape(stripped, quote=True)}" alt="{alt_escaped}">'
    if stripped.startswith(("http://", "https://")):
        return f'<img src="{_html.escape(stripped, quote=True)}" alt="{alt_escaped}">'
    resolved = _resolve_image(stripped, base_dir)
    if resolved is None:
        return f'<span class="fig-missing">[图片缺失：{_html.escape(stripped)}]</span>'
    mime = mimetypes.guess_type(resolved.name)[0] or "image/png"
    try:
        encoded = base64.b64encode(resolved.read_bytes()).decode("ascii")
    except Exception:
        return f'<span class="fig-missing">[图片无法读取：{_html.escape(stripped)}]</span>'
    src = f"data:{mime};base64,{encoded}"
    caption = f"<figcaption>{_html.escape(alt)}</figcaption>" if alt else ""
    return f'<figure><img src="{src}" alt="{alt_escaped}">{caption}</figure>'


def _inline(text: str, base_dir: Path) -> str:
    escaped = _html.escape(text, quote=False)
    tokens: dict[str, str] = {}
    counter = [0]

    def keep(html_fragment: str) -> str:
        counter[0] += 1
        token = f"\x00TOK{counter[0]}\x00"
        tokens[token] = html_fragment
        return token

    # inline code (protect before other inline formatting)
    escaped = re.sub(
        r"`([^`\n]+)`",
        lambda m: keep(f"<code>{m.group(1)}</code>"),
        escaped,
    )
    # images
    escaped = re.sub(
        r"!\[([^\]]*)\]\(([^)]+)\)",
        lambda m: keep(_image_html(m.group(2), m.group(1), base_dir)),
        escaped,
    )
    # links
    escaped = re.sub(
        r"\[([^\]]+)\]\(([^)]+)\)",
        lambda m: keep(
            f'<a href="{_html.escape(_html.unescape(m.group(2).strip()), quote=True)}">{m.group(1)}</a>'
        ),
        escaped,
    )
    # bold then italic
    escaped = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", escaped)
    escaped = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<em>\1</em>", escaped)
    escaped = re.sub(r"~~(.+?)~~", r"<del>\1</del>", escaped)
    # status colouring (tokens are unaffected; base64 never appears here)
    escaped = _decorate_status(escaped)
    for token, fragment in tokens.items():
        escaped = escaped.replace(token, fragment)
    return escaped


def _is_table_sep(line: str) -> bool:
    return bool(re.match(r"^\s*\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)+\|?\s*$", line))


def _split_table_row(line: str) -> list[str]:
    text = line.strip()
    if text.startswith("|"):
        text = text[1:]
    if text.endswith("|") and not text.endswith("\\|"):
        text = text[:-1]
    cells = re.split(r"(?<!\\)\|", text)
    return [cell.replace("\\|", "|").strip() for cell in cells]


def _render_table_html(lines: Sequence[str], base_dir: Path) -> str:
    header = _split_table_row(lines[0])
    aligns: list[str] = []
    if len(lines) > 1:
        for cell in _split_table_row(lines[1]):
            value = cell.strip()
            if value.startswith(":") and value.endswith(":"):
                aligns.append("center")
            elif value.endswith(":"):
                aligns.append("right")
            else:
                aligns.append("left")
    body = [_split_table_row(line) for line in lines[2:]]
    head_html = "".join(f"<th>{_inline(cell, base_dir)}</th>" for cell in header)
    rows_html: list[str] = []
    for row in body:
        cells_html: list[str] = []
        for index, cell in enumerate(row):
            align = aligns[index] if index < len(aligns) else "left"
            cells_html.append(f'<td class="align-{align}">{_inline(cell, base_dir)}</td>')
        rows_html.append("<tr>" + "".join(cells_html) + "</tr>")
    return (
        '<div class="table-wrap"><table><thead><tr>'
        + head_html
        + "</tr></thead><tbody>"
        + "".join(rows_html)
        + "</tbody></table></div>"
    )


def _is_block_start(line: str, next_line: str | None) -> bool:
    if re.match(r"^\s*```", line):
        return True
    if re.match(r"^#{1,6}\s+", line):
        return True
    if re.match(r"^\s*>", line):
        return True
    if re.match(r"^\s*([-*+]|\d+[.)])\s+", line):
        return True
    if re.match(r"^\s*([-*_])(\s*\1){2,}\s*$", line):
        return True
    if "|" in line and next_line is not None and _is_table_sep(next_line):
        return True
    return False


def _render_list_html(lines: Sequence[str], start: int, base_dir: Path) -> tuple[str, int]:
    ordered = bool(re.match(r"^\s*\d+[.)]\s+", lines[start]))
    tag = "ol" if ordered else "ul"
    items: list[str] = []
    index = start
    while index < len(lines):
        match = re.match(r"^(\s*)([-*+]|\d+[.)])\s+(.*)$", lines[index])
        if not match:
            break
        content = match.group(3)
        index += 1
        while (
            index < len(lines)
            and lines[index].strip()
            and not re.match(r"^(\s*)([-*+]|\d+[.)])\s+", lines[index])
            and lines[index].startswith((" ", "\t"))
        ):
            content += " " + lines[index].strip()
            index += 1
        items.append(f"<li>{_inline(content, base_dir)}</li>")
    return f"<{tag}>" + "".join(items) + f"</{tag}>", index


def _markdown_body(markdown: str, base_dir: Path) -> str:
    lines = markdown.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    output: list[str] = []
    index = 0
    total = len(lines)
    while index < total:
        line = lines[index]
        if not line.strip():
            index += 1
            continue
        next_line = lines[index + 1] if index + 1 < total else None
        # fenced code block
        fence = re.match(r"^\s*```(\w*)\s*$", line)
        if fence:
            language = fence.group(1)
            index += 1
            buffer: list[str] = []
            while index < total and not re.match(r"^\s*```", lines[index]):
                buffer.append(lines[index])
                index += 1
            index += 1  # closing fence
            code = _html.escape("\n".join(buffer))
            css_class = f' class="language-{language}"' if language else ""
            output.append(f"<pre><code{css_class}>{code}</code></pre>")
            continue
        # heading
        heading = re.match(r"^(#{1,6})\s+(.*?)\s*#*\s*$", line)
        if heading:
            level = len(heading.group(1))
            output.append(f"<h{level}>{_inline(heading.group(2), base_dir)}</h{level}>")
            index += 1
            continue
        # horizontal rule
        if re.match(r"^\s*([-*_])(\s*\1){2,}\s*$", line):
            output.append("<hr>")
            index += 1
            continue
        # table
        if "|" in line and next_line is not None and _is_table_sep(next_line):
            table_lines = [line, next_line]
            index += 2
            while index < total and "|" in lines[index] and lines[index].strip():
                table_lines.append(lines[index])
                index += 1
            output.append(_render_table_html(table_lines, base_dir))
            continue
        # blockquote
        if line.lstrip().startswith(">"):
            buffer = []
            while index < total and lines[index].lstrip().startswith(">"):
                buffer.append(re.sub(r"^\s*>\s?", "", lines[index]))
                index += 1
            inner = "<br>".join(_inline(item, base_dir) for item in buffer if item.strip())
            output.append(f"<blockquote>{inner}</blockquote>")
            continue
        # list
        if re.match(r"^\s*([-*+]|\d+[.)])\s+", line):
            list_html, index = _render_list_html(lines, index, base_dir)
            output.append(list_html)
            continue
        # paragraph
        buffer = []
        while index < total:
            current = lines[index]
            if not current.strip():
                break
            following = lines[index + 1] if index + 1 < total else None
            if buffer and _is_block_start(current, following):
                break
            if not buffer and _is_block_start(current, following):
                break
            buffer.append(current.strip())
            index += 1
        if buffer:
            if len(buffer) == 1 and re.fullmatch(r"!\[[^\]]*\]\([^)]+\)", buffer[0].strip()):
                output.append(_inline(buffer[0], base_dir))
            else:
                output.append(f"<p>{_inline(' '.join(buffer), base_dir)}</p>")
        else:  # pragma: no cover - safety net against infinite loops
            index += 1
    return "\n".join(output)


def render_markdown_to_html(
    md_text: str,
    out_path: Path,
    *,
    title: str = "PLINK 分析报告",
    base_dir: Path | str | None = None,
) -> Path:
    """Render *md_text* to a self-contained HTML file and return its path.

    Images with relative paths are resolved against *base_dir* (defaults to the
    directory of *out_path*) and embedded as base64 ``data:`` URIs.  The CSS from
    ``assets/report.css`` is inlined, so the resulting file opens directly in a
    browser without any external resources.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    base = Path(base_dir) if base_dir is not None else out.parent
    body = _markdown_body(str(md_text or ""), base)
    css = _load_css()
    document = (
        "<!DOCTYPE html>\n"
        '<html lang="zh-CN">\n'
        "<head>\n"
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{_html.escape(title or 'PLINK 分析报告')}</title>\n"
        "<style>\n"
        f"{css}\n"
        "</style>\n"
        "</head>\n"
        "<body>\n"
        '<main class="report">\n'
        f"{body}\n"
        "</main>\n"
        f'<footer class="report-footer">本报告由 plink-analysis skill 自动生成 · {_now_str()}</footer>\n'
        "</body>\n"
        "</html>\n"
    )
    out.write_text(document, encoding="utf-8", newline="\n")
    return out


__all__ = ["stage_report", "final_report", "render_markdown_to_html", "STAGE_DIRS", "STAGE_ORDER"]
