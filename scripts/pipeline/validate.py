#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""独立重算校验（DESIGN 第 8 节 validation.json；被 calibration stage 调用）。

设计原则
--------
validate 不看 stage 的自述数字，而是**从产物文件独立重算**再比对：

* 计数：``.fam``/``.bim``/``.bed`` 字节数 vs PLINK 日志 vs ``metrics.json``；
* 频率/缺失：``--freq``/``--missing`` vs 从 ``--recode A`` 剂量逐位点重算；
* 表型：编码（1/2 vs 0/1）、缺失、与 QC 后样本集合的匹配；
* PCA：用 ``.eigenvec``/``.eigenval`` 与 pruned 基因型独立做一次 PCA，
  比对特征向量方向与特征值尺度；并用「特征值之和 = 样本数 N」检验
  ``var_explained`` 是否算错；
* 亲缘：``.genome`` PI_HAT 范围/计数、``.rel`` 对称性与对角；
* GWAS：随机 20 个 SNP，用 statsmodels OLS/Logit 在同样本、同协变量下重算
  beta/SE，与 PLINK ``.assoc.linear/.logistic`` 比对（容差 1e-4）；
* 多重检验：从 P 值独立重算 λGC，与 PLINK ``--adjust`` 日志值和 metrics 比对；
* HWE：QC 后位点 p<1e-6 比例是否合理。

调用约定
--------
``run_validation(ctx_or_run_dir) -> dict``

* pipeline 内（calibration_runner 已这样调用）::

      from .validate import run_validation
      validation = run_validation(ctx)          # ctx: PipelineContext

* 单独对一个 run 目录执行（自动读 ``run.json``/``state.json`` 并定位 PLINK）::

      python -m pipeline.validate D:/plink_runs/demo1
      # 或
      from pipeline.validate import run_validation
      run_validation(Path("D:/plink_runs/demo1"))

写 ``<run_dir>/04_calibration/validation.json`` 与
``<run_dir>/04_calibration/tables/validation_checks.csv``，返回::

    {"checks": [{"name", "category", "status", "detail", "value", "threshold"}, ...],
     "overall": "pass|warn|fail",
     "n_checks": N, "n_pass": ..., "n_warn": ..., "n_fail": ...,
     "generated_at": "...", "run_id": "...", "artifacts": {...}}

每项 ``status`` ∈ pass|warn|fail；``overall`` = 有 fail 则 fail，有 warn 则 warn。
产物缺失时对应检查标记为 warn（detail 说明缺失路径），不会中断整个校验。
容差可用 ``run.json`` 的 ``validation`` 段覆盖（见 ``DEFAULT_TOL``）。
"""
from __future__ import annotations

import datetime as _dt
import logging
import math
import os
import re
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:  # pipeline 包内导入（calibration_runner 用法）
    from . import utils
except ImportError:  # pragma: no cover - 兼容 scripts/pipeline 直接执行
    _here = Path(__file__).resolve().parent
    sys.path.insert(0, str(_here.parent))
    sys.path.insert(0, str(_here))
    import utils  # type: ignore

# ---------------------------------------------------------------------------
# 常量 / 默认容差
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


DEFAULT_PLINK = _default_plink()
ALPHA_GW = 5e-8
VALID_ALLELES = {"A", "C", "G", "T", "0"}

DEFAULT_TOL: dict[str, Any] = {
    "maf": 1e-5,                # MAF 重算容差（DESIGN 第 14 节）
    "missing": 1e-6,            # 缺失率重算容差
    "beta": 1e-4,               # 线性模型 beta/SE 与 statsmodels 对照容差（DESIGN 第 14 节）
    "beta_logistic": 1e-3,      # logistic：PLINK 打印 OR/SE 只保留 4-5 位有效数字且迭代提前停止
    "lambda": 5e-3,             # λGC 相对容差（.assoc P 值仅 4-6 位有效数字）
    "pca_corr": 0.95,           # PCA 特征向量方向相关下限（全量剪枝位点时）
    "pca_corr_subsampled": 0.80,  # 位点过多而抽样时的下限
    "eigen_scale": 0.15,        # PCA 特征值尺度相对容差
    "hwe_warn": 0.005,          # QC 后 HWE p<1e-6 比例 warn 阈值
    "hwe_fail": 0.02,           # fail 阈值
    "n_beta_snps": 20,          # beta/SE 对照随机 SNP 数
    "n_maf_snps": 100,          # MAF 对照随机 SNP 数
    "n_missing_snps": 2000,     # 缺失率对照随机 SNP 数
    "pca_max_values": 40_000_000,   # PCA 重算剂量矩阵最大元素数（内存保护）
    "pca_max_flops": 5_000_000_000,  # PCA Gram 矩阵乘法预算（N^2*m）
    "seed": 20261023,
}


# ---------------------------------------------------------------------------
# 目标（run 目录 / ctx）解析
# ---------------------------------------------------------------------------
@dataclass
class _Target:
    run_dir: Path
    config: dict[str, Any] = field(default_factory=dict)
    state: dict[str, Any] = field(default_factory=dict)
    plink: Any = None
    logger: logging.Logger = field(default_factory=lambda: logging.getLogger("validate"))
    stage_dir: Path | None = None
    metrics: dict[str, Any] = field(default_factory=dict)
    artifacts: dict[str, Any] = field(default_factory=dict)
    tol: dict[str, Any] = field(default_factory=lambda: dict(DEFAULT_TOL))

    @property
    def cal_dir(self) -> Path:
        return Path(self.stage_dir) if self.stage_dir else (self.run_dir / "04_calibration")

    @property
    def data_dir(self) -> Path:
        p = self.cal_dir / "data"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def tables_dir(self) -> Path:
        p = self.cal_dir / "tables"
        p.mkdir(parents=True, exist_ok=True)
        return p


def _read_json(path: Path) -> dict[str, Any]:
    try:
        obj = utils.load_json(path, default={})
        return obj if isinstance(obj, dict) else {}
    except Exception:
        return {}


def _prefix_of(path: str | Path | None) -> Path | None:
    """把 ``*.bed/.bim/.fam`` 文件路径规范成 PLINK 前缀（给 ``--bfile`` 用）。

    state.json 里的 ``artifacts.bfile`` 可能存前缀（``.../clean``），glob 兜底时
    可能是 ``.../clean.bed``；两者都统一成前缀，避免 ``--bfile clean.bed`` 去找
    ``clean.bed.bed`` 的 ``Error: Failed to open``。
    """
    if path is None:
        return None
    p = Path(str(path))
    if p.suffix.lower() in {".bed", ".bim", ".fam"}:
        p = p.with_suffix("")
    return p


def _find_plink_exe(t: _Target) -> str | None:
    runtime = t.config.get("runtime") or {}
    for cand in (
        runtime.get("plink_exe"),
        (t.config.get("_derived") or {}).get("plink_exe"),
        os.environ.get("PLINK_EXE"),
        os.environ.get("PLINK"),
        DEFAULT_PLINK,
        shutil.which("plink"),
    ):
        if cand and Path(str(cand)).is_file():
            return str(cand)
    return None


class _MiniPlink:
    """没有 pipeline.utils 时的兜底 PLINK 封装（独立校验 / 单元测试用）。"""

    def __init__(self, exe: str, run_dir: Path, logger: logging.Logger):
        self.exe = exe
        self.run_dir = Path(run_dir)
        self.logger = logger
        self.version = "unknown"

    def run(self, args, *, stage: str = "calibration", name: str = "cmd",
            out_prefix=None, expect=(), timeout_s=None, check: bool = True):
        import subprocess
        import time

        argv = [self.exe, *[str(a) for a in args]]
        if out_prefix is not None:
            out_prefix = Path(out_prefix)
            out_prefix.parent.mkdir(parents=True, exist_ok=True)
            argv += ["--out", str(out_prefix)]
        started = time.time()
        proc = subprocess.run(argv, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", cwd=str(self.run_dir))
        log_path = Path(str(out_prefix) + ".log") if out_prefix else None
        log_text = (log_path.read_text(encoding="utf-8", errors="replace")
                    if log_path and log_path.exists() else (proc.stdout or "") + (proc.stderr or ""))
        parsed = utils.parse_plink_log(log_text)
        errors = list(parsed.get("errors", []))
        if proc.returncode != 0 and not errors:
            errors = [f"PLINK 返回码 {proc.returncode}"]
        missing = [str(Path(str(out_prefix) + s)) for s in expect
                   if out_prefix and not Path(str(out_prefix) + s).exists()]
        if missing:
            errors.append("缺少预期输出: " + ", ".join(missing))
        if check and errors:
            raise RuntimeError("; ".join(errors[:2]))
        return utils.PlinkResult(args=argv,
                                 out_prefix=str(out_prefix) if out_prefix else None,
                                 returncode=proc.returncode,
                                 elapsed_s=time.time() - started,
                                 log_path=str(log_path) if log_path else None,
                                 log_text=log_text, counts=parsed.get("counts", {}),
                                 warnings=parsed.get("warnings", []), errors=errors,
                                 stdout=proc.stdout or "")


def _resolve(obj: Any, logger: logging.Logger | None = None) -> _Target:
    if hasattr(obj, "run_dir") and not isinstance(obj, (str, Path)):
        logger = getattr(obj, "logger", None) or logger or logging.getLogger("validate")
        t = _Target(
            run_dir=Path(str(obj.run_dir)),
            config=dict(getattr(obj, "config", {}) or {}),
            state=dict(getattr(obj, "state", {}) or {}),
            plink=getattr(obj, "plink", None),
            logger=logger,
            stage_dir=Path(str(obj.stage_dir)) if getattr(obj, "stage_dir", None) else None,
        )
    else:
        run_dir = Path(str(obj))
        if not run_dir.is_absolute():
            run_dir = run_dir.resolve()
        logger = logger or logging.getLogger("validate")
        t = _Target(run_dir=run_dir, logger=logger)
        t.config = _read_json(run_dir / "run.json")
        t.state = _read_json(run_dir / "state.json")
    t.metrics = {
        "intake": _read_json(t.run_dir / "00_intake" / "metrics.json"),
        "qc": _read_json(t.run_dir / "01_qc" / "metrics.json"),
        "popgen": _read_json(t.run_dir / "02_popgen" / "metrics.json"),
        "gwas": _read_json(t.run_dir / "03_gwas" / "metrics.json"),
    }
    tol = dict(DEFAULT_TOL)
    override = t.config.get("validation") or {}
    if isinstance(override, dict):
        for key, val in override.items():
            if key in tol:
                tol[key] = val
    t.tol = tol
    return t


# ---------------------------------------------------------------------------
# 产物发现（state.artifacts 优先，glob 兜底）
# ---------------------------------------------------------------------------
def _stage_art(t: _Target, stage: str, key: str) -> Path | None:
    try:
        val = (((t.state.get("stages") or {}).get(stage) or {})
               .get("artifacts") or {}).get(key)
        if val and Path(str(val)).exists():
            return Path(str(val))
    except Exception:
        pass
    return None


def _rel(t: _Target, path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return str(Path(path).relative_to(t.run_dir)).replace("\\", "/")
    except ValueError:
        return str(path)


def discover_artifacts(t: _Target) -> dict[str, Any]:
    """定位校验所需全部产物；返回 {key: {"path": Path|None, "source": str}}。"""
    rd = t.run_dir
    found: dict[str, Any] = {}

    def put(key: str, path: Path | None, source: str) -> None:
        if path is not None:
            path = Path(path)
            if not path.exists():
                path = None
        if key not in found or (found[key]["path"] is None and path is not None):
            found[key] = {"path": path, "source": source}

    def first(patterns: list[str]) -> Path | None:
        for pat in patterns:
            hits = sorted(rd.glob(pat))
            if hits:
                return hits[0]
        return None

    put("clean", _stage_art(t, "qc", "bfile"), "state.qc.artifacts.bfile")
    if found.get("clean", {}).get("path") is None:
        put("clean", first(["01_qc/data/clean.bed"]), "01_qc/data/clean")
    put("base", _stage_art(t, "intake", "base_prefix"), "state.intake.artifacts.base_prefix")
    if found.get("base", {}).get("path") is None:
        put("base", first(["00_intake/data/base.bed"]), "00_intake/data/base")
    put("pruned", _stage_art(t, "qc", "pruned"), "state.qc.artifacts.pruned")
    if found.get("pruned", {}).get("path") is None:
        put("pruned", first(["01_qc/data/clean_pruned.prune.in",
                             "01_qc/data/pruned.prune.in",
                             "01_qc/data/*.prune.in"]),
            "glob 01_qc/data/*.prune.in")
    put("pheno", _stage_art(t, "intake", "pheno"), "state.intake.artifacts.pheno")
    if found.get("pheno", {}).get("path") is None:
        put("pheno", first(["00_intake/data/pheno_clean.txt",
                            "00_intake/data/pheno*.txt"]),
            "glob 00_intake/data/pheno_clean.txt")
    if found.get("pheno", {}).get("path") is None:
        raw = ((t.metrics.get("intake") or {}).get("phenotype") or {}).get("file")
        if raw:
            put("pheno", Path(str(raw)), "intake.metrics.phenotype.file")
    put("covar", _stage_art(t, "gwas", "covar"), "state.gwas.artifacts.covar")
    if found.get("covar", {}).get("path") is None:
        put("covar", first(["03_gwas/data/covar.txt"]), "03_gwas/data/covar.txt")
    put("assoc", _stage_art(t, "gwas", "assoc"), "state.gwas.artifacts.assoc")
    if found.get("assoc", {}).get("path") is None:
        put("assoc", first(["03_gwas/tables/gwas.assoc.linear",
                            "03_gwas/tables/gwas.assoc.logistic",
                            "03_gwas/tables/*.assoc.linear",
                            "03_gwas/tables/*.assoc.logistic"]),
            "glob 03_gwas/tables/*.assoc.*")
    put("adjusted", _stage_art(t, "gwas", "adjusted"), "state.gwas.artifacts.adjusted")
    if found.get("adjusted", {}).get("path") is None:
        cand = None
        if found.get("assoc", {}).get("path"):
            p = Path(str(found["assoc"]["path"]) + ".adjusted")
            cand = p if p.exists() else None
        put("adjusted", cand or first(["03_gwas/tables/*.adjusted"]),
            "assoc + .adjusted")
    put("results", _stage_art(t, "gwas", "results"), "state.gwas.artifacts.results")
    if found.get("results", {}).get("path") is None:
        put("results", first(["03_gwas/tables/gwas_results.csv"]),
            "03_gwas/tables/gwas_results.csv")
    put("eigenvec", _stage_art(t, "qc", "eigenvec"), "state.qc.artifacts.eigenvec")
    put("eigenval", _stage_art(t, "qc", "eigenval"), "state.qc.artifacts.eigenval")
    if found.get("eigenvec", {}).get("path") is None:
        put("eigenvec", first(["01_qc/tables/clean_pca.eigenvec", "01_qc/tables/pca.eigenvec",
                               "02_popgen/tables/pca.eigenvec"]), "glob *pca.eigenvec")
    if found.get("eigenval", {}).get("path") is None:
        put("eigenval", first(["01_qc/tables/clean_pca.eigenval", "01_qc/tables/pca.eigenval",
                               "02_popgen/tables/pca.eigenval"]), "glob *pca.eigenval")
    put("genome", _stage_art(t, "popgen", "genome"), "state.popgen.artifacts.genome")
    if found.get("genome", {}).get("path") is None:
        put("genome", first(["02_popgen/tables/genome.genome",
                             "01_qc/tables/genome.genome"]), "glob */genome.genome")
    put("rel", _stage_art(t, "popgen", "rel"), "state.popgen.artifacts.rel")
    if found.get("rel", {}).get("path") is None:
        put("rel", first(["02_popgen/tables/kinship.rel", "02_popgen/tables/*.rel"]),
            "glob 02_popgen/tables/kinship.rel")
    # GWAS 日志：优先含 lambda 行的那份（用于 λGC 对照）
    put("gwas_log", None, "n/a")
    for cand in sorted((rd / "03_gwas").glob("**/*.log")):
        try:
            if "lambda" in cand.read_text(encoding="utf-8", errors="replace").lower():
                put("gwas_log", cand, str(cand.relative_to(rd)))
                break
        except Exception:
            continue
    # 基因型产物统一规范为 PLINK 前缀（state 可能是前缀，glob 可能带 .bed）
    for key in ("clean", "base"):
        if found.get(key, {}).get("path") is not None:
            found[key]["path"] = _prefix_of(found[key]["path"])
    t.artifacts = {k: (str(v["path"]) if v["path"] else None) for k, v in found.items()}
    t.artifacts["_sources"] = {k: v["source"] for k, v in found.items()}
    return found


# ---------------------------------------------------------------------------
# PLINK 调用 / 文件读取
# ---------------------------------------------------------------------------
def _plink(t: _Target) -> Any:
    if t.plink is not None and hasattr(t.plink, "run"):
        return t.plink
    exe = _find_plink_exe(t)
    runtime = t.config.get("runtime") or {}
    if exe:
        try:
            t.plink = utils.Plink(
                exe, t.run_dir,
                threads=int(runtime.get("threads", 4) or 4),
                memory_mb=int(runtime.get("memory_mb", 4000) or 4000),
                logger=t.logger,
                timeout_s=int(runtime.get("timeout_s", 7200) or 7200),
            )
            return t.plink
        except Exception as exc:  # pragma: no cover
            t.logger.warning("构造 utils.Plink 失败，改用最小封装: %s", exc)
            t.plink = _MiniPlink(exe, t.run_dir, t.logger)
            return t.plink
    t.plink = None
    return None


def _run_plink(t: _Target, args: list[str], name: str, out_prefix: Path,
               expect: tuple[str, ...] = ()) -> Any:
    pk = _plink(t)
    if pk is None:
        raise FileNotFoundError("未找到 PLINK 可执行文件（设置 PLINK_EXE 或安装到 "
                                f"{DEFAULT_PLINK}）")
    return pk.run(args, stage="calibration", name=f"validate_{name}",
                  out_prefix=out_prefix, expect=list(expect), check=True)


def _read_tokens(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8-sig", errors="replace").split()


def _read_plink_table(path: Path) -> pd.DataFrame:
    """读取 PLINK 文本输出；自动识别空白 / 逗号分隔（stage 可能用 pandas 转存 CSV）。"""
    first = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()[0]
    sep = "," if ("," in first and "\t" not in first) else r"\s+"
    df = pd.read_csv(path, sep=sep, encoding="utf-8-sig")
    df.columns = [str(c).lstrip("#").strip() for c in df.columns]
    return df


@dataclass
class _Raw:
    """``--recode A`` 剂量矩阵（前两列 FID/IID，其余列 = SNP）。"""
    df: pd.DataFrame
    snps: list[str]
    alleles: dict[str, str]

    @property
    def samples(self) -> pd.DataFrame:
        return self.df[["FID", "IID"]]

    def matrix(self) -> np.ndarray:
        return self.df[self.snps].to_numpy(dtype=np.float64)


def _load_raw(raw_path: Path, bim: pd.DataFrame) -> _Raw:
    header = raw_path.read_text(encoding="utf-8-sig", errors="replace").split("\n", 1)[0].split()
    if len(header) < 7:
        raise ValueError(f".raw 文件列数异常: {raw_path}")
    dose_tokens = header[6:]
    # 用 .bim 的 SNP+A1 精确匹配列名（SNP ID 含下划线时也不歧义）
    exact = {}
    for i in range(len(bim)):
        exact.setdefault(f"{bim.iloc[i]['SNP']}_{bim.iloc[i]['A1']}", i)
    idx: list[int] = []
    used: set[int] = set()
    for tok in dose_tokens:
        i = exact.get(tok)
        if i is not None and i not in used:
            idx.append(i)
            used.add(i)
    if len(idx) != len(dose_tokens):
        # 兜底：PLINK recode 保持 bfile 位点顺序，按行对齐
        idx = list(range(len(dose_tokens)))
    sub = bim.iloc[idx[:min(len(idx), len(bim))]].reset_index(drop=True)
    df = pd.read_csv(raw_path, sep=r"\s+", na_values=["NA", "-9"], encoding="utf-8-sig")
    if df.shape[1] != len(header):
        raise ValueError(".raw 表头列数与数据行不一致")
    df.columns = header[:6] + [str(x) for x in sub["SNP"].tolist()]
    # 重复 SNP ID：保留第一个（QC 后不应出现，防御性处理）
    dup = df.columns.duplicated()
    if dup.any():
        df = df.loc[:, ~dup]
    snps = [str(x) for x in sub["SNP"].tolist()]
    seen: set[str] = set()
    snps = [s for s in snps if not (s in seen or seen.add(s))]
    alleles = {str(bim.iloc[i]["SNP"]): str(bim.iloc[i]["A1"]) for i in idx}
    return _Raw(df=df, snps=snps, alleles=alleles)


def _recode(t: _Target, clean: Path, snps: list[str], tag: str) -> _Raw:
    """对 clean bfile 的指定 SNP 子集执行 --recode A（结果按 tag 缓存）。"""
    clean = _prefix_of(clean)
    cache: dict[str, _Raw] = getattr(t, "_raw_cache", {})
    setattr(t, "_raw_cache", cache)
    if tag in cache:
        return cache[tag]
    list_path = t.data_dir / f"recode_{tag}.snplist"
    list_path.write_text("\n".join(map(str, snps)) + "\n", encoding="utf-8")
    prefix = t.data_dir / f"recode_{tag}"
    _run_plink(t, ["--bfile", str(clean), "--extract", str(list_path),
                   "--recode", "A", "--keep-allele-order"], f"recode_{tag}", prefix,
               expect=(".raw",))
    bim = utils.read_bim(Path(str(clean) + ".bim"))
    raw = _load_raw(Path(str(prefix) + ".raw"), bim)
    cache[tag] = raw
    return raw


def _read_eigenvec_robust(path: Path) -> pd.DataFrame:
    first = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()[0]
    tokens = first.split()
    if tokens and tokens[0].upper() in {"FID", "#FID"}:
        df = pd.read_csv(path, sep=r"\s+", header=0, encoding="utf-8-sig")
        df.columns = [str(c).lstrip("#").strip() for c in df.columns]
        return df
    df = pd.read_csv(path, sep=r"\s+", header=None, encoding="utf-8-sig")
    df.columns = ["FID", "IID"] + [f"PC{i + 1}" for i in range(df.shape[1] - 2)]
    return df


def _detect_header(path: Path) -> bool:
    first = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()[0].split()
    if len(first) < 2:
        return False
    return first[0].upper() in {"FID", "#FID"} or first[1].upper() in {"IID", "#IID"}


def _read_pheno(path: Path) -> pd.DataFrame:
    has_header = _detect_header(path)
    df = pd.read_csv(path, sep=r"\s+", header=0 if has_header else None,
                     encoding="utf-8-sig")
    if has_header:
        df.columns = [str(c).lstrip("#").strip() for c in df.columns]
        cols = list(df.columns)
        df = df[[cols[0], cols[1], cols[2]]].copy()
        df.columns = ["FID", "IID", "PHENO"]
    else:
        df = df.iloc[:, :3].copy()
        df.columns = ["FID", "IID", "PHENO"]
    df["FID"], df["IID"] = df["FID"].astype(str), df["IID"].astype(str)
    df["PHENO"] = pd.to_numeric(df["PHENO"], errors="coerce")
    df = df.loc[df["PHENO"].notna() & (df["PHENO"] != -9)]
    return df


# ---------------------------------------------------------------------------
# 单项检查
# ---------------------------------------------------------------------------
def _ck(name: str, category: str, status: str, detail: str,
        value: Any = None, threshold: Any = None) -> dict[str, Any]:
    return {"name": name, "category": category, "status": status,
            "detail": detail, "value": value, "threshold": threshold}


def _guard(records: list[dict[str, Any]], name: str, category: str, fn) -> None:
    try:
        rec = fn()
        if rec:
            records.append(rec)
    except Exception as exc:  # noqa: BLE001
        import traceback
        LOGGER.warning("校验 %s 无法执行: %s: %s", name, type(exc).__name__, exc)
        LOGGER.debug("%s", traceback.format_exc())
        records.append(_ck(name, category, "warn",
                           f"无法执行该项校验: {type(exc).__name__}: {exc}"))


def check_file_structure(t: _Target, found: dict) -> dict:
    clean = found["clean"]["path"]
    if clean is None:
        return _ck("structure_counts", "结构", "warn", "缺少 QC 后 bfile（01_qc/data/clean）")
    fam_p, bim_p, bed_p = (Path(str(clean) + s) for s in (".fam", ".bim", ".bed"))
    fam, bim = utils.read_fam(fam_p), utils.read_bim(bim_p)
    n_s, n_v = len(fam), len(bim)
    expected_bytes = 3 + n_v * math.ceil(n_s / 4)
    actual_bytes = bed_p.stat().st_size
    magic = bed_p.read_bytes()[:3]
    qc = t.metrics.get("qc") or {}
    gwas = t.metrics.get("gwas") or {}
    problems: list[str] = []
    if actual_bytes != expected_bytes:
        problems.append(f".bed 字节数 {actual_bytes} != 3 + {n_v}*ceil({n_s}/4) = {expected_bytes}")
    if magic != b"\x6c\x1b\x01":
        problems.append(f".bed 魔数异常 {magic!r}")
    for expect, actual, label, src in (
        (qc.get("n_samples_out"), n_s, "n_samples", "01_qc/metrics"),
        (gwas.get("n_samples"), n_s, "n_samples", "03_gwas/metrics"),
        (qc.get("n_variants_out"), n_v, "n_variants", "01_qc/metrics"),
        (gwas.get("n_tested"), n_v, "n_variants", "03_gwas/metrics"),
    ):
        if expect is not None and int(expect) != actual:
            problems.append(f"{src}.{label}={expect} 与文件实际 {actual} 不一致")
    # 与 PLINK 日志对照（clean.log 由 --make-bed 产生）
    for log in sorted(clean.parent.glob("*.log")):
        try:
            counts = (utils.parse_plink_log(
                log.read_text(encoding="utf-8", errors="replace")).get("counts") or {})
        except Exception:
            continue
        n_pass = counts.get("n_pass_qc")
        if isinstance(n_pass, (list, tuple)) and len(n_pass) == 2:
            v_log, s_log = int(n_pass[0]), int(n_pass[1])
            if (v_log, s_log) != (n_v, n_s):
                problems.append(f"{log.name} 记录 {s_log} 样本 x {v_log} 位点，与文件不符")
            break
    status = "pass" if not problems else "fail"
    detail = (f"{n_s} 样本 x {n_v} 位点（{len(bim['CHR'].unique())} 条染色体），"
              f".bed {actual_bytes} 字节，魔数 6c1b01" if not problems else "；".join(problems))
    return _ck("structure_counts", "结构", status, detail,
               value={"n_samples": n_s, "n_variants": n_v,
                      "n_chromosomes": int(len(bim["CHR"].unique())),
                      "bed_bytes": int(actual_bytes)},
               threshold=".bed 字节数 = 3 + n_v*ceil(n_s/4)，且与 metrics/PLINK 日志一致")


def check_duplicate_ids(t: _Target, found: dict) -> dict:
    clean = found["clean"]["path"]
    fam = utils.read_fam(Path(str(clean) + ".fam"))
    bim = utils.read_bim(Path(str(clean) + ".bim"))
    dup_samples = int(fam.duplicated(subset=["FID", "IID"]).sum())
    dup_snps = int(bim["SNP"].duplicated().sum())
    return _ck("duplicate_ids", "结构",
               "pass" if dup_samples == 0 and dup_snps == 0 else "fail",
               f"重复样本 ID {dup_samples} 个，重复 SNP ID {dup_snps} 个",
               value={"duplicate_samples": dup_samples, "duplicate_snps": dup_snps},
               threshold=0)


def check_alleles(t: _Target, found: dict) -> dict:
    bim = utils.read_bim(Path(str(found["clean"]["path"]) + ".bim"))
    bad = bim.loc[~bim["A1"].astype(str).isin(VALID_ALLELES) |
                  ~bim["A2"].astype(str).isin(VALID_ALLELES)]
    same = bim.loc[bim["A1"].astype(str) == bim["A2"].astype(str)]
    detail = f"{len(bim)} 个位点等位基因合法"
    if len(bad):
        detail += (f"；非法 {len(bad)} 个，如 "
                   f"{bad.head(3)[['SNP', 'A1', 'A2']].values.tolist()}")
    if len(same):
        detail += f"；A1==A2 的位点 {len(same)} 个"
    return _ck("allele_legality", "结构", "pass" if not len(bad) and not len(same) else "fail",
               detail, value={"invalid": int(len(bad)), "identical": int(len(same))},
               threshold="A1/A2 ∈ {A,C,G,T} 且 A1≠A2")


def check_qc_arithmetic(t: _Target, found: dict) -> dict:
    qc = t.metrics.get("qc") or {}
    clean = found["clean"]["path"]
    if not qc or clean is None:
        return _ck("qc_filter_arithmetic", "QC", "warn",
                   "缺少 01_qc/metrics.json 或 clean bfile")
    n_s = len(utils.read_fam(Path(str(clean) + ".fam")))
    n_v = len(utils.read_bim(Path(str(clean) + ".bim")))
    removed = qc.get("removed") or {}
    s_rm = int((removed.get("samples") or {}).get("total",
                                                   removed.get("samples_total", 0)) or 0)
    v_rm = int((removed.get("variants") or {}).get("total", 0) or 0)
    problems: list[str] = []
    if qc.get("n_samples_in") is not None and int(qc["n_samples_in"]) - s_rm != n_s:
        problems.append(f"样本: {qc.get('n_samples_in')} - {s_rm} != {n_s}")
    if qc.get("n_variants_in") is not None and int(qc["n_variants_in"]) - v_rm != n_v:
        problems.append(f"位点: {qc.get('n_variants_in')} - {v_rm} != {n_v}")
    if int(qc.get("n_samples_out", n_s)) != n_s:
        problems.append(f"metrics.n_samples_out={qc.get('n_samples_out')} != 文件 {n_s}")
    if int(qc.get("n_variants_out", n_v)) != n_v:
        problems.append(f"metrics.n_variants_out={qc.get('n_variants_out')} != 文件 {n_v}")
    intake = t.metrics.get("intake") or {}
    if intake.get("n_samples") is not None and qc.get("n_samples_in") is not None:
        if int(intake["n_samples"]) != int(qc["n_samples_in"]):
            problems.append(f"intake.n_samples={intake['n_samples']} != "
                            f"qc.n_samples_in={qc['n_samples_in']}")
    if intake.get("n_variants") is not None and qc.get("n_variants_in") is not None:
        if int(intake["n_variants"]) != int(qc["n_variants_in"]):
            problems.append(f"intake.n_variants={intake['n_variants']} != "
                            f"qc.n_variants_in={qc['n_variants_in']}")
    return _ck("qc_filter_arithmetic", "QC", "pass" if not problems else "fail",
               "样本/位点移除计数与 QC 后文件一致" if not problems else "；".join(problems),
               value={"n_samples_out": n_s, "n_variants_out": n_v,
                      "samples_removed": s_rm, "variants_removed": v_rm},
               threshold="n_in - removed == 文件实际计数")


def check_prune_counts(t: _Target, found: dict) -> dict:
    prune = found["pruned"]["path"]
    qc = t.metrics.get("qc") or {}
    if prune is None:
        return _ck("ld_prune_counts", "QC", "warn", "缺少 LD 剪枝列表 *.prune.in")
    n_in = len(_read_tokens(prune))
    out_path = Path(str(prune).replace(".prune.in", ".prune.out"))
    n_out = len(_read_tokens(out_path)) if out_path.exists() else None
    n_total = int(qc.get("n_variants_out") or 0)
    expected = int(((qc.get("pruned") or {}).get("n_variants")) or n_in)
    problems = []
    if n_in != expected:
        problems.append(f"prune.in 行数 {n_in} != metrics.pruned.n_variants {expected}")
    if n_out is not None and n_total and n_in + n_out != n_total:
        problems.append(f"prune.in {n_in} + prune.out {n_out} != n_variants_out {n_total}")
    return _ck("ld_prune_counts", "QC", "pass" if not problems else "fail",
               f"LD 剪枝保留 {n_in} 个位点"
               + (f"，剔除 {n_out} 个" if n_out is not None else "")
               + ("" if not problems else "；" + "；".join(problems)),
               value={"n_pruned": n_in, "n_removed": n_out},
               threshold="prune.in 计数 == metrics.pruned.n_variants")


def check_maf(t: _Target, found: dict) -> dict:
    """MAF / 等位基因计数独立重算。

    ``--frq`` 的 MAF 只打印 4 位小数（比较误差可达 5e-5），所以用 ``--freqx``
    的**整数等位基因计数**做主对照：A1 计数 / 缺失数逐位点精确相等，再由计数
    重算 MAF，可达到 DESIGN 要求的 <1e-5。
    """
    clean = _prefix_of(found["clean"]["path"])
    bim = utils.read_bim(Path(str(clean) + ".bim"))
    fam = utils.read_fam(Path(str(clean) + ".fam"))
    prefix = t.data_dir / "clean_freqx"
    # --keep-allele-order：否则 PLINK 会把 A1 重排为次等位基因，计数无法与 .bim 对齐
    _run_plink(t, ["--bfile", str(clean), "--freqx", "--keep-allele-order"],
               "clean_freqx", prefix, expect=(".frqx",))
    freqx = pd.read_csv(Path(str(prefix) + ".frqx"), sep="\t", encoding="utf-8-sig")
    freqx.columns = [str(c).strip() for c in freqx.columns]
    rng = np.random.default_rng(int(t.tol["seed"]))
    n_pick = min(int(t.tol["n_maf_snps"]), len(bim))
    pick = bim.iloc[np.sort(rng.choice(len(bim), n_pick, replace=False))]
    raw = _recode(t, clean, pick["SNP"].astype(str).tolist(), "maf")
    # 统计口径：PLINK --freqx 默认只统计 founders（PAT/MAT=0）
    samples = raw.df[["FID", "IID"]].copy()
    samples["FID"], samples["IID"] = samples["FID"].astype(str), samples["IID"].astype(str)
    key_founder = {(str(f), str(i)): (str(pa) == "0" and str(ma) == "0")
                   for f, i, pa, ma in zip(fam["FID"], fam["IID"], fam["PAT"], fam["MAT"])}
    all_founders = all(key_founder.values())
    x = raw.matrix()
    if not all_founders:
        mask = np.array([key_founder.get((f, i), False)
                         for f, i in zip(samples["FID"], samples["IID"])])
        x = x[mask, :]
    n_ok = np.sum(~np.isnan(x), axis=0).astype(np.int64)
    c1 = np.nansum(x, axis=0).astype(np.int64)
    n_miss = np.sum(np.isnan(x), axis=0).astype(np.int64)
    mine = pd.DataFrame({"SNP": raw.snps, "c_a1": c1, "n_obs": n_ok, "n_miss": n_miss})
    col_c1, col_het, col_a2, col_miss = "C(HOM A1)", "C(HET)", "C(HOM A2)", "C(MISSING)"
    for c in (col_c1, col_het, col_a2, col_miss):
        if c not in freqx.columns:
            raise ValueError(f".frqx 缺少列 {c}（实际列: {freqx.columns.tolist()}）")
    freqx["c_a1"] = 2 * freqx[col_c1].astype(np.int64) + freqx[col_het].astype(np.int64)
    freqx["c_a2"] = 2 * freqx[col_a2].astype(np.int64) + freqx[col_het].astype(np.int64)
    freqx["n_miss"] = freqx[col_miss].astype(np.int64)
    fq = freqx[["SNP", "A1", "c_a1", "c_a2", "n_miss"]].rename(columns={
        "c_a1": "c_a1_plink", "c_a2": "c_a2_plink", "n_miss": "n_miss_plink"})
    py = mine[["SNP", "c_a1", "n_obs", "n_miss"]].rename(columns={
        "c_a1": "c_a1_py", "n_miss": "n_miss_py"})
    merged = fq.merge(py, on="SNP", how="inner")
    if merged.empty:
        return _ck("maf_recompute", "频率", "warn", "无法对齐 .frqx 与 .raw 的 SNP 列表")
    bad_a1 = int(np.sum(merged["c_a1_plink"].to_numpy() != merged["c_a1_py"].to_numpy()))
    bad_miss = int(np.sum(merged["n_miss_plink"].to_numpy() != merged["n_miss_py"].to_numpy()))
    tot_plink = (merged["c_a1_plink"] + merged["c_a2_plink"]).to_numpy(dtype=float)
    tot_py = 2.0 * merged["n_obs"].to_numpy(dtype=float)
    maf_plink = np.minimum(merged["c_a1_plink"], merged["c_a2_plink"]).to_numpy(dtype=float) / np.maximum(tot_plink, 1)
    maf_py = np.minimum(merged["c_a1_py"].to_numpy(dtype=float),
                        tot_py - merged["c_a1_py"].to_numpy(dtype=float)) / np.maximum(tot_py, 1)
    max_maf_diff = float(np.nanmax(np.abs(maf_plink - maf_py)))
    n_allele_mismatch = int(np.sum(merged["A1"].astype(str).to_numpy() !=
                                   np.array([raw.alleles.get(s, "") for s in merged["SNP"]])))
    # 附带：--frq 打印值（4 位小数）一致性，仅提示
    frq_path = Path(str(prefix) + ".frq")
    frq_disp = None
    if frq_path.exists():
        try:
            frq = _read_plink_table(frq_path)[["SNP", "MAF"]]
            frq["MAF"] = pd.to_numeric(frq["MAF"], errors="coerce")
            m2 = merged[["SNP"]].assign(maf_py=maf_py).merge(frq, on="SNP", how="inner")
            frq_disp = float(np.nanmax(np.abs(m2["MAF"] - m2["maf_py"])))
        except Exception:
            frq_disp = None
    problems: list[str] = []
    if bad_a1 or bad_miss or n_allele_mismatch:
        problems.append(f"等位基因计数/缺失数与 PLINK 不一致：A1 计数 {bad_a1} 个位点、"
                        f"缺失计数 {bad_miss} 个位点、等位基因符号 {n_allele_mismatch} 个")
    if max_maf_diff > float(t.tol["maf"]):
        problems.append(f"重算 MAF 最大绝对差 {max_maf_diff:.2e} > {t.tol['maf']:.0e}")
    detail = (f"随机 {len(merged)} 个位点：A1 等位基因计数/缺失数逐位点与 --freqx 精确比对，"
              f"不一致 {bad_a1 + bad_miss} 项；重算 MAF 最大绝对差 {max_maf_diff:.2e}"
              f"（容差 {t.tol['maf']:.0e}）")
    if frq_disp is not None:
        detail += f"；--frq 打印值差异 {frq_disp:.2e}（4 位小数舍入）"
    if problems:
        detail += "；" + "；".join(problems)
    return _ck("maf_recompute", "频率", "pass" if not problems else "fail", detail,
               value={"max_abs_maf_diff": max_maf_diff,
                      "allele_count_mismatch": bad_a1 + bad_miss,
                      "n_snps": int(len(merged)),
                      "n_allele_symbol_mismatch": n_allele_mismatch,
                      "frq_display_diff": frq_disp},
               threshold=f"等位基因计数精确一致；max|ΔMAF| <= {t.tol['maf']}")


def check_missingness(t: _Target, found: dict) -> dict:
    """缺失率独立重算：同一 SNP 子集上 ``--missing`` 的整数 N_MISS vs .raw 的 NaN 计数。"""
    clean = _prefix_of(found["clean"]["path"])
    bim = utils.read_bim(Path(str(clean) + ".bim"))
    rng = np.random.default_rng(int(t.tol["seed"]) + 1)
    n_pick = min(int(t.tol["n_missing_snps"]), len(bim))
    pick = bim.iloc[np.sort(rng.choice(len(bim), n_pick, replace=False))]
    raw = _recode(t, clean, pick["SNP"].astype(str).tolist(), "missing")
    list_path = t.data_dir / "recode_missing.snplist"
    prefix = t.data_dir / "clean_missing"
    _run_plink(t, ["--bfile", str(clean), "--extract", str(list_path), "--missing"],
               "clean_missing", prefix, expect=(".lmiss", ".imiss"))
    lmiss = _read_plink_table(Path(str(prefix) + ".lmiss"))
    imiss = _read_plink_table(Path(str(prefix) + ".imiss"))
    x = raw.matrix()
    var_miss = pd.DataFrame({"SNP": raw.snps,
                             "py_nmiss": np.sum(np.isnan(x), axis=0).astype(np.int64)})
    merged = lmiss[["SNP", "N_MISS"]].merge(var_miss, on="SNP", how="inner")
    merged["N_MISS"] = pd.to_numeric(merged["N_MISS"], errors="coerce").astype(np.int64)
    bad_var = int(np.sum(merged["N_MISS"].to_numpy() != merged["py_nmiss"].to_numpy()))
    d1 = 0.0 if not bad_var else float(np.nanmax(
        np.abs(merged["N_MISS"].to_numpy() - merged["py_nmiss"].to_numpy())
        / max(len(raw.df), 1)))
    sample_miss = pd.DataFrame({"IID": raw.df["IID"].astype(str),
                                "py_nmiss": np.sum(np.isnan(x), axis=1).astype(np.int64)})
    sm = imiss[["IID", "N_MISS"]].merge(sample_miss, on="IID", how="inner")
    sm["N_MISS"] = pd.to_numeric(sm["N_MISS"], errors="coerce").astype(np.int64)
    bad_sample = int(np.sum(sm["N_MISS"].to_numpy() != sm["py_nmiss"].to_numpy()))
    d2 = 0.0 if not bad_sample else float(np.nanmax(
        np.abs(sm["N_MISS"].to_numpy() - sm["py_nmiss"].to_numpy())
        / max(len(raw.snps), 1)))
    tol = float(t.tol["missing"])
    max_diff = float(max(d1, d2))
    problems: list[str] = []
    if bad_var or bad_sample:
        problems.append(f"N_MISS 不一致：位点 {bad_var}/{len(merged)} 个、"
                        f"样本 {bad_sample}/{len(sm)} 个")
    if max_diff > tol:
        problems.append(f"缺失率最大差 {max_diff:.2e} > {tol:.0e}")
    mean_rate = float(np.mean(np.isnan(x)))
    qc = t.metrics.get("qc") or {}
    mr = qc.get("missing_rate") or {}
    detail = (f"{len(merged)} 个位点、{len(sm)} 个样本的缺失数分别与 --missing 的 N_MISS "
              f"精确比对（不一致 {bad_var + bad_sample} 项）；子集总体缺失率 {mean_rate:.4f}")
    if mr.get("variant_mean") is not None and abs(mean_rate - float(mr["variant_mean"])) > 5e-3:
        detail += f"；注意：与 qc.missing_rate.variant_mean={mr['variant_mean']} 差异较大" \
                  "（子集抽样或质控顺序不同）"
    if problems:
        detail += "；" + "；".join(problems)
    return _ck("missingness_recompute", "频率", "pass" if not problems else "fail", detail,
               value={"max_fraction_diff": max_diff,
                      "variant_mismatch": bad_var, "sample_mismatch": bad_sample,
                      "n_snps": int(len(merged)), "n_samples": int(len(sm)),
                      "subset_missing_rate": mean_rate},
               threshold=f"整数 N_MISS 精确一致；max|ΔF_MISS| <= {tol:.0e}")


def check_hwe(t: _Target, found: dict) -> dict:
    """QC 后 HWE 极值比例。

    二元表型按 PLINK ``--hwe`` 标准语义看**对照（UNAFF）**内的 HWE；数量性状看 ALL。
    """
    clean = _prefix_of(found["clean"]["path"])
    pheno_type = ((t.config.get("_derived") or {}).get("phenotype_type")
                  or (t.metrics.get("gwas") or {}).get("phenotype_type")
                  or ((t.metrics.get("intake") or {}).get("phenotype") or {}).get("type")
                  or "quantitative")
    pheno_p = found["pheno"]["path"]
    prefix = t.data_dir / "clean_hwe"
    args = ["--bfile", str(clean)]
    if pheno_type == "binary" and pheno_p is not None:
        args += ["--pheno", str(pheno_p)]
    args += ["--hardy", "midp"]
    _run_plink(t, args, "clean_hwe", prefix, expect=(".hwe",))
    hwe = _read_plink_table(Path(str(prefix) + ".hwe"))
    if "TEST" in hwe.columns:
        test_col = hwe["TEST"].astype(str)
        if pheno_type == "binary" and test_col.str.startswith("UNAFF").any():
            hwe = hwe.loc[test_col.str.startswith("UNAFF")].copy()
            scope = "对照组(UNAFF)"
        else:
            hwe = hwe.loc[test_col.str.startswith("ALL")].copy()
            scope = "全部样本(ALL)"
    else:
        scope = "全部样本"
    p = pd.to_numeric(hwe.get("P"), errors="coerce").dropna().to_numpy(dtype=float)
    frac = float(np.mean(p < 1e-6)) if p.size else float("nan")
    warn_frac, fail_frac = float(t.tol["hwe_warn"]), float(t.tol["hwe_fail"])
    status = "pass" if frac <= warn_frac else ("warn" if frac <= fail_frac else "fail")
    return _ck("hwe_extremes", "QC", status,
               f"QC 后 {p.size} 个位点（{scope}）HWE(midp) p<1e-6 比例 {frac:.4%}"
               f"（warn>{warn_frac:.1%}，fail>{fail_frac:.1%}）",
               value=frac, threshold={"warn": warn_frac, "fail": fail_frac})


def check_phenotype(t: _Target, found: dict) -> dict:
    pheno_path, clean = found["pheno"]["path"], found["clean"]["path"]
    if pheno_path is None or clean is None:
        return _ck("phenotype_coding", "表型", "warn", "缺少表型文件或 clean bfile")
    df = _read_pheno(pheno_path)
    pheno_type = ((t.config.get("_derived") or {}).get("phenotype_type")
                  or (t.metrics.get("gwas") or {}).get("phenotype_type")
                  or ((t.metrics.get("intake") or {}).get("phenotype") or {}).get("type")
                  or "quantitative")
    fam = utils.read_fam(Path(str(clean) + ".fam"))
    fam_keys = set(zip(fam["FID"].astype(str), fam["IID"].astype(str)))
    in_fam = np.array([(f, i) in fam_keys for f, i in zip(df["FID"], df["IID"])])
    present = df.loc[in_fam]
    problems: list[str] = []
    value: dict[str, Any] = {"type": pheno_type, "n_rows": int(len(df)),
                             "n_in_fam": int(in_fam.sum()), "n_usable": int(len(present))}
    if int(in_fam.sum()) < len(fam) - 5:
        problems.append(f"表型仅覆盖 QC 后样本 {int(in_fam.sum())}/{len(fam)}")
    elif int(in_fam.sum()) < len(fam):
        problems.append(f"{len(fam) - int(in_fam.sum())} 个 QC 后样本无表型记录")
    if pheno_type == "binary":
        values = set(np.unique(present["PHENO"].to_numpy()))
        value["prevalence"] = float(np.mean(present["PHENO"] == 2)) if len(present) else None
        if values <= {1, 2}:
            pass
        elif values <= {0, 1}:
            problems.append("二元表型编码为 0/1；PLINK 1.9 --logistic 默认 1=对照、2=病例，"
                            "需转换为 1/2")
        else:
            problems.append(f"二元表型取值异常 {sorted(values)[:6]}")
        prev = value.get("prevalence")
        if prev is not None and not (0.02 <= prev <= 0.6):
            problems.append(f"病例比例 {prev:.3f} 不合理（应在 0.02~0.6）")
    else:
        n_unique = int(present["PHENO"].nunique())
        value["n_unique"] = n_unique
        if n_unique < 20:
            problems.append(f"数量性状只有 {n_unique} 个不同取值，疑似编码错误")
        if not np.isfinite(present["PHENO"].to_numpy()).all():
            problems.append("表型含 NaN/Inf")
    detail = (f"表型类型={pheno_type}，{len(df)} 行，匹配 QC 后样本 {int(in_fam.sum())}，"
              f"可用 {len(present)}")
    if pheno_type == "binary":
        detail += f"，病例比例 {value.get('prevalence')}"
    if problems:
        detail += "；" + "；".join(problems)
    return _ck("phenotype_coding", "表型", "pass" if not problems else "fail", detail,
               value=value,
               threshold="二元 1/2 编码且患病率合理；数量性状取值多样、无缺失异常")


def check_pca(t: _Target, found: dict) -> dict:
    vec_p, val_p = found["eigenvec"]["path"], found["eigenval"]["path"]
    clean, prune_p = _prefix_of(found["clean"]["path"]), found["pruned"]["path"]
    if vec_p is None or val_p is None:
        return _ck("pca_recompute", "群体结构", "warn", "缺少 PCA eigenvec/eigenval 产物")
    vec = _read_eigenvec_robust(vec_p)
    eigval = utils.read_eigenval(val_p)
    n_samples = len(vec)
    pc_cols = [str(c) for c in vec.columns if str(c).startswith("PC")]
    k = min(5, len(pc_cols), len(eigval))
    problems: list[str] = []
    value: dict[str, Any] = {"n_samples": n_samples, "n_pcs": int(len(eigval)),
                             "k_checked": int(k),
                             "eigenvector_source": _rel(t, vec_p),
                             "eigenvalue_source": _rel(t, val_p)}
    # (a) 方差解释比例：PLINK 1.9 特征值之和 = 样本数 N（关系矩阵的迹），不是 top-k 之和
    expected_var = eigval[:k] / float(n_samples)
    value["eigenvalues"] = [round(float(v), 6) for v in eigval[:k]]
    value["var_explained_expected"] = [round(float(v), 6) for v in expected_var]
    metrics_var = None
    for src in ("popgen", "qc"):
        pca = (t.metrics.get(src) or {}).get("pca") or {}
        if pca.get("var_explained"):
            metrics_var = np.asarray(pca["var_explained"][:k], dtype=float)
            value["var_explained_source"] = src
            break
    if metrics_var is not None:
        value["var_explained_metrics"] = [round(float(v), 6) for v in metrics_var]
        diff = float(np.max(np.abs(metrics_var - expected_var)))
        value["var_explained_max_abs_diff"] = diff
        if diff > 1e-3:
            problems.append(
                f"var_explained 与 eigval/N 不符（最大差 {diff:.3f}）：metrics 前 {k} PC="
                f"{np.round(metrics_var, 4).tolist()}，正确应为 "
                f"{np.round(expected_var, 4).tolist()}（分母是样本数 N={n_samples}，"
                "不是 top-k 特征值之和）")
    # (b) 特征向量范数约定（PLINK 输出单位范数 => var*N ≈ 1）
    if pc_cols:
        varN = np.var(vec[pc_cols].to_numpy(dtype=float), axis=0) * n_samples
        norm_dev = float(np.max(np.abs(varN[:k] - 1.0)))
        value["eigenvec_norm_deviation"] = norm_dev
        if norm_dev > 0.05:
            problems.append(f"eigenvec 列范数异常（var*N 最大偏差 {norm_dev:.3f}）")
    # (c) 独立重算 PCA（pruned 位点上再做一次特征分解）
    n_pruned = len(_read_tokens(prune_p)) if prune_p is not None else 0
    max_values = int(t.tol["pca_max_values"])
    max_flops = float(t.tol["pca_max_flops"])
    if clean is not None and prune_p is not None and n_pruned > 10:
        pruned = _read_tokens(prune_p)
        rng = np.random.default_rng(int(t.tol["seed"]) + 2)
        max_by_mem = max(200, max_values // max(n_samples, 1))
        max_by_flops = max(200, int(max_flops / max(n_samples ** 2, 1)))
        m_use = min(len(pruned), max_by_mem, max_by_flops)
        subsampled = m_use < len(pruned)
        if subsampled:
            sel = np.sort(rng.choice(len(pruned), m_use, replace=False))
            used = [pruned[i] for i in sel]
        else:
            used = pruned
        try:
            raw = _recode(t, clean, used, "pca")
            m = len(raw.snps)
            x = raw.matrix()
            p = np.clip(np.nanmean(x, axis=0) / 2.0, 1e-6, 1 - 1e-6)
            x = np.where(np.isnan(x), 2.0 * p[None, :], x)
            z = ((x - 2.0 * p[None, :]) /
                 np.sqrt(2.0 * p * (1.0 - p))[None, :]).astype(np.float64)
            gram = (z @ z.T) / float(m)          # ≈ PLINK 的方差标准化关系矩阵
            w, V = np.linalg.eigh(gram)
            order = np.argsort(w)[::-1]
            w, V = w[order], V[:, order]
            # gram = (1/m) Z Z^T 与 PLINK 的方差标准化关系矩阵同尺度，
            # 其特征值即与 .eigenval 同尺度（不要再除以 m）
            my_eig = w
            corr: dict[str, float] = {}
            for i in range(k):
                u = vec[pc_cols[i]].to_numpy(dtype=float)
                u = u - u.mean()
                u = u / max(np.linalg.norm(u), 1e-12)
                v = V[:, i]
                v = v - v.mean()
                v = v / max(np.linalg.norm(v), 1e-12)
                corr[f"PC{i + 1}"] = round(abs(float(np.dot(u, v))), 6)
            value["pca_direction_corr"] = corr
            value["pca_n_snps_used"] = m
            value["pca_subsampled"] = bool(subsampled)
            value["pca_eigen_ratio"] = [round(float(my_eig[i] / max(eigval[i], 1e-12)), 4)
                                        for i in range(min(k, len(my_eig)))]
            corr_min = float(t.tol["pca_corr_subsampled"] if subsampled
                             else t.tol["pca_corr"])
            if corr.get("PC1", 0) < corr_min:
                problems.append(f"PC1 方向与独立重算相关 {corr.get('PC1'):.3f} < {corr_min}"
                                f"（{m} 个位点{'(抽样)' if subsampled else ''}）")
            ratio1 = value["pca_eigen_ratio"][0] if value["pca_eigen_ratio"] else float("nan")
            if np.isfinite(ratio1) and abs(ratio1 - 1) > float(t.tol["eigen_scale"]):
                problems.append(f"PC1 特征值尺度比 {ratio1:.3f}（容差 ±{t.tol['eigen_scale']}）")
        except Exception as exc:  # noqa: BLE001
            value["recompute_error"] = f"{type(exc).__name__}: {exc}"
            problems.append(f"PCA 独立重算未完成: {exc}")
    else:
        value["recompute_skipped"] = "缺少 clean bfile 或 prune.in"

    status = "pass" if not problems else "fail"
    detail = (f"{n_samples} 样本 x {int(len(eigval))} 个 PC；PC1 特征值 {eigval[0]:.4f}，"
              f"方差解释 {expected_var[0]:.3%}" if len(eigval) else "无特征值")
    if value.get("pca_direction_corr"):
        detail += (f"；独立重算方向相关 PC1="
                   f"{value['pca_direction_corr'].get('PC1')}"
                   f"（{value.get('pca_n_snps_used')} 位点）")
    if problems:
        detail += "；" + "；".join(problems)
    return _ck("pca_recompute", "群体结构", status, detail, value=value,
               threshold={"var_explained": "eigval / N_samples",
                          "pca_corr": t.tol["pca_corr"],
                          "eigen_scale": t.tol["eigen_scale"]})


def check_relatedness(t: _Target, found: dict) -> dict:
    genome_p = found["genome"]["path"]
    if genome_p is None:
        return _ck("relatedness_range", "亲缘", "warn", "缺少 --genome 产物")
    g = _read_plink_table(genome_p)
    if "PI_HAT" not in g.columns:
        return _ck("relatedness_range", "亲缘", "warn", f"{genome_p.name} 缺少 PI_HAT 列")
    pi = pd.to_numeric(g["PI_HAT"], errors="coerce").to_numpy(dtype=float)
    pi = pi[np.isfinite(pi)]
    if pi.size == 0:
        return _ck("relatedness_range", "亲缘", "warn", "PI_HAT 全为空")
    problems: list[str] = []
    lo, hi = float(pi.min()), float(pi.max())
    if lo < -0.05 or hi > 1.05:
        problems.append(f"PI_HAT 越界 [{lo:.3f}, {hi:.3f}]（合理范围 [-0.05, 1.05]）")
    n_gt = int(np.sum(pi > 0.25))
    n_dup = int(np.sum(pi > 0.90))
    value = {"n_pairs": int(pi.size), "min": lo, "max": hi, "mean": float(pi.mean()),
             "n_pairs_gt_025": n_gt, "n_pairs_gt_090": n_dup}
    pop_rel = (t.metrics.get("popgen") or {}).get("relatedness") or {}
    if pop_rel.get("n_pairs_gt_025") is not None and int(pop_rel["n_pairs_gt_025"]) != n_gt:
        problems.append(f"metrics.n_pairs_gt_025={pop_rel['n_pairs_gt_025']} != 重算 {n_gt}")
    if pop_rel.get("max_pihat") is not None and abs(float(pop_rel["max_pihat"]) - hi) > 1e-6:
        problems.append(f"metrics.max_pihat={pop_rel['max_pihat']} != 重算 {hi:.4f}")
    pair_keys: set[tuple] = set()
    dup_pairs = 0
    for r in g[["FID1", "IID1", "FID2", "IID2"]].itertuples(index=False):
        key = tuple(sorted([(str(r[0]), str(r[1])), (str(r[2]), str(r[3]))]))
        if key in pair_keys:
            dup_pairs += 1
        pair_keys.add(key)
    if dup_pairs:
        problems.append(f"--genome 输出有 {dup_pairs} 个重复样本对")
    return _ck("relatedness_range", "亲缘", "pass" if not problems else "fail",
               f"{pi.size} 个样本对：PI_HAT ∈ [{lo:.3f}, {hi:.3f}]，"
               f">0.25 有 {n_gt} 对，>0.90 有 {n_dup} 对（近重复）"
               + ("" if not problems else "；" + "；".join(problems)),
               value=value, threshold="PI_HAT ∈ [-0.05, 1.05]，计数与 metrics 一致")


def check_kinship(t: _Target, found: dict) -> dict:
    rel_p = found["rel"]["path"]
    if rel_p is None:
        return _ck("kinship_matrix", "亲缘", "warn", "缺少 --make-rel 关系矩阵")
    mat = np.loadtxt(rel_p, dtype=float)
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        return _ck("kinship_matrix", "亲缘", "warn", f"关系矩阵维度异常 {mat.shape}")
    asym = float(np.max(np.abs(mat - mat.T)))
    diag = np.diag(mat)
    diag_mean_dev = float(abs(np.mean(diag) - 1.0))
    diag_max_dev = float(np.max(np.abs(diag - 1.0)))
    off = mat[~np.eye(mat.shape[0], dtype=bool)]
    problems: list[str] = []
    if asym > 1e-6:
        problems.append(f"矩阵不对称，最大 |A-A^T|={asym:.2e}")
    # --make-rel 的对角 = 1 + 近交系数，个体间可有数个百分点波动；看均值与极端值
    if diag_mean_dev > 0.05 or diag_max_dev > 0.5:
        problems.append(f"对角偏离 1（均值 {np.mean(diag):.3f}，最大偏差 {diag_max_dev:.3f}）")
    if off.size and (off.min() < -2.0 or off.max() > 2.0):
        problems.append(f"非对角元素越界 [{off.min():.3f}, {off.max():.3f}]")
    return _ck("kinship_matrix", "亲缘", "pass" if not problems else "fail",
               f"{mat.shape[0]}x{mat.shape[1]} 关系矩阵：对称性 {asym:.1e}，"
               f"对角均值 {diag.mean():.4f}，非对角范围 [{off.min():.3f}, {off.max():.3f}]"
               + ("" if not problems else "；" + "；".join(problems)),
               value={"n": int(mat.shape[0]), "asymmetry": asym,
                      "diag_mean": float(diag.mean()), "diag_max_dev": diag_max_dev,
                      "offdiag_min": float(off.min()), "offdiag_max": float(off.max())},
               threshold="对称（<1e-6）、对角≈1、非对角 ∈ [-2, 2]")


def _assoc_add(t: _Target, found: dict) -> pd.DataFrame:
    assoc_p = found["assoc"]["path"]
    if assoc_p is None:
        raise FileNotFoundError("缺少 GWAS .assoc 关联结果文件")
    df = _read_plink_table(assoc_p)
    if "TEST" in df.columns:
        df = df.loc[df["TEST"].astype(str).str.upper() == "ADD"].copy()
    for col in ("BETA", "OR", "STAT", "P", "NMISS", "SE", "L95", "U95"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def check_assoc_table(t: _Target, found: dict) -> dict:
    df = _assoc_add(t, found)
    gwas = t.metrics.get("gwas") or {}
    problems: list[str] = []
    n_tested = gwas.get("n_tested")
    if n_tested is not None and int(n_tested) != len(df):
        problems.append(f"metrics.n_tested={n_tested} != .assoc ADD 行数 {len(df)}")
    n_s = gwas.get("n_samples")
    if n_s and "NMISS" in df.columns:
        nmiss = df["NMISS"].dropna()
        if len(nmiss) and (nmiss.max() > int(n_s) or nmiss.min() < int(n_s) * 0.5):
            problems.append(f"NMISS 范围异常 [{nmiss.min():.0f}, {nmiss.max():.0f}]（n={n_s}）")
    if "P" in df.columns:
        # 口径：metrics.n_significant 对应 Bonferroni 阈值（DESIGN 第 8 节）；
        # 若没有 thresholds 则退回 5e-8
        bonf = (gwas.get("thresholds") or {}).get("bonferroni")
        alpha = float(bonf) if bonf else ALPHA_GW
        n_sig = int(np.sum(df["P"] < alpha))
        if gwas.get("n_significant") is not None and int(gwas["n_significant"]) != n_sig:
            problems.append(f"metrics.n_significant={gwas['n_significant']} != "
                            f"重算 {n_sig}（p<{alpha:.3g}）")
        value_alpha = alpha
    else:
        n_sig, value_alpha = None, None
    return _ck("gwas_assoc_table", "GWAS", "pass" if not problems else "fail",
               f"{len(df)} 个位点的加性检验结果"
               + ("" if not problems else "；" + "；".join(problems)),
               value={"n_rows": int(len(df)), "n_significant": n_sig,
                      "alpha_used": value_alpha},
               threshold="行数 == metrics.n_tested，NMISS/显著性计数一致")


def _load_covar(t: _Target, found: dict) -> tuple[pd.DataFrame, list[str]]:
    covar_p = found["covar"]["path"]
    if covar_p is None:
        raise FileNotFoundError("缺少 03_gwas/data/covar.txt")
    cols = list((t.metrics.get("gwas") or {}).get("covariates") or [])
    df = pd.read_csv(covar_p, sep=r"\s+", header=None, encoding="utf-8-sig",
                     dtype={0: str, 1: str})
    if cols and df.shape[1] == len(cols) + 2:
        df.columns = ["FID", "IID"] + cols
    else:
        df.columns = ["FID", "IID"] + [f"COV{i + 1}" for i in range(df.shape[1] - 2)]
        cols = list(df.columns[2:])
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["FID"], df["IID"] = df["FID"].astype(str), df["IID"].astype(str)
    return df, cols


def check_gwas_beta_se(t: _Target, found: dict) -> dict:
    """随机 N 个 SNP：statsmodels 在同协变量下重算 beta/SE 与 PLINK 结果对照。"""
    assoc = _assoc_add(t, found)
    if assoc.empty:
        return _ck("gwas_beta_se_recompute", "GWAS", "warn", "GWAS 关联结果为空")
    clean, pheno_p = _prefix_of(found["clean"]["path"]), found["pheno"]["path"]
    if clean is None or pheno_p is None:
        return _ck("gwas_beta_se_recompute", "GWAS", "warn", "缺少 clean bfile 或表型文件")
    pheno = _read_pheno(pheno_p)
    covar, cov_cols = _load_covar(t, found)
    gwas = t.metrics.get("gwas") or {}
    model = str(gwas.get("model") or
                ("logistic" if gwas.get("phenotype_type") == "binary" else "linear"))
    bim = utils.read_bim(Path(str(clean) + ".bim"))
    bim_keys = set(bim["SNP"].astype(str))
    candidates = [s for s in assoc["SNP"].astype(str) if s in bim_keys]
    if not candidates:
        return _ck("gwas_beta_se_recompute", "GWAS", "warn", ".assoc 与 .bim 的 SNP 无法对齐")
    rng = np.random.default_rng(int(t.tol["seed"]) + 1)
    n_pick = min(int(t.tol["n_beta_snps"]), len(candidates))
    picks = [str(x) for x in rng.choice(np.array(candidates), n_pick, replace=False)]
    raw = _recode(t, clean, picks, "beta")
    # 样本合并：raw 样本 ∩ pheno ∩ covar（与 PLINK --pheno/--covar 完整样本集合一致）
    base = raw.df[["FID", "IID"]].copy()
    base["FID"], base["IID"] = base["FID"].astype(str), base["IID"].astype(str)
    dat = (base.merge(pheno, on=["FID", "IID"], how="inner")
               .merge(covar, on=["FID", "IID"], how="inner"))
    for c in cov_cols:
        dat = dat.loc[dat[c].notna()]
    if len(dat) < 20:
        return _ck("gwas_beta_se_recompute", "GWAS", "warn",
                   f"可用样本只有 {len(dat)} 个，无法重算")
    y = dat["PHENO"].to_numpy(dtype=float)
    if model == "logistic":
        # PLINK 的二元表型是 1=对照/2=病例；statsmodels.Logit 要求 0/1
        vals = set(np.unique(y[np.isfinite(y)]))
        if vals <= {1.0, 2.0}:
            y = np.where(y == 2.0, 1.0, 0.0)
        elif not vals <= {0.0, 1.0}:
            raise ValueError(f"二元表型取值异常 {sorted(vals)[:5]}，无法用 Logit 重算")
    C = dat[cov_cols].to_numpy(dtype=float)
    row_lookup = {iid: k for k, iid in enumerate(raw.df["IID"].astype(str))}
    pos = np.array([row_lookup[i] for i in dat["IID"].astype(str)])
    dose = raw.df[raw.snps].to_numpy(dtype=float)[pos, :]
    snp_col = {s: j for j, s in enumerate(raw.snps)}
    assoc_idx = assoc.set_index(assoc["SNP"].astype(str))
    import statsmodels.api as sm
    deltas = {"beta": [], "se": []}
    compared: list[dict[str, Any]] = []
    failed: list[str] = []
    for snp in picks:
        if snp not in snp_col:
            failed.append(snp)
            continue
        row = assoc_idx.loc[snp]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        g = dose[:, snp_col[snp]]
        ok = np.isfinite(g)
        if ok.sum() < 20:
            failed.append(snp)
            continue
        X = np.column_stack([np.ones(int(ok.sum())), C[ok], g[ok]])
        yy = y[ok]
        try:
            if model == "logistic":
                res = sm.Logit(yy, X).fit(disp=0, method="newton", maxiter=100)
                b_hat, se_hat = float(res.params[-1]), float(res.bse[-1])
                if "OR" in assoc.columns and np.isfinite(float(row.get("OR", np.nan))):
                    or_val = float(row["OR"])
                    b_plink = math.log(or_val) if or_val > 0 else float("nan")
                else:
                    b_plink = float(row.get("BETA", np.nan))
                if "SE" in assoc.columns:
                    se_plink = float(row["SE"])
                else:
                    lo = float(row.get("L95", np.nan))
                    hi = float(row.get("U95", np.nan))
                    se_plink = ((math.log(hi) - math.log(lo)) / (2 * 1.959964)
                                if lo > 0 and hi > 0 else float("nan"))
            else:
                res = sm.OLS(yy, X).fit()
                b_hat, se_hat = float(res.params[-1]), float(res.bse[-1])
                b_plink = float(row.get("BETA", np.nan))
                se_plink = float(row.get("SE", np.nan))
        except Exception:
            failed.append(snp)
            continue
        if not (np.isfinite(b_plink) and np.isfinite(se_plink)):
            failed.append(snp)
            continue
        db, dse = abs(b_hat - b_plink), abs(se_hat - se_plink)
        deltas["beta"].append(db)
        deltas["se"].append(dse)
        compared.append({"snp": snp, "beta_plink": round(b_plink, 8),
                         "beta_recompute": round(b_hat, 8),
                         "se_plink": round(se_plink, 8),
                         "se_recompute": round(se_hat, 8)})
    if not compared:
        return _ck("gwas_beta_se_recompute", "GWAS", "warn",
                   f"未完成任何 SNP 的 beta/SE 重算（失败 {len(failed)} 个，模型 {model}）")
    max_b, max_s = max(deltas["beta"]), max(deltas["se"])
    tol = float(t.tol["beta"] if model != "logistic"
                else t.tol.get("beta_logistic", 1e-3))
    status = "pass" if max_b <= tol and max_s <= tol else (
        "warn" if max_b <= tol * 10 and max_s <= tol * 10 else "fail")
    detail = (f"{model} 模型，{len(compared)} 个随机 SNP：max|Δbeta|={max_b:.2e}，"
              f"max|ΔSE|={max_s:.2e}（容差 {tol:.0e}）")
    if failed:
        detail += f"；{len(failed)} 个 SNP 无法对照"
    return _ck("gwas_beta_se_recompute", "GWAS", status, detail,
               value={"model": model, "n_compared": len(compared),
                      "max_abs_beta_diff": max_b, "max_abs_se_diff": max_s,
                      "examples": compared[:3]},
               threshold=f"max|Δ| <= {tol}")


def check_lambda(t: _Target, found: dict) -> dict:
    assoc = _assoc_add(t, found)
    if assoc.empty or "P" not in assoc.columns:
        return _ck("lambda_gc_recompute", "GWAS", "warn", "缺少 GWAS P 值")
    p = assoc["P"].to_numpy(dtype=float)
    lam_indep = utils.lambda_gc(p)
    gwas = t.metrics.get("gwas") or {}
    comparisons: dict[str, float] = {}
    if gwas.get("lambda_gc") is not None:
        comparisons["metrics.lambda_gc"] = float(gwas["lambda_gc"])
    if gwas.get("lambda_plink") is not None:
        comparisons["metrics.lambda_plink"] = float(gwas["lambda_plink"])
    log_p = found.get("gwas_log", {}).get("path")
    if log_p is not None:
        text = Path(log_p).read_text(encoding="utf-8", errors="replace")
        m = re.search(r"lambda[^0-9]{0,40}([0-9]+\.[0-9]+)", text)
        if m:
            comparisons["plink.log.lambda"] = float(m.group(1))
    adj_p = found["adjusted"]["path"]
    if adj_p is not None:
        try:
            adj = _read_plink_table(adj_p)
            if "UNADJ" in adj.columns:
                lam_adj = utils.lambda_gc(pd.to_numeric(adj["UNADJ"],
                                                        errors="coerce").dropna())
                if np.isfinite(lam_adj):
                    comparisons["adjusted.recompute"] = float(lam_adj)
        except Exception:
            pass
    rel_tol = float(t.tol["lambda"])
    problems: list[str] = []
    for label, val in comparisons.items():
        rel = abs(float(val) - lam_indep) / max(abs(lam_indep), 1e-12)
        if rel > rel_tol:
            problems.append(f"{label}={val:.5f} 与独立重算 {lam_indep:.5f} "
                            f"相对差 {rel:.2e}（容差 {rel_tol:.0e}）")
    n_sig = int(np.sum(p < ALPHA_GW))
    detail = (f"独立重算 λGC={lam_indep:.5f}"
              + "".join(f"，{k}={v:.5f}" for k, v in comparisons.items())
              + f"；全基因组显著位点 {n_sig} 个")
    if problems:
        detail += "；" + "；".join(problems)
    return _ck("lambda_gc_recompute", "GWAS", "pass" if not problems else "fail", detail,
               value={"lambda_independent": lam_indep, "comparisons": comparisons,
                      "n_significant": n_sig},
               threshold=f"相对容差 {rel_tol}")


def check_gwas_results_table(t: _Target, found: dict) -> dict:
    results_p = found["results"]["path"]
    if results_p is None:
        return _ck("gwas_results_table", "GWAS", "warn", "缺少 gwas_results.csv")
    res = pd.read_csv(results_p, encoding="utf-8-sig")
    assoc = _assoc_add(t, found)
    leads = (t.metrics.get("gwas") or {}).get("lead_snps") or []
    problems: list[str] = []
    need = {"SNP", "BETA", "SE", "P"}
    if not need <= set(res.columns):
        problems.append(f"缺少列 {sorted(need - set(res.columns))}")
    if len(res) != len(assoc):
        problems.append(f"行数 {len(res)} != .assoc ADD 行数 {len(assoc)}")
    elif "P" in res.columns and "P" in assoc.columns:
        a = assoc.sort_values("SNP")["P"].to_numpy(dtype=float)
        b = res.sort_values("SNP")["P"].to_numpy(dtype=float)
        d = float(np.nanmax(np.abs(a - b)))
        if not np.isfinite(d) or d > 1e-4:
            problems.append(f"P 值与 .assoc 最大差 {d:.2e}")
    snp_set = set(res["SNP"].astype(str)) if "SNP" in res.columns else set()
    missing_leads = [s.get("snp") for s in leads if s.get("snp") not in snp_set]
    if missing_leads:
        problems.append(f"lead SNP 不在结果表中: {missing_leads[:3]}")
    return _ck("gwas_results_table", "GWAS", "pass" if not problems else "fail",
               f"gwas_results.csv {len(res)} 行，含 {len(leads)} 个 lead SNP"
               + ("" if not problems else "；" + "；".join(problems)),
               value={"n_rows": int(len(res)), "n_lead_snps": len(leads)},
               threshold="列齐全、行数一致、P 与 .assoc 一致、lead SNP 可追溯")


# ---------------------------------------------------------------------------
# 汇总入口
# ---------------------------------------------------------------------------
LOGGER = logging.getLogger("validate")


def run_validation(ctx_or_run_dir: Any, *, logger: logging.Logger | None = None,
                   out_path: str | Path | None = None) -> dict[str, Any]:
    """执行独立重算校验并写 ``04_calibration/validation.json``。

    Parameters
    ----------
    ctx_or_run_dir:
        ``PipelineContext``（含 run_dir/config/plink/logger/state/stage_dir），
        或 run 目录路径（``str``/``Path``）。
    logger:
        可选 logger；ctx 模式下默认用 ``ctx.logger``。
    out_path:
        可选，覆盖 validation.json 路径（默认 ``<run_dir>/04_calibration/validation.json``）。

    Returns
    -------
    dict
        ``{"checks": [...], "overall": ..., "n_checks": ..., ...}``，同时落盘。
    """
    logger = logger or LOGGER
    t = _resolve(ctx_or_run_dir, logger)
    if not t.run_dir.exists():
        raise FileNotFoundError(f"run 目录不存在: {t.run_dir}")
    logger.info("独立重算校验开始：run_dir=%s", t.run_dir)
    found = discover_artifacts(t)
    _plink(t)

    records: list[dict[str, Any]] = []
    _guard(records, "structure_counts", "结构", lambda: check_file_structure(t, found))
    _guard(records, "duplicate_ids", "结构", lambda: check_duplicate_ids(t, found))
    _guard(records, "allele_legality", "结构", lambda: check_alleles(t, found))
    _guard(records, "qc_filter_arithmetic", "QC", lambda: check_qc_arithmetic(t, found))
    _guard(records, "ld_prune_counts", "QC", lambda: check_prune_counts(t, found))
    _guard(records, "maf_recompute", "频率", lambda: check_maf(t, found))
    _guard(records, "missingness_recompute", "频率", lambda: check_missingness(t, found))
    _guard(records, "hwe_extremes", "QC", lambda: check_hwe(t, found))
    _guard(records, "phenotype_coding", "表型", lambda: check_phenotype(t, found))
    _guard(records, "pca_recompute", "群体结构", lambda: check_pca(t, found))
    _guard(records, "relatedness_range", "亲缘", lambda: check_relatedness(t, found))
    _guard(records, "kinship_matrix", "亲缘", lambda: check_kinship(t, found))
    _guard(records, "gwas_assoc_table", "GWAS", lambda: check_assoc_table(t, found))
    _guard(records, "gwas_beta_se_recompute", "GWAS", lambda: check_gwas_beta_se(t, found))
    _guard(records, "lambda_gc_recompute", "GWAS", lambda: check_lambda(t, found))
    _guard(records, "gwas_results_table", "GWAS", lambda: check_gwas_results_table(t, found))

    n_pass = sum(1 for c in records if c["status"] == "pass")
    n_warn = sum(1 for c in records if c["status"] == "warn")
    n_fail = sum(1 for c in records if c["status"] == "fail")
    overall = "fail" if n_fail else ("warn" if n_warn else "pass")
    run_id = t.state.get("run_id") or t.config.get("run_id") or t.run_dir.name
    artifacts = {k: _rel(t, Path(v)) if v else None
                 for k, v in t.artifacts.items() if not str(k).startswith("_")}
    result: dict[str, Any] = {
        "run_id": run_id,
        "generated_at": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "overall": overall,
        "n_checks": len(records),
        "n_pass": n_pass,
        "n_warn": n_warn,
        "n_fail": n_fail,
        "tolerances": t.tol,
        "artifacts": artifacts,
        "artifact_sources": t.artifacts.get("_sources") or {},
        "checks": records,
    }
    target = Path(out_path) if out_path else t.cal_dir / "validation.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    utils.write_json(target, result)
    try:
        csv_path = target.parent / "tables" / "validation_checks.csv"
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        utils.write_csv(pd.DataFrame(records), csv_path)
    except Exception as exc:  # noqa: BLE001
        logger.warning("写 validation_checks.csv 失败: %s", exc)
    for c in records:
        log = {"pass": logger.info, "warn": logger.warning, "fail": logger.error}[c["status"]]
        log("[validate] %-4s %s: %s", c["status"].upper(), c["name"], c["detail"])
    logger.info("独立重算校验完成：%d 项（pass %d / warn %d / fail %d），总体 %s",
                len(records), n_pass, n_warn, n_fail, overall)
    return result


# 兼容别名
validate = run_validation


if __name__ == "__main__":  # pragma: no cover - 手工对某个 run 目录执行
    import argparse
    import json

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description="对已有 run 目录执行独立重算校验")
    ap.add_argument("run_dir", help="plink_runs/<run_id> 目录")
    ap.add_argument("--json", action="store_true", help="打印完整 JSON")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    res = run_validation(a.run_dir)
    print(f"overall={res['overall']}  checks={res['n_checks']} "
          f"pass={res['n_pass']} warn={res['n_warn']} fail={res['n_fail']}")
    for c in res["checks"]:
        print(f"  [{c['status'].upper():4s}] {c['name']}: {c['detail']}")
    if a.json:
        print(json.dumps(res, ensure_ascii=False, indent=2))
