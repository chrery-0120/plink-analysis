"""Shared utilities for the PLINK analysis pipeline.

Everything that touches the filesystem, the PLINK executable, or the PLINK log
format lives here.  Stage modules must call :class:`Plink` instead of running
PLINK themselves so that every command is logged and every failure carries a
diagnosis.

Contract: see DESIGN.md section 7.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

ENCODING = "utf-8"
STAGE_DIRS = {
    "intake": "00_intake",
    "qc": "01_qc",
    "popgen": "02_popgen",
    "gwas": "03_gwas",
    "calibration": "04_calibration",
    "reports": "05_reports",
}


# ---------------------------------------------------------------------------
# basic helpers
# ---------------------------------------------------------------------------
def ensure_utf8_stdout() -> None:
    """Windows consoles are often GBK; force UTF-8 so Chinese log lines survive."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")  # type: ignore[attr-defined]
        except Exception:  # pragma: no cover - non-standard stream
            pass


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def write_json(path: str | Path, obj: Any) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(_sanitize(obj), ensure_ascii=False, indent=2, default=_json_default)
    path.write_text(payload, encoding=ENCODING)
    return path


def _sanitize(value: Any) -> Any:
    """Make an object strictly JSON-serialisable (NaN/Inf -> null)."""
    if isinstance(value, dict):
        return {str(k): _sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_sanitize(v) for v in value]
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, (np.floating,)):
        number = float(value)
        return number if np.isfinite(number) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.ndarray):
        return _sanitize(value.tolist())
    if isinstance(value, (pd.Timestamp, datetime)):
        return value.isoformat()
    if isinstance(value, Path):
        return str(value)
    return value


def _json_default(obj: Any) -> Any:
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    if isinstance(obj, (pd.Timestamp, datetime)):
        return obj.isoformat()
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)


def load_json(path: str | Path, default: Any = None) -> Any:
    path = Path(path)
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding=ENCODING))


def append_jsonl(path: str | Path, obj: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding=ENCODING) as handle:
        handle.write(json.dumps(obj, ensure_ascii=False, default=_json_default) + "\n")


def write_csv(df: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return path


def read_csv(path: str | Path, **kwargs: Any) -> pd.DataFrame:
    return pd.read_csv(path, encoding="utf-8-sig", **kwargs)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_text_robust(path: str | Path) -> str:
    """Read a text file that may be UTF-8 (with/without BOM), GBK or Latin-1."""
    raw = Path(path).read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "gb18030", "latin-1"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def setup_logger(log_path: str | Path, level: int = logging.INFO) -> logging.Logger:
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"plink_pipeline.{log_path.stem}.{id(log_path)}")
    logger.setLevel(level)
    logger.propagate = False
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%Y-%m-%d %H:%M:%S")
    file_handler = logging.FileHandler(log_path, encoding=ENCODING)
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(stream_handler)
    return logger


def stage_dir(run_dir: str | Path, stage: str) -> Path:
    """Return (and create) the directory for a stage, including subfolders."""
    run_dir = Path(run_dir)
    folder = STAGE_DIRS.get(stage, stage)
    path = run_dir / folder
    for sub in ("tables", "plots", "data"):
        (path / sub).mkdir(parents=True, exist_ok=True)
    return path


def rel_to_run(run_dir: str | Path, path: str | Path) -> str:
    try:
        return Path(path).resolve().relative_to(Path(run_dir).resolve()).as_posix()
    except ValueError:
        return str(path)


def fmt_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, sec = divmod(seconds, 60)
    return f"{int(minutes)}m{sec:.0f}s"


# ---------------------------------------------------------------------------
# PLINK invocation
# ---------------------------------------------------------------------------
@dataclass
class PlinkResult:
    args: list[str]
    out_prefix: str | None
    returncode: int
    elapsed_s: float
    log_path: str | None
    log_text: str
    counts: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    stdout: str = ""

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and not self.errors


class PlinkError(RuntimeError):
    def __init__(self, message: str, *, command: str, returncode: int,
                 log_path: str | None, log_tail: str, diagnosis: list[str]):
        super().__init__(message)
        self.command = command
        self.returncode = returncode
        self.log_path = log_path
        self.log_tail = log_tail
        self.diagnosis = diagnosis

    def user_message(self) -> str:
        lines = [f"PLINK 执行失败: {self}"]
        if self.command:
            lines.append(f"命令: {self.command}")
        if self.log_path:
            lines.append(f"日志: {self.log_path}")
        if self.diagnosis:
            lines.append("可能原因与修复建议:")
            lines.extend(f"  - {item}" for item in self.diagnosis)
        if self.log_tail:
            lines.append("日志末尾:")
            lines.extend("  " + line for line in self.log_tail.strip().splitlines()[-8:])
        return "\n".join(lines)


_COUNT_PATTERNS: dict[str, list[re.Pattern[str]]] = {
    "n_variants_loaded": [re.compile(r"(\d+) variants? loaded from \.bim file"),
                          re.compile(r"(\d+) markers? loaded from \.map file")],
    "n_samples_loaded": [re.compile(r"(\d+) (?:people|samples) .*loaded from \.fam file"),
                         re.compile(r"(\d+) (?:people|samples) loaded from \.fam file")],
    "n_pass_qc": [re.compile(r"(\d+) variants? and (\d+) (?:people|samples) pass filters and QC")],
    "n_samples_pass_qc": [re.compile(r"(\d+) (?:people|samples) pass filters and QC")],
    "n_variants_remaining": [re.compile(r"(\d+) variants? remaining after main filters")],
    "n_samples_remaining": [re.compile(r"(\d+) (?:people|samples) remaining after main filters")],
    "genotyping_rate": [re.compile(r"Total genotyping rate is ([\d.]+)")],
    "n_variants_removed_missing": [re.compile(r"(\d+) variants? removed due to missing genotype data")],
    "n_samples_removed_missing": [re.compile(r"(\d+) (?:people|samples) removed due to missing genotype data")],
    "n_variants_removed_hwe": [re.compile(r"(\d+) variants? removed due to Hardy-Weinberg exact test")],
    "n_variants_removed_maf": [re.compile(r"(\d+) variants? removed due to minor allele threshold")],
    "n_samples_removed_het": [re.compile(r"(\d+) (?:people|samples) removed due to heterozygosity")],
    "n_cases": [re.compile(r"(\d+) cases? and (\d+) controls?")],
    "lambda_gc": [re.compile(r"Genomic control lambda(?: estimate)?:? ([\d.]+)", re.IGNORECASE),
                  re.compile(r"Genomic inflation est\.? lambda[^=]*=\s*([\d.]+)", re.IGNORECASE),
                  re.compile(r"lambda(?: GC)?[^0-9]{0,20}([\d.]+)", re.IGNORECASE)],
}


def parse_plink_log(text: str) -> dict[str, Any]:
    """Extract counts, warnings, notes and errors from a PLINK .log body."""
    counts: dict[str, Any] = {}
    for key, patterns in _COUNT_PATTERNS.items():
        for pattern in patterns:
            match = pattern.search(text)
            if match:
                groups = match.groups()
                if len(groups) == 1:
                    value: Any = groups[0]
                else:
                    value = tuple(groups)
                if key in {"n_pass_qc", "n_cases"}:
                    counts[key] = [int(v) for v in groups]
                else:
                    try:
                        counts[key] = float(value) if "." in str(value) else int(value)
                    except (TypeError, ValueError):
                        counts[key] = value
                break

    warnings: list[str] = []
    errors: list[str] = []
    notes: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        low = line.lower()
        if low.startswith("error:") or low.startswith("error "):
            errors.append(line)
        elif low.startswith("warning:") or low.startswith("warning "):
            warnings.append(line)
        elif low.startswith("note:") or low.startswith("note "):
            notes.append(line)
    return {"counts": counts, "warnings": _dedupe(warnings),
            "errors": _dedupe(errors), "notes": _dedupe(notes), "text": text}


def _dedupe(items: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


# keyword -> list of Chinese suggestions
_DIAGNOSIS_RULES: list[tuple[tuple[str, ...], list[str]]] = [
    (("invalid chromosome code",),
     ["染色体命名超出 PLINK 默认 chr-set：用 --chr-set <该物种常染色体数>（如牛 29、鸡 38）",
      "若含 Z/W 或 scaffold 名称，再加 --allow-extra-chr；两者可同时使用（PLINK 1.9 实测可用）",
      "也可先用 --allow-extra-chr 读入，再用 --update-chr 统一命名"]),
    (("duplicate id", "duplicate sample", "duplicated id"),
     ["样本或位点 ID 重复：检查 .fam/.bim；VCF 导入可加 --double-id 或 --const-fid",
      "重复位点用 --list-duplicate-vars 查看后 --exclude 去除"]),
    (("no variants remaining", "no people remaining", "no samples remaining"),
     ["过滤过严或输入本身为空：先看 .log 中三条过滤计数，再放宽 --geno/--mind/--maf/--hwe",
      "若位点全被 HWE 过滤，检查是否把 case/control 或品种混在一起做 HWE（HWE 只适合在对照/单一群体内）"]),
    (("failed to open", "cannot open", "no such file"),
     ["路径/权限/文件占用：确认文件存在、未被 Excel 打开、路径无非法字符；Windows 反斜杠用引号包裹",
      "相对路径以 run 目录为基准，建议全部传绝对路径"]),
    (("magic number", "not a valid .bed", "bed file has wrong"),
     ["不是 PLINK 1 二进制 .bed（可能把 .ped 或 VCF 改名为 .bed）：重新用 --make-bed 生成",
      "样本主序 .bed 不被支持，需用 PLINK 1.9 重新转换"]),
    (("invalid phenotype", "phenotype file", "phenotype values"),
     ["表型文件格式/编码问题：PLINK 期望 FID IID PHENO；二元性状必须编码为 1=对照 2=病例，缺失用 -9",
      "检查缺失值写法（NA/0）与分隔符（空格/Tab），并确认 FID/IID 与 .fam 完全一致"]),
    (("--covar", "covariate"),
     ["协变量必须全为数值且缺失用 -9；分类变量需先哑变量编码",
      "确认 --covar-number 列号与文件列顺序一致（FID/IID 之后从 1 开始）"]),
    (("multiallelic", "more than 2 alleles", "non-biallelic"),
     ["VCF 含多等位位点：加 --biallelic-only strict 或先 --vcf-half-call m 再拆分处理",
      "如需保留多等位，改用 PLINK 2 或 bcftools 拆分后再分析"]),
    (("memory", "out of memory", "cannot allocate"),
     ["内存不足：降低 --memory 或分染色体/分块运行；减少 --r2/--genome 的窗口",
      "关闭其他占用内存的程序"]),
    (("thread", "--threads"),
     ["线程参数不被该命令支持：去掉 --threads 重试（部分 PLINK 1.9 命令只支持单线程）"]),
    (("half-call", "vcf"),
     ["VCF 基因型含 0/. 等半缺失：加 --vcf-half-call m（按缺失处理）或 r（按参考等位基因）"]),
    (("split chromosome",),
     ["位点在 .bim 中未按染色体/位置排序：先用 --make-bed 重排，或检查输入是否混有不同排序", ]),
    (("sex", "y chromosome"),
     ["性染色体与性别不一致：哺乳动物用 --check-sex 检查；家禽（ZW）不支持 PLINK 的 XY 性别判定，应关闭性染色体步骤",
      "若数据来自混合物种，先按物种拆分再分析"]),
]


def diagnose_log(text: str) -> list[str]:
    low = text.lower()
    suggestions: list[str] = []
    for keywords, items in _DIAGNOSIS_RULES:
        if any(keyword in low for keyword in keywords):
            suggestions.extend(items)
    return _dedupe(suggestions)


class Plink:
    """Thin, logging wrapper around a PLINK executable."""

    def __init__(self, exe: str | Path, run_dir: str | Path, *,
                 threads: int = 4, memory_mb: int = 4000,
                 logger: logging.Logger | None = None, timeout_s: int = 7200):
        self.exe = str(Path(exe))
        self.run_dir = Path(run_dir)
        self.threads = int(threads or 0)
        self.memory_mb = int(memory_mb or 0)
        self.logger = logger or logging.getLogger("plink")
        self.timeout_s = int(timeout_s)
        self.commands_log = self.run_dir / "logs" / "commands.jsonl"
        self.plink_log_dir = self.run_dir / "logs" / "plink"
        self.plink_log_dir.mkdir(parents=True, exist_ok=True)
        self.version = self._probe_version()

    # -- public API ------------------------------------------------------
    def run(self, args: Sequence[str], *, stage: str, name: str,
            out_prefix: str | Path | None = None,
            expect: Sequence[str] = (), timeout_s: int | None = None,
            check: bool = True, use_thread_flags: bool = True) -> PlinkResult:
        argv: list[str] = [self.exe, *[str(a) for a in args]]
        if out_prefix is not None:
            out_prefix = Path(out_prefix)
            out_prefix.parent.mkdir(parents=True, exist_ok=True)
            if "--out" not in argv:
                argv.extend(["--out", str(out_prefix)])
        if use_thread_flags and out_prefix is not None:
            if self.threads > 0 and "--threads" not in argv:
                argv.extend(["--threads", str(self.threads)])
            if self.memory_mb > 0 and "--memory" not in argv:
                argv.extend(["--memory", str(int(self.memory_mb))])

        started = time.time()
        try:
            proc = subprocess.run(argv, cwd=str(self.run_dir), capture_output=True,
                                  text=True, encoding=ENCODING, errors="replace",
                                  timeout=timeout_s or self.timeout_s)
            stdout, stderr = proc.stdout or "", proc.stderr or ""
            returncode = proc.returncode
        except subprocess.TimeoutExpired as exc:
            elapsed = time.time() - started
            self._log_command(stage, name, argv, out_prefix, returncode=-1, elapsed=elapsed,
                              ok=False, error="timeout")
            raise PlinkError("超时", command=_quote(argv), returncode=-1,
                             log_path=None, log_tail="",
                             diagnosis=[f"命令超过 {timeout_s or self.timeout_s}s：减少数据量、"
                                        "分染色体运行或调大 runtime.timeout_s"]) from exc
        elapsed = time.time() - started

        log_path = Path(str(out_prefix) + ".log") if out_prefix is not None else None
        log_text = ""
        if log_path and log_path.exists():
            log_text = log_path.read_text(encoding=ENCODING, errors="replace")
        if not log_text:
            log_text = stdout + "\n" + stderr

        parsed = parse_plink_log(log_text)
        errors = list(parsed["errors"])
        if returncode != 0 and not errors:
            tail = (stderr or stdout).strip().splitlines()[-3:]
            errors = [f"PLINK 返回码 {returncode}: " + " / ".join(tail)] if tail else [f"PLINK 返回码 {returncode}"]
        diagnosis = diagnose_log(log_text)

        result = PlinkResult(args=argv, out_prefix=str(out_prefix) if out_prefix else None,
                             returncode=returncode, elapsed_s=elapsed,
                             log_path=str(log_path) if log_path else None,
                             log_text=log_text, counts=parsed["counts"],
                             warnings=parsed["warnings"], errors=errors, stdout=stdout)

        missing = [str(Path(str(out_prefix) + suffix)) for suffix in expect
                   if not Path(str(out_prefix) + suffix).exists()]
        if missing and check:
            errors.append("缺少预期输出: " + ", ".join(missing))
        ok = returncode == 0 and not errors
        self._log_command(stage, name, argv, out_prefix, returncode=returncode,
                          elapsed=elapsed, ok=ok, error="; ".join(errors) or None,
                          warnings=len(result.warnings), counts=parsed["counts"])
        if not ok and check:
            tail = "\n".join(log_text.strip().splitlines()[-10:])
            raise PlinkError("; ".join(errors[:3]) or "PLINK 命令失败",
                             command=_quote(argv), returncode=returncode,
                             log_path=str(log_path) if log_path else None,
                             log_tail=tail, diagnosis=diagnosis)
        for warning in result.warnings[:20]:
            if "BLAS/LAPACK" in warning:
                self.logger.debug("PLINK 提示 [%s/%s]: %s", stage, name, warning)
                continue
            self.logger.warning("PLINK 警告 [%s/%s]: %s", stage, name, warning)
        return result

    # -- internals -------------------------------------------------------
    def _probe_version(self) -> str:
        try:
            proc = subprocess.run([self.exe, "--version"], capture_output=True, text=True,
                                  timeout=60, encoding=ENCODING, errors="replace")
            text = (proc.stdout or proc.stderr).strip().splitlines()
            return text[0] if text else "unknown"
        except Exception:  # pragma: no cover
            return "unknown"

    def _log_command(self, stage: str, name: str, argv: Sequence[str],
                     out_prefix: Path | None, *, returncode: int, elapsed: float,
                     ok: bool, error: str | None = None, warnings: int = 0,
                     counts: dict[str, Any] | None = None) -> None:
        entry = {
            "time": now_iso(), "stage": stage, "name": name,
            "command": _quote(argv), "out_prefix": str(out_prefix) if out_prefix else None,
            "returncode": returncode, "elapsed_s": round(elapsed, 3), "ok": ok,
            "error": error, "n_warnings": warnings, "counts": counts or {},
        }
        append_jsonl(self.commands_log, entry)


def _quote(argv: Sequence[str]) -> str:
    out = []
    for token in argv:
        if any(ch.isspace() for ch in token) or '"' in token:
            out.append('"' + token.replace('"', '\\"') + '"')
        else:
            out.append(token)
    return " ".join(out)


# ---------------------------------------------------------------------------
# statistics helpers shared by stages / calibration
# ---------------------------------------------------------------------------
def lambda_gc(pvals: Iterable[float]) -> float:
    """Genomic inflation factor from p-values (median chi-square / 0.4549)."""
    arr = np.asarray([float(p) for p in pvals if p is not None and np.isfinite(float(p)) and 0 < float(p) <= 1])
    if arr.size == 0:
        return float("nan")
    from scipy.stats import chi2
    chi2_values = chi2.isf(arr, 1)
    return float(np.median(chi2_values) / 0.4549364)


def bh_fdr(pvals: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values)."""
    arr = np.asarray([float(p) for p in pvals], dtype=float)
    n = arr.size
    if n == 0:
        return arr
    order = np.argsort(arr)
    ranked = arr[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def safe_div(numerator: float, denominator: float, default: float = 0.0) -> float:
    try:
        if denominator in (0, None) or (isinstance(denominator, float) and np.isnan(denominator)):
            return default
        return float(numerator) / float(denominator)
    except Exception:
        return default


# ---------------------------------------------------------------------------
# PLINK file readers
# ---------------------------------------------------------------------------
def read_bim(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", header=None,
                       names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                       dtype={"CHR": str, "SNP": str, "A1": str, "A2": str})


def read_fam(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", header=None,
                       names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
                       dtype={"FID": str, "IID": str, "PAT": str, "MAT": str})


def read_frq(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", comment=None)
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def read_lmiss(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+")
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def read_imiss(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+")
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def read_hwe(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+")
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def read_het(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+")
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def read_raw(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", na_values=["NA", "-9"])


def read_assoc(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+")
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def read_eigenvec(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None)
    n_pc = df.shape[1] - 2
    df.columns = ["FID", "IID"] + [f"PC{i + 1}" for i in range(n_pc)]
    return df


def read_eigenval(path: str | Path) -> np.ndarray:
    return pd.read_csv(path, sep=r"\s+", header=None)[0].to_numpy(dtype=float)


def resolve_base_prefix(run_dir: str | Path) -> Path:
    """Read the working PLINK prefix written by the intake stage."""
    pointer = Path(run_dir) / "00_intake" / "base_prefix.txt"
    if not pointer.exists():
        raise FileNotFoundError(f"缺少 {pointer}：请先运行 intake 阶段")
    prefix = pointer.read_text(encoding=ENCODING).strip()
    if not prefix:
        raise ValueError("base_prefix.txt 为空")
    return Path(prefix)
