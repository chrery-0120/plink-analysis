"""Configuration defaults, species presets and chromosome-flag inference.

The whole pipeline is parameterised by a single dict persisted to ``run.json``;
thresholds never live in the stage modules.  See DESIGN.md section 4.
"""
from __future__ import annotations

import copy
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from . import __version__ as PIPELINE_VERSION
from . import utils

# ---------------------------------------------------------------------------
# species presets
# ---------------------------------------------------------------------------
# autosome_count = number of autosomes (PLINK --chr-set value when > 22)
# sex_system = "XY" (mammal) or "ZW" (birds)
SPECIES_PRESETS: dict[str, dict[str, Any]] = {
    "human": {"label": "人", "autosome_count": 22, "sex_system": "XY"},
    "pig": {"label": "猪", "autosome_count": 18, "sex_system": "XY"},
    "cattle": {"label": "牛", "autosome_count": 29, "sex_system": "XY"},
    "sheep": {"label": "绵羊", "autosome_count": 26, "sex_system": "XY"},
    "goat": {"label": "山羊", "autosome_count": 29, "sex_system": "XY"},
    "horse": {"label": "马", "autosome_count": 31, "sex_system": "XY"},
    "dog": {"label": "狗", "autosome_count": 38, "sex_system": "XY"},
    "mouse": {"label": "小鼠", "autosome_count": 19, "sex_system": "XY"},
    "chicken": {"label": "鸡", "autosome_count": 38, "sex_system": "ZW"},
    "duck": {"label": "鸭", "autosome_count": 40, "sex_system": "ZW"},
    "other": {"label": "其他/未知", "autosome_count": None, "sex_system": "unknown"},
}

NUMERIC_LIMIT = 95  # PLINK 1.9 --chr-set upper practical limit


def species_preset(name: str | None) -> dict[str, Any]:
    if not name:
        return copy.deepcopy(SPECIES_PRESETS["other"])
    return copy.deepcopy(SPECIES_PRESETS.get(name.lower(), SPECIES_PRESETS["other"]))


# ---------------------------------------------------------------------------
# chromosome analysis
# ---------------------------------------------------------------------------
_CHR_PREFIX = re.compile(r"^(chr)?0*(?P<num>\d+)$", re.IGNORECASE)


def normalize_chr(code: Any) -> str:
    text = str(code).strip()
    if text.upper().startswith("CHR"):
        text = text[3:]
    if text.upper() in {"M", "MT", "MITO"}:
        return "MT"
    if text.upper() in {"X", "Y", "XY", "Z", "W"}:
        return text.upper()
    match = _CHR_PREFIX.match(text)
    if match:
        return str(int(match.group("num")))
    return text


def analyze_chromosomes(chromosomes: Iterable[Any]) -> dict[str, Any]:
    """Classify the chromosome codes present in a genotype file."""
    codes = [normalize_chr(c) for c in chromosomes if str(c).strip() not in {"", "0"}]
    numeric: list[int] = []
    nonstandard: list[str] = []
    for code in codes:
        if code.isdigit():
            numeric.append(int(code))
        elif code in {"X", "Y", "Z", "W", "XY", "MT"}:
            continue
        else:
            nonstandard.append(code)
    return {
        "n_codes": len(set(codes)),
        "numeric_max": max(numeric) if numeric else None,
        "numeric_min": min(numeric) if numeric else None,
        "n_numeric": len(set(str(n) for n in numeric)),
        "has_x": "X" in codes,
        "has_y": "Y" in codes,
        "has_z": "Z" in codes,
        "has_w": "W" in codes,
        "has_mt": "MT" in codes,
        "nonstandard_codes": sorted(set(nonstandard))[:20],
        "has_nonstandard": bool(nonstandard),
        "has_chr_prefix": any(str(c).strip().lower().startswith("chr") for c in chromosomes),
    }


def infer_species(chromosomes: Iterable[Any], explicit: str = "auto") -> dict[str, Any]:
    """Infer species / chr-set / sex system from chromosome codes."""
    info = analyze_chromosomes(chromosomes)
    if explicit and explicit.lower() not in {"auto", ""}:
        preset = species_preset(explicit)
        autosomes = preset.get("autosome_count") or info.get("numeric_max")
        return {
            "name": explicit.lower(),
            "label": preset.get("label", explicit),
            "autosome_count": autosomes,
            "sex_system": preset.get("sex_system", "unknown"),
            "chr_set": autosomes if autosomes and autosomes > 22 else None,
            "allow_extra_chr": info["has_nonstandard"] or info["has_z"] or info["has_w"],
            "evidence": "user-specified",
            "chromosome_info": info,
        }

    if info["has_z"] or info["has_w"]:
        autosomes = info["numeric_max"] or species_preset("chicken")["autosome_count"]
        return {
            "name": "chicken" if autosomes >= 30 else "bird",
            "label": "鸡/家禽(ZW)",
            "autosome_count": autosomes,
            "sex_system": "ZW",
            "chr_set": autosomes if autosomes and autosomes > 22 else None,
            "allow_extra_chr": True,
            "evidence": "Z/W chromosome present",
            "chromosome_info": info,
        }

    if info["has_x"] or info["has_y"]:
        autosomes = info["numeric_max"]
        match = None
        for name, preset in SPECIES_PRESETS.items():
            if preset["autosome_count"] == autosomes:
                match = name
                break
        if match is None:
            match = "other"
        return {
            "name": match,
            "label": SPECIES_PRESETS[match]["label"],
            "autosome_count": autosomes,
            "sex_system": "XY",
            "chr_set": autosomes if autosomes and autosomes > 22 else None,
            "allow_extra_chr": info["has_nonstandard"] or autosomes is None or (autosomes > NUMERIC_LIMIT),
            "evidence": "X/Y chromosome present",
            "chromosome_info": info,
        }

    autosomes = info["numeric_max"]
    match = None
    for name, preset in SPECIES_PRESETS.items():
        if preset["autosome_count"] == autosomes:
            match = name
            break
    if match and match != "other":
        preset = SPECIES_PRESETS[match]
        return {
            "name": match,
            "label": preset["label"],
            "autosome_count": autosomes,
            "sex_system": preset["sex_system"],
            "chr_set": autosomes if autosomes and autosomes > 22 else None,
            "allow_extra_chr": info["has_nonstandard"] or (autosomes > NUMERIC_LIMIT),
            "evidence": f"常染色体数 {autosomes} 与 {preset['label']} 匹配（无性染色体）",
            "chromosome_info": info,
        }
    return {
        "name": "other",
        "label": "未知/无性染色体",
        "autosome_count": autosomes,
        "sex_system": "unknown",
        "chr_set": autosomes if autosomes and autosomes > 22 else None,
        "allow_extra_chr": info["has_nonstandard"] or (autosomes is not None and autosomes > NUMERIC_LIMIT),
        "evidence": "no sex chromosome detected",
        "chromosome_info": info,
    }


def chr_flags(species: dict[str, Any]) -> list[str]:
    """PLINK flags needed to accept this dataset's chromosome coding."""
    flags: list[str] = []
    autosomes = species.get("autosome_count") or species.get("chr_set")
    if autosomes and int(autosomes) > 22:
        if int(autosomes) <= NUMERIC_LIMIT:
            flags.extend(["--chr-set", str(int(autosomes))])
    if species.get("allow_extra_chr"):
        flags.append("--allow-extra-chr")
    return flags


# ---------------------------------------------------------------------------
# default config
# ---------------------------------------------------------------------------
def default_config() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "pipeline_version": PIPELINE_VERSION,
        "run_id": None,
        "created_at": None,
        "input": {
            "genotype": {"type": "auto", "prefix": None, "path": None},
            "phenotype": {"type": "auto", "path": None, "missing": ["-9", "NA", "NaN", "nan", "."]},
            "covariates": {"path": None, "columns": None, "categorical": []},
            "groups": {"path": None, "name": None},
            "pedigree": {"path": None},
        },
        "species": {"name": "auto", "chr_set": None, "autosome_count": None,
                    "sex_system": "auto", "nonautosomal": []},
        "qc": {
            "mind": 0.10, "geno": 0.10, "maf": 0.01, "hwe": 1e-6,
            "het_sd": 3.0, "duplicate_pihat": 0.90, "related_pihat": 0.25,
            "sex_check": "auto", "pca_outlier_sd": 6.0,
            "ld_prune": {"window": 50, "step": 5, "r2": 0.2},
            "remove_related": False, "remove_pca_outliers": False,
        },
        "popgen": {"pca": 20, "ld_window_kb": 1000, "ld_bin_kb": 50,
                   "roh": {"snp": 50, "kb": 1000, "density": 50, "gap": 1000,
                            "het_missing": 0.0}},
        "gwas": {"model": "auto", "covariates": "auto", "n_pcs": 5, "ci": 0.95,
                 "maf": 0.01, "geno": 0.10, "hwe": 1e-6,
                 "clump": {"p1": 5e-8, "p2": 1e-4, "r2": 0.2, "kb": 1000},
                 "genomic_control": "none"},
        "calibration": {"enabled": True, "n_replicates": 10, "n_causal": 20,
                        "h2": 0.5, "n_snps_subset": 20000, "alpha": 5e-8,
                        "seed": 20261023, "truth_path": None, "model": "auto"},
        "runtime": {"threads": 4, "memory_mb": 4000, "timeout_s": 7200,
                    "keep_intermediate": True},
        "report": {"language": "zh-CN", "title": None},
    }


def default_run_id() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M") + "_run"


# ---------------------------------------------------------------------------
# config loading / merging / validation
# ---------------------------------------------------------------------------
def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def validate_config(cfg: dict[str, Any]) -> tuple[list[str], list[str]]:
    """Return (errors, warnings) for a run configuration."""
    errors: list[str] = []
    warnings: list[str] = []
    qc = cfg.get("qc", {})
    gwas = cfg.get("gwas", {})
    calib = cfg.get("calibration", {})

    for key in ("mind", "geno", "maf", "hwe"):
        value = qc.get(key)
        if value is None or not (0 <= float(value) < 1):
            errors.append(f"qc.{key} 必须在 [0,1) 内，当前 {value!r}")
    if not (0 < float(qc.get("duplicate_pihat", 0.9)) <= 1):
        errors.append("qc.duplicate_pihat 必须在 (0,1] 内")
    if float(qc.get("related_pihat", 0.25)) > float(qc.get("duplicate_pihat", 0.9)):
        warnings.append("qc.related_pihat 大于 duplicate_pihat，亲缘阈值逻辑可能反了")
    if int(gwas.get("n_pcs", 0)) < 0:
        errors.append("gwas.n_pcs 不能为负")
    if calib.get("enabled", True):
        if int(calib.get("n_replicates", 0)) < 2:
            warnings.append("calibration.n_replicates < 2，校准结论会不稳定（建议 ≥ 10）")
    runtime = cfg.get("runtime", {})
    if int(runtime.get("threads", 0)) < 1:
        errors.append("runtime.threads 必须 ≥ 1")
    if int(runtime.get("memory_mb", 0)) < 256:
        errors.append("runtime.memory_mb 必须 ≥ 256")
    return errors, warnings


def load_config(path: str | Path) -> dict[str, Any]:
    return utils.load_json(path)


def save_config(cfg: dict[str, Any], path: str | Path) -> Path:
    return utils.write_json(path, cfg)


# ---------------------------------------------------------------------------
# auto-detection helpers
# ---------------------------------------------------------------------------
def detect_chromosomes_from_prefix(prefix: str | Path) -> list[str]:
    prefix = Path(prefix)
    bim = prefix.with_suffix(".bim") if prefix.suffix != ".bim" else prefix
    if bim.exists():
        return utils.read_bim(bim)["CHR"].astype(str).tolist()
    fam = prefix.with_suffix(".fam")
    if fam.exists():
        warnings = []
        return []
    return []


def detect_chromosomes_from_vcf(path: str | Path, limit: int = 200000) -> list[str]:
    import gzip
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" or path.name.endswith(".vcf.gz") else open
    chroms: list[str] = []
    seen: set[str] = set()
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            if chrom not in seen:
                seen.add(chrom)
                chroms.append(chrom)
            if len(chroms) > 5000:
                break
            limit -= 1
            if limit <= 0:
                break
    return chroms


def count_vcf_variants(path: str | Path) -> int:
    import gzip
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" or path.name.endswith(".vcf.gz") else open
    count = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("#"):
                count += 1
    return count


def vcf_sample_ids(path: str | Path) -> list[str]:
    import gzip
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" or path.name.endswith(".vcf.gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    return []


def resolve_species(cfg: dict[str, Any], chromosomes: Iterable[Any]) -> dict[str, Any]:
    explicit = cfg.get("species", {}).get("name", "auto")
    return infer_species(chromosomes, explicit if explicit else "auto")


def species_chr_flags(cfg: dict[str, Any]) -> list[str]:
    info = cfg.get("species", {})
    if not info.get("autosome_count") and not info.get("chr_set"):
        return []
    return chr_flags(info)


def describe_species(cfg: dict[str, Any]) -> str:
    sp = cfg.get("species", {})
    return (f"{sp.get('label', sp.get('name', '未知'))} "
            f"(常染色体 {sp.get('autosome_count', '未知')}, 性染色体系统 {sp.get('sex_system', '未知')})")
