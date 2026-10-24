"""Stage 0: input discovery, format conversion, preflight checks and manifest.

Outputs (relative to the run directory):
    00_intake/manifest.json
    00_intake/preflight.json
    00_intake/metrics.json
    00_intake/base_prefix.txt          # absolute PLINK prefix used by later stages
    00_intake/data/base.{bed,bim,fam}  # normalised genotype (only when conversion needed)
    00_intake/data/pheno_clean.txt     # FID IID PHENO (binary recoded 1/2, missing -9)
    00_intake/data/covar_clean.txt     # FID IID numeric covariates
    00_intake/data/group_clean.txt     # FID IID GROUP
    00_intake/data/within.txt          # FID IID CLUSTER (for PLINK --within)
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import utils
from .config import (chr_flags, detect_chromosomes_from_vcf, infer_species,
                     resolve_species)

MAX_HASH_BYTES = 1024 ** 3  # do not hash files larger than 1 GiB
TEXT_ENCODINGS = ("utf-8-sig", "utf-8", "gb18030", "latin-1")


# ---------------------------------------------------------------------------
# discovery
# ---------------------------------------------------------------------------
def discover_inputs(path: str | Path) -> dict[str, Any]:
    """Locate a PLINK-readable genotype dataset from a file, prefix or directory."""
    path = Path(path)
    if not path.exists():
        # allow passing a PLINK prefix without extension
        bed = Path(str(path) + ".bed")
        if bed.exists():
            path = bed
        else:
            raise FileNotFoundError(f"输入路径不存在: {path}")

    candidates: list[dict[str, Any]] = []

    def add(kind: str, files: dict[str, Path]):
        if files:
            candidates.append({"type": kind, "files": {k: str(v) for k, v in files.items()}})

    def bfile_candidates(folder: Path):
        for bed in sorted(folder.glob("*.bed")):
            prefix = bed.with_suffix("")
            bim, fam = prefix.with_suffix(".bim"), prefix.with_suffix(".fam")
            if bim.exists() and fam.exists():
                add("bed", {"bed": bed, "bim": bim, "fam": fam})

    if path.is_dir():
        bfile_candidates(path)
        for ped in sorted(path.glob("*.ped")):
            map_file = ped.with_suffix(".map")
            if map_file.exists():
                add("ped", {"ped": ped, "map": map_file})
        for tped in sorted(path.glob("*.tped")):
            tfam = tped.with_suffix(".tfam")
            if tfam.exists():
                add("tped", {"tped": tped, "tfam": tfam})
        for vcf in sorted(list(path.glob("*.vcf")) + list(path.glob("*.vcf.gz"))):
            add("vcf", {"vcf": vcf})
        for gen in sorted(path.glob("*.gen")) + sorted(path.glob("*.gen.gz")):
            sample = gen.with_suffix(gen.suffix + ".sample") if gen.suffix == ".gz" else gen.with_suffix(".sample")
            if not sample.exists():
                sample = gen.with_name(gen.name.split(".")[0] + ".sample")
            if sample.exists():
                add("gen", {"gen": gen, "sample": sample})
        for txt in sorted(path.glob("*.txt")):
            try:
                first = txt.read_text(encoding="utf-8", errors="replace").splitlines()[:2]
            except Exception:
                continue
            if first and "rsid" in first[0].lower() and ("genotype" in first[0].lower()):
                add("23andme", {"23andme": txt})
    else:
        suffix = path.suffix.lower()
        if suffix == ".bed":
            prefix = path.with_suffix("")
            add("bed", {"bed": path, "bim": prefix.with_suffix(".bim"), "fam": prefix.with_suffix(".fam")})
        elif suffix == ".vcf":
            add("vcf", {"vcf": path})
        elif path.name.endswith(".vcf.gz"):
            add("vcf", {"vcf": path})
        elif suffix == ".ped":
            add("ped", {"ped": path, "map": path.with_suffix(".map")})
        elif suffix == ".tped":
            add("tped", {"tped": path, "tfam": path.with_suffix(".tfam")})
        elif suffix in {".gen", ".gz"} and path.name.endswith(".gen"):
            add("gen", {"gen": path, "sample": path.with_suffix(".sample")})
        else:
            # treat as a basename prefix
            prefix = path
            bed = Path(str(prefix) + ".bed")
            if bed.exists():
                add("bed", {"bed": bed, "bim": Path(str(prefix) + ".bim"), "fam": Path(str(prefix) + ".fam")})
            else:
                raise FileNotFoundError(f"无法从 {path} 识别基因型格式（支持 bed/ped/tped/vcf/gen/23andMe）")

    if not candidates:
        raise FileNotFoundError(f"在 {path} 下没有找到可识别的基因型数据（bed/bim/fam、ped/map、vcf 等）")
    # prefer PLINK binary, then VCF, then text formats
    order = {"bed": 0, "vcf": 1, "ped": 2, "tped": 3, "gen": 4, "23andme": 5}
    candidates.sort(key=lambda item: order.get(item["type"], 9))
    chosen = candidates[0]
    chosen["all_candidates"] = [
        {"type": c["type"], "files": c["files"]} for c in candidates
    ]
    for key, file in chosen["files"].items():
        if not Path(file).exists():
            raise FileNotFoundError(f"缺少配套文件: {file}")
    return chosen


# ---------------------------------------------------------------------------
# text helpers
# ---------------------------------------------------------------------------
def read_table_auto(path: str | Path, min_cols: int = 2) -> tuple[pd.DataFrame, bool]:
    """Read a whitespace/tab table, detect whether the first row is a header."""
    text = utils.read_text_robust(path)
    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"文件为空: {path}")
    delimiter = "\t" if "\t" in lines[0] else r"\s+"
    tokens = re.split(delimiter, lines[0].strip())
    header = False
    if len(tokens) >= min_cols:
        header = any(tok.upper() in {"FID", "IID", "#FID", "ID", "PHENO", "PHENOTYPE", "SEX"}
                     or not re.fullmatch(r"[-+]?\d*\.?\d+([eE][-+]?\d+)?", tok)
                     for tok in tokens[:3])
        # A first row whose 3rd token is numeric is data, not a header.
        if len(tokens) >= 3:
            try:
                float(tokens[2])
                header = False
            except ValueError:
                pass
    if header:
        df = pd.read_csv(path, sep=delimiter, comment="#", dtype=str,
                         encoding="utf-8-sig", engine="python")
    else:
        df = pd.read_csv(path, sep=delimiter, header=None, dtype=str,
                         comment="#", encoding="utf-8-sig", engine="python")
        df.columns = [f"V{i + 1}" for i in range(df.shape[1])]
    df.columns = [str(c).strip().lstrip("#") for c in df.columns]
    return df, header


def _find_col(df: pd.DataFrame, names: tuple[str, ...]) -> str | None:
    upper = {str(c).upper(): c for c in df.columns}
    for name in names:
        if name.upper() in upper:
            return upper[name.upper()]
    return None


def _to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str).str.strip(), errors="coerce")


# ---------------------------------------------------------------------------
# phenotype / covariate / group loading
# ---------------------------------------------------------------------------
def load_phenotype(path: Path, fam: pd.DataFrame,
                   missing_tokens: list[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    df, header = read_table_auto(path)
    fid_col = _find_col(df, ("FID", "FAMILY", "FAMILYID"))
    iid_col = _find_col(df, ("IID", "ID", "INDIVIDUAL"))
    value_cols = [c for c in df.columns if c not in {fid_col, iid_col}]
    if iid_col is None:
        # assume FID IID PHENO order (or IID PHENO)
        if df.shape[1] >= 3:
            fid_col, iid_col = df.columns[0], df.columns[1]
            value_cols = list(df.columns[2:])
        else:
            iid_col = df.columns[0]
            value_cols = list(df.columns[1:])
            fid_col = None
    if not value_cols:
        raise ValueError("表型文件没有表型列")
    pheno_col = value_cols[0]
    if len(value_cols) > 1:
        pass

    table = pd.DataFrame({
        "FID": df[fid_col].astype(str).str.strip() if fid_col else None,
        "IID": df[iid_col].astype(str).str.strip(),
        "PHENO_RAW": df[pheno_col].astype(str).str.strip(),
    })
    if fid_col is None:
        table["FID"] = None

    fam_ids = set(zip(fam["FID"].astype(str), fam["IID"].astype(str)))
    fam_iid = {iid: fid for fid, iid in fam_ids}
    matched = 0
    unmatched_rows = 0
    if table["FID"].notna().all() and table["FID"].ne("").all():
        table["key"] = list(zip(table["FID"], table["IID"]))
        matched = sum(key in fam_ids for key in table["key"])
        unmatched_rows = len(table) - matched
        if matched == 0:
            # fall back to IID-only matching
            table["FID"] = table["IID"].map(fam_iid)
            table["key"] = list(zip(table["FID"], table["IID"]))
            matched = int(table["key"].isin(fam_ids).sum())
    else:
        table["FID"] = table["IID"].map(fam_iid)
        table["key"] = list(zip(table["FID"], table["IID"]))
        matched = int(table["key"].isin(fam_ids).sum())

    missing_tokens_lower = {t.lower() for t in missing_tokens} | {"", "na", "nan", "none", "null"}
    values = table["PHENO_RAW"].str.lower().isin(missing_tokens_lower)
    numeric = pd.to_numeric(table["PHENO_RAW"].where(~values), errors="coerce")
    numeric = numeric.mask(numeric <= -9)
    raw_values = numeric.dropna()

    info: dict[str, Any] = {
        "file": str(path), "header": header, "column": pheno_col,
        "n_rows": int(len(table)), "n_matched": int(matched),
        "n_unmatched": int(len(table) - matched), "n_missing": int(values.sum()),
        "coding": "unknown", "type": "unknown",
        "unique_values": sorted({str(v) for v in raw_values.round(6).tolist()})[:10],
    }

    pheno = pd.Series(np.nan, index=range(len(fam)), dtype=float)
    fam_keys = list(zip(fam["FID"].astype(str), fam["IID"].astype(str)))
    lookup = {}
    for _, row in table.iterrows():
        raw = row["PHENO_RAW"].strip().lower()
        if raw in missing_tokens_lower:
            continue
        try:
            value = float(row["PHENO_RAW"])
        except ValueError:
            continue
        lookup[(str(row["FID"]), str(row["IID"]))] = value

    for idx, key in enumerate(fam_keys):
        if key in lookup:
            pheno.iloc[idx] = lookup[key]

    valid = pheno.dropna()
    if valid.empty:
        info["type"] = "none"
        info["coding"] = "all-missing"
    else:
        uniques = set(np.round(valid.to_numpy(), 6).tolist())
        if uniques <= {1.0, 2.0}:
            info["type"] = "binary"
            info["coding"] = "1=control,2=case"
        elif uniques <= {0.0, 1.0}:
            info["type"] = "binary"
            info["coding"] = "0/1 -> 1=control,2=case (auto-converted)"
            pheno = pheno.map(lambda v: 1.0 if v == 0 else (2.0 if v == 1 else v))
            valid = pheno.dropna()
        elif len(uniques) == 2:
            info["type"] = "binary"
            info["coding"] = f"two-level numeric {sorted(uniques)} treated as case/control"
        else:
            info["type"] = "quantitative"
            info["coding"] = "continuous numeric"

    out = pd.DataFrame({"FID": fam["FID"].astype(str), "IID": fam["IID"].astype(str),
                        "PHENO": pheno.fillna(-9).astype(float)})
    if info["type"] == "binary":
        counts = out["PHENO"].value_counts().to_dict()
        info["n_cases"] = int(counts.get(2.0, 0))
        info["n_controls"] = int(counts.get(1.0, 0))
        info["n_missing"] = int((out["PHENO"] == -9).sum())
    return out, info


def load_covariates(path: Path, fam: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    df, header = read_table_auto(path)
    fid_col = _find_col(df, ("FID", "FAMILY"))
    iid_col = _find_col(df, ("IID", "ID", "INDIVIDUAL"))
    if iid_col is None:
        if df.shape[1] >= 3:
            fid_col, iid_col = df.columns[0], df.columns[1]
        else:
            iid_col, fid_col = df.columns[0], None
    value_cols = [c for c in df.columns if c not in {fid_col, iid_col}]

    fam_ids = set(zip(fam["FID"].astype(str), fam["IID"].astype(str)))
    fam_iid = {iid: fid for fid, iid in fam_ids}
    table = pd.DataFrame({"IID": df[iid_col].astype(str).str.strip()})
    table["FID"] = (df[fid_col].astype(str).str.strip() if fid_col is not None
                    else table["IID"].map(fam_iid))
    table["key"] = list(zip(table["FID"], table["IID"]))
    matched = int(table["key"].isin(fam_ids).sum())

    out = pd.DataFrame({"FID": fam["FID"].astype(str), "IID": fam["IID"].astype(str)})
    info: dict[str, Any] = {"file": str(path), "n_matched": matched,
                            "n_rows": int(len(table)), "columns": [], "categorical": [],
                            "missing_filled": {}}
    series_by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for idx, row in df.iterrows():
        key = (str(table.loc[idx, "FID"]), str(table.loc[idx, "IID"]))
        series_by_key[key] = {col: row[col] for col in value_cols}

    for col in value_cols:
        values = df[col].astype(str).str.strip()
        numeric = _to_numeric(values)
        n_numeric = int(numeric.notna().sum())
        name = str(col).strip()
        if n_numeric >= max(2, int(0.8 * len(values))):
            out[name] = [series_by_key.get((f, i), {}).get(col) for f, i in zip(out["FID"], out["IID"])]
            out[name] = _to_numeric(out[name])
            missing = out[name].isna()
            if missing.any():
                out[name] = out[name].fillna(out[name].mean())
                info["missing_filled"][name] = int(missing.sum())
            info["columns"].append(name)
        else:
            # categorical -> dummy code, drop first level
            categories = list(pd.unique(values))
            for level in categories[1:]:
                dummy_name = f"{name}_{level}"
                out[dummy_name] = [
                    1.0 if str(series_by_key.get((f, i), {}).get(col)) == str(level) else 0.0
                    for f, i in zip(out["FID"], out["IID"])
                ]
                info["columns"].append(dummy_name)
            info["categorical"].append(name)
    return out, info


def load_groups(path: Path, fam: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    df, header = read_table_auto(path)
    fid_col = _find_col(df, ("FID", "FAMILY"))
    iid_col = _find_col(df, ("IID", "ID", "INDIVIDUAL"))
    if iid_col is None:
        if df.shape[1] >= 3:
            fid_col, iid_col = df.columns[0], df.columns[1]
            group_col = df.columns[2]
        else:
            iid_col, group_col = df.columns[0], df.columns[1]
            fid_col = None
    else:
        candidates = [c for c in df.columns if c not in {fid_col, iid_col}]
        if not candidates:
            raise ValueError("分组文件缺少分组列")
        group_col = candidates[0]

    fam_ids = set(zip(fam["FID"].astype(str), fam["IID"].astype(str)))
    fam_iid = {iid: fid for fid, iid in fam_ids}
    table = pd.DataFrame({"IID": df[iid_col].astype(str).str.strip(),
                          "GROUP": df[group_col].astype(str).str.strip()})
    table["FID"] = (df[fid_col].astype(str).str.strip() if fid_col is not None
                    else table["IID"].map(fam_iid))
    table["key"] = list(zip(table["FID"], table["IID"]))
    matched = int(table["key"].isin(fam_ids).sum())
    lookup = {key: group for key, group in zip(table["key"], table["GROUP"])}
    out = pd.DataFrame({"FID": fam["FID"].astype(str), "IID": fam["IID"].astype(str),
                        "GROUP": [lookup.get((f, i), "") for f, i in zip(fam["FID"], fam["IID"])]})
    levels = out.loc[out["GROUP"] != "", "GROUP"].value_counts()
    info = {"file": str(path), "column": str(group_col), "n_matched": matched,
            "n_levels": int(len(levels)), "levels": [
                {"name": str(k), "n": int(v)} for k, v in levels.items()]}
    return out, info


# ---------------------------------------------------------------------------
# genotype conversion / normalisation
# ---------------------------------------------------------------------------
def convert_to_bed(ctx, discovery: dict[str, Any]) -> Path:
    """Convert any supported input to a PLINK 1 binary fileset in 00_intake/data."""
    cfg = ctx.config
    data_dir = ctx.stage_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    plink: utils.Plink = ctx.plink
    kind = discovery["type"]
    files = discovery["files"]

    if kind == "bed":
        return Path(files["bed"]).with_suffix("")

    out = data_dir / "converted"
    if kind == "vcf":
        vcf = Path(files["vcf"])
        chroms = detect_chromosomes_from_vcf(vcf)
        species = resolve_species(cfg, chroms)
        flags = []
        if species.get("chr_set"):
            flags += ["--chr-set", str(species["chr_set"])]
        if species.get("allow_extra_chr"):
            flags.append("--allow-extra-chr")
        base_args = ["--vcf", str(vcf), "--double-id", "--vcf-half-call", "m",
                     "--set-missing-var-ids", "@:#", "--keep-allele-order",
                     "--biallelic-only", "strict", *flags]
        try:
            result = plink.run([*base_args, "--make-bed"], stage="intake", name="vcf_to_bed",
                               out_prefix=out, expect=[".bed", ".bim", ".fam"])
            n_multi = _log_count(result.log_text, r"(\d+) variants? (?:removed|excluded).*allele")
            if n_multi:
                ctx.logger.warning("VCF 多等位位点 %d 个已按 strict 模式剔除（避免静默丢等位基因）", n_multi)
        except utils.PlinkError as exc:
            if "allele" in str(exc).lower() or "multiallelic" in str(exc).lower():
                ctx.logger.warning("VCF 多等位位点处理失败，改用 --biallelic-only strict 重试")
                plink.run([*base_args, "--make-bed"], stage="intake", name="vcf_to_bed_strict",
                          out_prefix=out, expect=[".bed", ".bim", ".fam"])
            else:
                raise
    elif kind in {"ped", "tped"}:
        if kind == "ped":
            args = ["--file", str(Path(files["ped"]).with_suffix(""))]
        else:
            args = ["--tfile", str(Path(files["tped"]).with_suffix(""))]
        plink.run([*args, "--make-bed", "--keep-allele-order"], stage="intake",
                  name=f"{kind}_to_bed", out_prefix=out,
                  expect=[".bed", ".bim", ".fam"])
    elif kind == "gen":
        gen = files["gen"]
        sample = files["sample"]
        try:
            plink.run(["--gen", str(gen), "--sample", str(sample), "--make-bed"],
                      stage="intake", name="gen_to_bed", out_prefix=out,
                      expect=[".bed", ".bim", ".fam"])
        except utils.PlinkError as exc:
            raise utils.PlinkError(
                "GEN 格式导入失败（常见于多样染色体/多文件 GEN）",
                command=exc.command, returncode=exc.returncode, log_path=exc.log_path,
                log_tail=exc.log_tail,
                diagnosis=["GEN 输入建议用 PLINK 2 的 --gen/--sample 或先转 VCF",
                           "确认为单染色体 GEN 时可加 --oxford-single-chr 重试"]) from exc
    elif kind == "23andme":
        plink.run(["--23file", str(files["23andme"]), "--make-bed"], stage="intake",
                  name="23andme_to_bed", out_prefix=out, expect=[".bed", ".bim", ".fam"])
        ctx.logger.warning("23andMe 原始数据通常只有 1 个样本，无法做群体/关联分析，仅适合个体 QC")
    else:
        raise ValueError(f"不支持的输入类型: {kind}")
    return out


def normalize_bed(ctx, prefix: Path) -> Path:
    """Strip chr prefixes, fix missing/duplicate variant IDs and sort the fileset."""
    data_dir = ctx.stage_dir / "data"
    bim_path = prefix.with_suffix(".bim")
    bim = utils.read_bim(bim_path)
    changed = False

    # 1) strip 'chr' prefix / normalise numeric codes
    new_chr = [str(c).strip() for c in bim["CHR"]]
    normalised = []
    for code in new_chr:
        text = str(code)
        if text.lower().startswith("chr"):
            text = text[3:]
        normalised.append(text)
    if normalised != new_chr:
        bim["CHR"] = normalised
        changed = True

    # 2) missing variant IDs
    missing_ids = bim["SNP"].astype(str).isin({".", ""})
    if missing_ids.any():
        ctx.logger.info("发现 %d 个缺失位点 ID，将用 CHR:BP 生成", int(missing_ids.sum()))

    # 3) duplicate variant IDs -> suffix duplicates
    duplicated = bim["SNP"].astype(str).duplicated(keep=False)
    if duplicated.any():
        counts: dict[str, int] = {}
        new_ids = []
        for snp in bim["SNP"].astype(str):
            idx = counts.get(snp, 0)
            counts[snp] = idx + 1
            new_ids.append(snp if idx == 0 else f"{snp}_dup{idx}")
        bim["SNP"] = new_ids
        changed = True
        ctx.logger.warning("发现 %d 个重复位点 ID，已加 _dup 后缀", int(duplicated.sum()))

    if changed:
        bim.to_csv(bim_path, sep="\t", header=False, index=False, encoding="utf-8")

    out = data_dir / "base"
    args = ["--bfile", str(prefix), "--make-bed", "--keep-allele-order"]
    if missing_ids.any():
        args += ["--set-missing-var-ids", "@:#"]
    ctx.plink.run(args, stage="intake", name="normalize_bed", out_prefix=out,
                  expect=[".bed", ".bim", ".fam"])
    return out


# ---------------------------------------------------------------------------
# preflight checks
# ---------------------------------------------------------------------------
def _check(name: str, status: str, detail: str, **extra: Any) -> dict[str, Any]:
    item = {"name": name, "status": status, "detail": detail}
    item.update(extra)
    return item


def _collect_target_ids(path: Path | str | None) -> tuple[set[tuple[str, str]], dict[str, str]]:
    """Read FID/IID pairs (and IID->FID map) from a phenotype/group style table."""
    pairs: set[tuple[str, str]] = set()
    iid_to_fid: dict[str, str] = {}
    if not path:
        return pairs, iid_to_fid
    try:
        df, _ = read_table_auto(Path(path))
    except Exception:  # noqa: BLE001
        return pairs, iid_to_fid
    fid_col = _find_col(df, ("FID", "FAMILY", "FAMILYID"))
    iid_col = _find_col(df, ("IID", "ID", "INDIVIDUAL"))
    if iid_col is None:
        if df.shape[1] >= 3:
            fid_col, iid_col = df.columns[0], df.columns[1]
        elif df.shape[1] >= 2:
            fid_col, iid_col = None, df.columns[0]
        else:
            return pairs, iid_to_fid
    for idx, row in df.iterrows():
        iid = str(row[iid_col]).strip()
        fid = str(row[fid_col]).strip() if fid_col is not None else ""
        if iid:
            pairs.add((fid, iid))
            if fid:
                iid_to_fid.setdefault(iid, fid)
    return pairs, iid_to_fid


def _remap_sample_ids(base: Path, cfg: dict[str, Any], ctx) -> bool:
    """VCF sample names are often ``FID_IID`` or ``FID:IID`` (PLINK --double-id
    keeps the whole string as both FID and IID).  When phenotype/group IDs do
    not match, try to recover the original FID/IID by splitting on a delimiter."""
    pheno_path = cfg["input"]["phenotype"].get("path")
    group_path = cfg["input"]["groups"].get("path")
    pairs: set[tuple[str, str]] = set()
    iid_to_fid: dict[str, str] = {}
    for path in (pheno_path, group_path):
        file_pairs, file_map = _collect_target_ids(path)
        pairs |= file_pairs
        for iid, fid in file_map.items():
            iid_to_fid.setdefault(iid, fid)
    if not pairs:
        return False

    fam = utils.read_fam(base.with_suffix(".fam"))
    n = len(fam)
    if n == 0:
        return False
    direct = sum((str(f), str(i)) in pairs for f, i in zip(fam["FID"], fam["IID"]))
    if direct >= 0.9 * n:
        return False

    delimiters = ["_", ":", "-", ".", "|"]
    new_ids: list[tuple[str, str]] = []
    changed = 0
    for fid, iid in zip(fam["FID"].astype(str), fam["IID"].astype(str)):
        if (fid, iid) in pairs:
            new_ids.append((fid, iid))
            continue
        mapped: tuple[str, str] | None = None
        if iid in pairs and iid_to_fid.get(iid):
            mapped = (iid_to_fid[iid], iid)
        if mapped is None:
            for delimiter in delimiters:
                if delimiter not in iid:
                    continue
                left, right = iid.split(delimiter, 1)
                candidate = (left.strip(), right.strip())
                if candidate in pairs:
                    mapped = candidate
                    break
                candidate_rev = (right.strip(), left.strip())
                if candidate_rev in pairs:
                    mapped = candidate_rev
                    break
        if mapped is None and iid in iid_to_fid:
            mapped = (iid_to_fid[iid], iid)
        if mapped is None:
            new_ids.append((fid, iid))
        else:
            new_ids.append(mapped)
            if mapped != (fid, iid):
                changed += 1

    if changed == 0 or changed < 0.5 * n:
        return False
    fam["FID"] = [item[0] for item in new_ids]
    fam["IID"] = [item[1] for item in new_ids]
    fam.to_csv(base.with_suffix(".fam"), sep="\t", header=False, index=False, encoding="utf-8")
    ctx.logger.warning("已把 %d/%d 个样本 ID 按分隔符重映射为表型/分组 ID（请核对 .fam）", changed, n)
    return True


def _log_count(text: str, pattern: str) -> int:
    import re
    match = re.search(pattern, text)
    return int(match.group(1)) if match else 0


def _hash_or_size(path: Path) -> dict[str, Any]:
    size = path.stat().st_size
    digest = utils.sha256_file(path) if size <= MAX_HASH_BYTES else None
    return {"path": str(path), "size": size, "sha256": digest}


def run(ctx) -> dict[str, Any]:
    cfg = ctx.config
    stage = ctx.stage_dir
    data_dir = stage / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    plink: utils.Plink = ctx.plink
    checks: list[dict[str, Any]] = []

    # ---- discover + convert --------------------------------------------
    discovery = discover_inputs(cfg["input"]["genotype"]["path"] or cfg["input"]["genotype"]["prefix"])
    ctx.logger.info("识别到输入类型: %s", discovery["type"])
    if len(discovery.get("all_candidates", [])) > 1:
        checks.append(_check("input_candidates", "warn",
                             f"发现 {len(discovery['all_candidates'])} 个候选数据集，已选用 {discovery['type']}；"
                             "如需指定请用 --input 指向具体前缀或目录"))
    raw_prefix = convert_to_bed(ctx, discovery)
    base = normalize_bed(ctx, raw_prefix)
    if _remap_sample_ids(base, cfg, ctx):
        checks.append(_check("sample_id_remap", "warn",
                             "VCF/文本样本名与表型/分组 ID 不一致，已按分隔符自动拆分 FID/IID；"
                             "请核对 00_intake/manifest.json 与 .fam 前几行"))

    fam = utils.read_fam(base.with_suffix(".fam"))
    bim = utils.read_bim(base.with_suffix(".bim"))
    n_samples, n_variants = len(fam), len(bim)

    # ---- species inference ---------------------------------------------
    species = infer_species(bim["CHR"].astype(str).tolist(), cfg.get("species", {}).get("name", "auto"))
    cfg["species"] = {**species,
                      "name": species["name"], "label": species["label"],
                      "autosome_count": species["autosome_count"],
                      "sex_system": species["sex_system"],
                      "chr_set": species["chr_set"],
                      "allow_extra_chr": species["allow_extra_chr"]}
    ctx.logger.info("推断物种: %s (%s)", species["label"], species["evidence"])

    # ---- basic genotype checks -----------------------------------------
    has_bed = base.with_suffix(".bed").exists()
    checks.append(_check("genotype_files", "pass" if has_bed else "fail",
                         f"PLINK 前缀 {base}；样本 {n_samples}，位点 {n_variants}"))
    if n_samples == 0 or n_variants == 0:
        checks.append(_check("genotype_nonempty", "fail", "样本数或位点数为 0，无法继续"))
    dup_samples = int(fam["IID"].duplicated().sum())
    checks.append(_check("sample_ids_unique", "pass" if dup_samples == 0 else "fail",
                         f"重复 IID 数量 {dup_samples}"))
    dup_variants = int(bim["SNP"].duplicated().sum())
    checks.append(_check("variant_ids_unique", "pass" if dup_variants == 0 else "warn",
                         f"重复位点 ID 数量 {dup_variants}（normalize 已自动加后缀）"))
    bad_alleles = int((~bim["A1"].astype(str).str.fullmatch(r"[ACGTNacgtn0]")).sum() +
                      (~bim["A2"].astype(str).str.fullmatch(r"[ACGTNacgtn0]")).sum())
    checks.append(_check("allele_codes", "pass" if bad_alleles == 0 else "warn",
                         f"非 A/C/G/T/N 等位基因记录 {bad_alleles} 条"))
    ambiguous = int((bim["A1"].str.upper().isin({"A", "T"}) & bim["A2"].str.upper().isin({"A", "T"})).sum() +
                    (bim["A1"].str.upper().isin({"C", "G"}) & bim["A2"].str.upper().isin({"C", "G"})).sum())
    checks.append(_check("ambiguous_strand_snps", "warn" if ambiguous > 0 else "pass",
                         f"AT/GC 模糊位点 {ambiguous} 个（链方向无法从数据判断，跨平台合并时需注意）"))
    zero_bp = int((bim["BP"].astype(float) == 0).sum())
    checks.append(_check("variant_positions", "warn" if zero_bp > 0 else "pass",
                         f"BP=0 的位点 {zero_bp} 个"))
    sex_codes = fam["SEX"].value_counts().to_dict()
    checks.append(_check("sex_codes", "pass" if any(int(k) > 0 for k in sex_codes) else "warn",
                         f"性别编码分布 {sex_codes}（0=未知/缺失）；未知性别样本过多时性别检查关闭"))
    checks.append(_check("chromosome_codes", "pass",
                         f"染色体 {sorted(bim['CHR'].unique().tolist(), key=str)[:15]}；"
                         f"推断 {species['label']}，chr-set={species['chr_set']}"))

    # ---- phenotype ------------------------------------------------------
    pheno_info: dict[str, Any] = {"present": False}
    pheno_table = None
    pheno_path = cfg["input"]["phenotype"].get("path")
    if pheno_path:
        try:
            pheno_table, pheno_info = load_phenotype(Path(pheno_path), fam,
                                                     cfg["input"]["phenotype"]["missing"])
            pheno_info["present"] = True
            pheno_table.to_csv(data_dir / "pheno_clean.txt", sep="\t", header=False,
                               index=False, encoding="utf-8")
            match_rate = utils.safe_div(pheno_info["n_matched"], len(fam))
            status = "pass" if match_rate >= 0.95 else ("warn" if match_rate >= 0.5 else "fail")
            checks.append(_check("phenotype_file", status,
                                 f"匹配 {pheno_info['n_matched']}/{len(fam)} 个样本"
                                 f"（{pheno_info['n_unmatched']} 行未匹配）；类型={pheno_info['type']}，"
                                 f"编码={pheno_info['coding']}",
                                 n_matched=pheno_info["n_matched"], match_rate=match_rate))
            if pheno_info["type"] == "binary":
                checks.append(_check("phenotype_coding", "pass",
                                     f"病例 {pheno_info.get('n_cases', 0)}，对照 {pheno_info.get('n_controls', 0)}"))
            elif pheno_info["type"] == "quantitative":
                series = pheno_table["PHENO"].replace(-9, np.nan).dropna()
                checks.append(_check("phenotype_coding", "pass",
                                     f"数量性状 n={len(series)}，均值={series.mean():.3f}，"
                                     f"标准差={series.std():.3f}，缺失={len(pheno_table) - len(series)}"))
            else:
                checks.append(_check("phenotype_coding", "warn", "表型全部缺失或无法识别类型"))
        except Exception as exc:  # noqa: BLE001 - surface as a preflight failure
            checks.append(_check("phenotype_file", "fail", f"读取失败: {exc}"))
    else:
        checks.append(_check("phenotype_file", "pass", "未提供表型文件，GWAS 阶段将跳过"))

    # ---- covariates -----------------------------------------------------
    covar_info: dict[str, Any] = {"present": False}
    covar_table = None
    covar_path = cfg["input"]["covariates"].get("path")
    if covar_path:
        try:
            covar_table, covar_info = load_covariates(Path(covar_path), fam)
            covar_info["present"] = True
            covar_table.to_csv(data_dir / "covar_clean.txt", sep="\t", header=False,
                               index=False, encoding="utf-8")
            match_rate = utils.safe_div(covar_info["n_matched"], len(fam))
            checks.append(_check("covariate_file", "pass" if match_rate >= 0.95 else "warn",
                                 f"匹配 {covar_info['n_matched']}/{len(fam)}；数值列 "
                                 f"{covar_info['columns']}；分类列自动哑变量编码 {covar_info['categorical']}"))
        except Exception as exc:  # noqa: BLE001
            checks.append(_check("covariate_file", "fail", f"读取失败: {exc}"))
    else:
        checks.append(_check("covariate_file", "pass", "未提供协变量文件，将只使用 PCA 主成分"))

    # ---- groups ---------------------------------------------------------
    group_info: dict[str, Any] = {"present": False}
    group_path = cfg["input"]["groups"].get("path")
    if group_path:
        try:
            group_table, group_info = load_groups(Path(group_path), fam)
            group_info["present"] = True
            group_table.to_csv(data_dir / "group_clean.txt", sep="\t", header=False,
                               index=False, encoding="utf-8")
            within = group_table.loc[group_table["GROUP"] != "", ["FID", "IID", "GROUP"]]
            within.to_csv(data_dir / "within.txt", sep="\t", header=False, index=False,
                          encoding="utf-8")
            match_rate = utils.safe_div(group_info["n_matched"], len(fam))
            checks.append(_check("group_file", "pass" if match_rate >= 0.9 else "warn",
                                 f"匹配 {group_info['n_matched']}/{len(fam)}；分组数 "
                                 f"{group_info['n_levels']}：{group_info['levels'][:8]}"))
        except Exception as exc:  # noqa: BLE001
            checks.append(_check("group_file", "fail", f"读取失败: {exc}"))
    else:
        checks.append(_check("group_file", "pass", "未提供分组文件，FST 等分组分析将跳过"))

    # ---- manifest -------------------------------------------------------
    files: list[dict[str, Any]] = []
    for key, file in discovery["files"].items():
        try:
            files.append({"role": key, **_hash_or_size(Path(file))})
        except OSError as exc:
            checks.append(_check("file_readable", "fail", f"{file}: {exc}"))
    for extras in (pheno_path, covar_path, group_path):
        if extras and Path(extras).exists():
            files.append({"role": "aux", **_hash_or_size(Path(extras))})
    manifest = {
        "discovered_type": discovery["type"],
        "files": files,
        "n_samples": int(n_samples),
        "n_variants": int(n_variants),
        "chromosomes": sorted(set(bim["CHR"].astype(str)))[:60],
        "species": species,
        "plink": {"exe": plink.exe, "version": plink.version},
        "python": {"executable": utils.sys.executable},
    }
    utils.write_json(stage / "manifest.json", manifest)
    utils.write_json(stage / "preflight.json",
                     {"checks": checks,
                      "overall": _overall_status(checks)})
    (stage / "base_prefix.txt").write_text(str(base), encoding=utils.ENCODING)

    # ---- config bookkeeping --------------------------------------------
    cfg["species"].update({k: species[k] for k in
                           ("name", "label", "autosome_count", "sex_system", "chr_set")})
    cfg["input"]["genotype"].update({
        "type": discovery["type"], "prefix": str(base), "converted": discovery["type"] != "bed",
    })
    cfg.setdefault("_derived", {}).update({
        "base_prefix": str(base),
        "n_samples": int(n_samples), "n_variants": int(n_variants),
        "chromosomes": sorted(set(bim["CHR"].astype(str))),
        "chr_flags": chr_flags(species),
        "phenotype_type": pheno_info.get("type", "none"),
        "phenotype_present": bool(pheno_info.get("present")),
        "covariate_columns": covar_info.get("columns", []),
        "covariate_present": bool(covar_info.get("present")),
        "group_present": bool(group_info.get("present")),
        "group_levels": group_info.get("levels", []),
    })

    status = _overall_status(checks)
    metrics = {
        "stage": "intake",
        "status": status,
        "genotype_type": discovery["type"],
        "n_samples": int(n_samples),
        "n_variants": int(n_variants),
        "chromosomes": sorted(set(bim["CHR"].astype(str)))[:60],
        "species_inferred": species,
        "files": files,
        "preflight": {"status": status, "checks": checks},
        "phenotype": pheno_info,
        "covariates": covar_info,
        "groups": group_info,
        "tables": {},
        "figures": [],
        "summary": _summary(n_samples, n_variants, species, pheno_info, covar_info,
                            group_info, checks),
        "warnings": [c["detail"] for c in checks if c["status"] == "warn"],
        "errors": [c["detail"] for c in checks if c["status"] == "fail"],
    }
    utils.write_json(stage / "metrics.json", metrics)
    ctx.logger.info("intake 完成：%d 样本 × %d 位点，状态 %s", n_samples, n_variants, status)
    return metrics


def _overall_status(checks: list[dict[str, Any]]) -> str:
    if any(c["status"] == "fail" for c in checks):
        return "fail"
    if any(c["status"] == "warn" for c in checks):
        return "warn"
    return "pass"


def _summary(n_samples: int, n_variants: int, species: dict[str, Any],
             pheno: dict[str, Any], covar: dict[str, Any], groups: dict[str, Any],
             checks: list[dict[str, Any]]) -> str:
    lines = [
        f"数据接入完成：{n_samples} 个样本 × {n_variants} 个位点。",
        f"物种推断：{species['label']}（常染色体数 {species['autosome_count']}，"
        f"性染色体系统 {species['sex_system']}，依据：{species['evidence']}）。",
    ]
    if pheno.get("present"):
        lines.append(f"表型：{pheno.get('type')}，匹配 {pheno.get('n_matched')}/{n_samples}，"
                     f"编码 {pheno.get('coding')}。")
    else:
        lines.append("表型：未提供，后续仅做群体遗传分析。")
    if covar.get("present"):
        lines.append(f"协变量：{len(covar.get('columns', []))} 列（含哑变量）。")
    if groups.get("present"):
        lines.append(f"分组：{groups.get('n_levels')} 组。")
    failed = [c["name"] for c in checks if c["status"] == "fail"]
    warned = [c["name"] for c in checks if c["status"] == "warn"]
    if failed:
        lines.append("阻断性问题：" + "、".join(failed) + "，请先修复再继续。")
    elif warned:
        lines.append("需注意：" + "、".join(warned) + "（详见 00_intake/preflight.json）。")
    else:
        lines.append("预检未发现阻断性问题。")
    return "\n".join(lines)
