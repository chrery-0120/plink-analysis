#!/usr/bin/env python
"""PLINK 全流程分析 —— 命令行入口.

用法示例::

    python plink_pipeline.py auto --input D:/data/pig --pheno pheno.txt --covar covar.txt \\
        --group group.txt --species auto --out-root D:/plink_runs

    python plink_pipeline.py run --run-dir D:/plink_runs/20261023_2200_pig
    python plink_pipeline.py status --run-dir ...
    python plink_pipeline.py report --run-dir ...

每个阶段结束后会打印中文摘要；完整过程写入 ``logs/pipeline.log``。
"""
from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from pipeline import __version__ as PIPELINE_VERSION
from pipeline import config as config_mod  # noqa: E402
from pipeline import utils  # noqa: E402

STAGES = ["intake", "qc", "popgen", "gwas", "calibration"]
STAGE_MODULES = {
    "intake": "pipeline.intake",
    "qc": "pipeline.stage_qc",
    "popgen": "pipeline.stage_popgen",
    "gwas": "pipeline.stage_gwas",
    "calibration": "pipeline.calibration_runner",
}


# ---------------------------------------------------------------------------
# context
# ---------------------------------------------------------------------------
@dataclass
class PipelineContext:
    run_dir: Path
    config: dict[str, Any]
    plink: utils.Plink
    logger: Any
    state: dict[str, Any]
    stage: str = "intake"

    @property
    def stage_dir(self) -> Path:
        return utils.stage_dir(self.run_dir, self.stage)


def build_context(run_dir: Path, cfg: dict[str, Any], logger, plink_exe: str | None = None) -> PipelineContext:
    runtime = cfg.get("runtime", {})
    plink_path = _find_plink(plink_exe)
    plink = utils.Plink(plink_path, run_dir,
                        threads=int(runtime.get("threads", 4)),
                        memory_mb=int(runtime.get("memory_mb", 4000)),
                        logger=logger, timeout_s=int(runtime.get("timeout_s", 7200)))
    state = utils.load_json(run_dir / "state.json", default={"run_id": cfg.get("run_id"),
                                                             "stages": {}})
    return PipelineContext(run_dir=run_dir, config=cfg, plink=plink, logger=logger, state=state)


def _workspace_dir() -> Path:
    """<workspace> when the skill lives in <workspace>/.agents/skills/plink-analysis."""
    return Path(__file__).resolve().parents[4]


def _plink_candidates() -> list[Path]:
    import shutil
    names = ("plink.exe",) if os.name == "nt" else ("plink",)
    candidates: list[Path] = []
    for env_name in ("PLINK_EXE", "PLINK"):
        env = os.environ.get(env_name)
        if env:
            candidates.append(Path(env))
    for base in (_workspace_dir(), Path.cwd()):
        for name in names:
            candidates.append(base / "tools" / "plink" / name)
            candidates.append(base / name)
    for name in names:
        found = shutil.which(name)
        if found:
            candidates.append(Path(found))
    return candidates


def _find_plink(explicit: str | None) -> str:
    candidates = [Path(explicit)] if explicit else []
    candidates.extend(_plink_candidates())
    for path in candidates:
        try:
            resolved = path.expanduser().resolve()
        except OSError:
            continue
        if resolved.exists():
            return str(resolved)
    raise FileNotFoundError(
        "找不到 PLINK 可执行文件。请运行 `python scripts/install_plink.py` 安装，"
        "或用 --plink-exe 指定 plink.exe 路径，或设置环境变量 PLINK_EXE。")


# ---------------------------------------------------------------------------
# stage runner
# ---------------------------------------------------------------------------
def run_stage(name: str, ctx: PipelineContext, *, force: bool = False) -> dict[str, Any]:
    ctx.stage = name
    stage_start = time.time()
    entry = ctx.state.setdefault("stages", {}).setdefault(name, {})
    entry.update({"status": "running", "started": utils.now_iso(), "error": None})
    utils.write_json(ctx.run_dir / "state.json", ctx.state)
    ctx.logger.info("=== 开始阶段 %s ===", name)

    metrics: dict[str, Any] = {}
    status = "error"
    error_message = None
    try:
        module = importlib.import_module(STAGE_MODULES[name])
        metrics = module.run(ctx)
        status = str(metrics.get("status", "pass"))
        metrics.setdefault("overall", status)
        # DESIGN section 8 names the aggregate field "overall"; keep both names
        try:
            utils.write_json(ctx.stage_dir / "metrics.json", metrics)
        except Exception:  # noqa: BLE001 - never fail a stage just for the alias
            pass
    except Exception as exc:  # noqa: BLE001 - record, report and keep going
        error_message = f"{type(exc).__name__}: {exc}"
        status = "error"
        metrics = {"stage": name, "status": "error", "summary": f"阶段执行失败: {error_message}",
                   "warnings": [], "errors": [error_message], "tables": {}, "figures": []}
        ctx.logger.error("阶段 %s 失败: %s", name, error_message)
        ctx.logger.debug(traceback.format_exc())
        utils.write_json(ctx.stage_dir / "metrics.json", metrics)

    elapsed = time.time() - stage_start
    report_path = _write_stage_report(ctx, name, metrics)
    entry.update({
        "status": status,
        "finished": utils.now_iso(),
        "elapsed_s": round(elapsed, 3),
        "metrics_path": utils.rel_to_run(ctx.run_dir, ctx.stage_dir / "metrics.json"),
        "report_path": utils.rel_to_run(ctx.run_dir, report_path) if report_path else None,
        "artifacts": _artifacts(name, metrics, ctx),
        "error": error_message,
        "summary": metrics.get("summary", ""),
    })
    ctx.state["updated_at"] = utils.now_iso()
    utils.write_json(ctx.run_dir / "state.json", ctx.state)
    # persist config learned during the stage (species, chr flags, covariate columns...)
    utils.write_json(ctx.run_dir / "run.json", ctx.config)

    _print_stage(ctx, name, metrics, status, elapsed)
    return metrics


def _write_stage_report(ctx: PipelineContext, name: str, metrics: dict[str, Any]) -> Path | None:
    try:
        report_mod = importlib.import_module("pipeline.report")
        return report_mod.stage_report(ctx.run_dir, name, metrics, ctx.config)
    except Exception as exc:  # noqa: BLE001 - fallback keeps the pipeline usable
        ctx.logger.warning("阶段报告生成失败（使用简易报告）: %s", exc)
        return _fallback_report(ctx, name, metrics)


def _fallback_report(ctx: PipelineContext, name: str, metrics: dict[str, Any]) -> Path:
    stage_name_cn = {"intake": "数据接入", "qc": "质量控制", "popgen": "群体结构",
                     "gwas": "关联分析", "calibration": "校准与校验"}.get(name, name)
    lines = [f"# {stage_name_cn}报告（简易）", "",
             f"- 状态: **{metrics.get('status', 'unknown')}**",
             f"- 摘要: {metrics.get('summary', '')}", ""]
    if metrics.get("tables"):
        lines.append("| 表格 | 路径 |")
        lines.append("|---|---|")
        for key, value in metrics["tables"].items():
            if value:
                lines.append(f"| {key} | {value} |")
        lines.append("")
    if metrics.get("figures"):
        lines.append("## 图")
        for figure in metrics["figures"]:
            lines.append(f"![]({figure})")
        lines.append("")
    if metrics.get("warnings"):
        lines.append("## 警告")
        lines.extend(f"- {item}" for item in metrics["warnings"])
    if metrics.get("errors"):
        lines.append("## 错误")
        lines.extend(f"- {item}" for item in metrics["errors"])
    path = ctx.stage_dir / "report.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _artifacts(name: str, metrics: dict[str, Any], ctx: PipelineContext) -> dict[str, str]:
    run_dir = ctx.run_dir
    if name == "intake":
        return {"base_prefix": str(run_dir / "00_intake" / "base_prefix.txt"),
                "pheno": str(run_dir / "00_intake" / "data" / "pheno_clean.txt"),
                "covar": str(run_dir / "00_intake" / "data" / "covar_clean.txt"),
                "group": str(run_dir / "00_intake" / "data" / "group_clean.txt")}
    if name == "qc":
        return {"bfile": str(run_dir / "01_qc" / "data" / "clean"),
                "pruned": str(run_dir / "01_qc" / "data" / "clean_pruned.prune.in"),
                "eigenvec": str(run_dir / "01_qc" / "tables" / "clean_pca.eigenvec"),
                "eigenval": str(run_dir / "01_qc" / "tables" / "clean_pca.eigenval"),
                "pca_covar": str(run_dir / "01_qc" / "data" / "pca_covar.txt")}
    if name == "popgen":
        return {"eigenvec": str(run_dir / "02_popgen" / "tables" / "pca.eigenvec"),
                "eigenval": str(run_dir / "02_popgen" / "tables" / "pca.eigenval"),
                "genome": str(run_dir / "02_popgen" / "tables" / "genome.genome"),
                "rel": str(run_dir / "02_popgen" / "tables" / "kinship.rel")}
    if name == "gwas":
        return {"assoc": str(metrics.get("assoc_file") or ""),
                "adjusted": str(metrics.get("adjusted_file") or ""),
                "covar": str(metrics.get("covar_file") or ""),
                "results": str(run_dir / "03_gwas" / "tables" / "gwas_results.csv")}
    if name == "calibration":
        return {"validation": str(run_dir / "04_calibration" / "validation.json"),
                "metrics": str(run_dir / "04_calibration" / "metrics.json")}
    return {}


def _print_stage(ctx: PipelineContext, name: str, metrics: dict[str, Any], status: str,
                 elapsed: float) -> None:
    stage_cn = {"intake": "数据接入", "qc": "质量控制", "popgen": "群体结构",
                "gwas": "关联分析", "calibration": "校准与校验"}.get(name, name)
    bar = "=" * 24
    print(f"\n{bar} 阶段 {name} · {stage_cn} | {status.upper()} | {utils.fmt_elapsed(elapsed)} {bar}")
    summary = str(metrics.get("summary") or "").strip()
    if summary:
        print(summary)
    for table_key, table_path in (metrics.get("tables") or {}).items():
        if table_path:
            print(f"  表格: {table_key} -> {table_path}")
    for figure in metrics.get("figures") or []:
        print(f"  图: {figure}")
    if metrics.get("errors"):
        print("  错误:")
        for item in metrics["errors"][:10]:
            print(f"    - {item}")
    print(f"  阶段目录: {utils.rel_to_run(ctx.run_dir, ctx.stage_dir)}")


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------
def cmd_auto(args: argparse.Namespace) -> int:
    input_path = Path(args.input).expanduser().resolve()
    if not (input_path.exists() or Path(str(input_path) + ".bed").exists()
            or Path(str(input_path) + ".vcf").exists()):
        print(f"ERROR: 输入路径不存在: {input_path}")
        return 2
    cfg = config_mod.default_config()
    cfg["run_id"] = args.run_id or (utils.datetime.now().strftime("%Y%m%d_%H%M") + "_" +
                                    input_path.stem[:20])
    cfg["created_at"] = utils.now_iso()
    cfg["input"]["genotype"] = {"type": "auto", "prefix": None, "path": str(input_path)}
    auto_files = _autodetect_aux_files(input_path)
    pheno_path = args.pheno or auto_files.get("pheno")
    covar_path = args.covar or auto_files.get("covar")
    group_path = args.group or auto_files.get("group")
    for label, value in (("表型", pheno_path), ("协变量", covar_path), ("分组", group_path)):
        if value and not args.__dict__.get({"表型": "pheno", "协变量": "covar", "分组": "group"}[label]):
            print(f"自动识别{label}文件: {value}")
    if pheno_path:
        cfg["input"]["phenotype"]["path"] = str(Path(pheno_path).expanduser().resolve())
    if covar_path:
        cfg["input"]["covariates"]["path"] = str(Path(covar_path).expanduser().resolve())
    if group_path:
        cfg["input"]["groups"]["path"] = str(Path(group_path).expanduser().resolve())
    if args.species:
        cfg["species"]["name"] = args.species
    if args.title:
        cfg["report"]["title"] = args.title
    if args.threads:
        cfg["runtime"]["threads"] = int(args.threads)
    if args.memory:
        cfg["runtime"]["memory_mb"] = int(args.memory)
    if args.no_calibration:
        cfg["calibration"]["enabled"] = False

    errors, warnings = config_mod.validate_config(cfg)
    if errors:
        print("配置错误:\n  - " + "\n  - ".join(errors))
        return 2
    for warning in warnings:
        print(f"配置提示: {warning}")

    out_root = Path(args.out_root).expanduser().resolve() if args.out_root else Path.cwd() / "plink_runs"
    run_dir = out_root / cfg["run_id"]
    if run_dir.exists() and any(run_dir.iterdir()):
        if not args.force:
            print(f"ERROR: 运行目录已存在且非空: {run_dir}（用 --force 覆盖或换 --run-id）")
            return 2

    logger = utils.setup_logger(run_dir / "logs" / "pipeline.log")
    utils.ensure_utf8_stdout()
    cfg["_paths"] = {"run_dir": str(run_dir), "out_root": str(out_root)}
    utils.write_json(run_dir / "run.json", cfg)
    utils.write_json(run_dir / "state.json", {"run_id": cfg["run_id"], "stages": {},
                                              "created_at": cfg["created_at"]})
    print(f"运行目录: {run_dir}")
    print(f"输入: {input_path}")

    ctx = build_context(run_dir, cfg, logger, args.plink_exe)
    stages = _parse_stages(args.stages)
    if args.dry_run:
        print("DRY-RUN: 将执行阶段 " + ", ".join(stages))
        return 0
    return _run_stages(ctx, stages, force=args.force)


def cmd_run(args: argparse.Namespace) -> int:
    run_dir = Path(args.run_dir).expanduser().resolve()
    run_json = run_dir / "run.json"
    if not run_json.exists():
        print(f"ERROR: 找不到 {run_json}")
        return 2
    cfg = utils.load_json(run_json)
    logger = utils.setup_logger(run_dir / "logs" / "pipeline.log")
    utils.ensure_utf8_stdout()
    ctx = build_context(run_dir, cfg, logger, args.plink_exe)
    stages = _parse_stages(args.stage if args.stage else "all")
    return _run_stages(ctx, stages, force=args.force)


def _run_stages(ctx: PipelineContext, stages: list[str], *, force: bool) -> int:
    statuses: dict[str, str] = {}
    for stage in stages:
        metrics = run_stage(stage, ctx, force=force)
        status = str(metrics.get("status", "error"))
        statuses[stage] = status
        if status in {"fail", "error"} and not force:
            print(f"\n阶段 {stage} 状态为 {status.upper()}，已停止。修复后可用 "
                  f"`run --run-dir {ctx.run_dir} --stage {stage}` 续跑，或加 --force 强制继续。")
            break
    final_path = _write_final_report(ctx)
    _print_overall(ctx, statuses, final_path)
    if any(s in {"fail", "error"} for s in statuses.values()):
        return 2
    if any(s == "warn" for s in statuses.values()):
        return 1
    return 0


def _write_final_report(ctx: PipelineContext) -> Path | None:
    try:
        report_mod = importlib.import_module("pipeline.report")
        result = report_mod.final_report(ctx.run_dir, ctx.config, ctx.state)
        md = result.get("md") if isinstance(result, dict) else result
        return Path(md) if md else None
    except Exception as exc:  # noqa: BLE001
        ctx.logger.warning("汇总报告生成失败（使用简易汇总）: %s", exc)
        reports_dir = ctx.run_dir / "05_reports"
        reports_dir.mkdir(parents=True, exist_ok=True)
        lines = ["# PLINK 全流程分析报告（简易汇总）", "",
                 f"- 运行: {ctx.config.get('run_id')}",
                 f"- 生成时间: {utils.now_iso()}", ""]
        for stage, entry in ctx.state.get("stages", {}).items():
            lines.append(f"## {stage} — {entry.get('status')}")
            lines.append("")
            summary = entry.get("summary") or ""
            if summary:
                lines.append(summary)
            lines.append("")
        path = reports_dir / "final_report.md"
        path.write_text("\n".join(lines), encoding="utf-8")
        return path


def _print_overall(ctx: PipelineContext, statuses: dict[str, str], final_path: Path | None) -> None:
    print("\n" + "=" * 30 + " 流程结束 " + "=" * 30)
    for stage, status in statuses.items():
        print(f"  {stage:<12} {status.upper()}")
    if final_path:
        print(f"  汇总报告: {final_path}")
    print("  阶段产物目录: " + str(ctx.run_dir))
    print("SUMMARY_JSON: " + json.dumps({"run_id": ctx.config.get("run_id"),
                                         "stages": statuses,
                                         "final_report": str(final_path) if final_path else None},
                                        ensure_ascii=False))


def cmd_status(args: argparse.Namespace) -> int:
    run_dir = Path(args.run_dir).expanduser().resolve()
    state = utils.load_json(run_dir / "state.json")
    if not state:
        print(f"ERROR: 找不到 {run_dir / 'state.json'}")
        return 2
    utils.ensure_utf8_stdout()
    print(f"运行: {state.get('run_id')}  更新: {state.get('updated_at', '-')}")
    for stage in STAGES:
        entry = state.get("stages", {}).get(stage)
        if not entry:
            print(f"  {stage:<12} {'-':<8}")
            continue
        elapsed = entry.get("elapsed_s")
        print(f"  {stage:<12} {entry.get('status', '-'):<8} {entry.get('finished', '-'):<20} "
              f"{str(elapsed) + 's' if elapsed else '':<8} {entry.get('report_path') or ''}")
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    run_dir = Path(args.run_dir).expanduser().resolve()
    run_json = run_dir / "run.json"
    if not run_json.exists():
        print(f"ERROR: 找不到 {run_json}")
        return 2
    cfg = utils.load_json(run_json)
    logger = utils.setup_logger(run_dir / "logs" / "pipeline.log")
    utils.ensure_utf8_stdout()
    ctx = build_context(run_dir, cfg, logger, args.plink_exe)
    path = _write_final_report(ctx)
    print(f"汇总报告: {path}")
    return 0


def cmd_selftest(args: argparse.Namespace) -> int:
    script = SCRIPT_DIR / "selftest.py"
    if not script.exists():
        print(f"ERROR: 找不到自检脚本 {script}")
        return 2
    import subprocess
    return subprocess.call([sys.executable, str(script), *args.rest])


# ---------------------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------------------
def _autodetect_aux_files(input_path: Path) -> dict[str, str]:
    """Find common phenotype/covariate/group files sitting next to the genotype data."""
    folder = input_path if input_path.is_dir() else input_path.parent
    if not folder.is_dir():
        return {}
    patterns = {
        "pheno": ("pheno.txt", "phenotype.txt", "pheno.tsv", "phenotype.tsv",
                  "phenotypes.txt", "trait.txt", "traits.txt"),
        "covar": ("covar.txt", "covariates.txt", "covar.tsv", "covariates.tsv",
                  "cov.txt", "covariates.csv"),
        "group": ("group.txt", "groups.txt", "breed.txt", "population.txt",
                  "pop.txt", "cluster.txt"),
    }
    found: dict[str, str] = {}
    available: dict[str, Path] = {p.name.lower(): p for p in folder.iterdir() if p.is_file()}
    for key, names in patterns.items():
        for name in names:
            if name in available:
                found[key] = str(available[name])
                break
    return found


def _parse_stages(value: str | None) -> list[str]:
    if not value or value == "all":
        return list(STAGES)
    out = [item.strip() for item in value.split(",") if item.strip()]
    unknown = [item for item in out if item not in STAGES]
    if unknown:
        raise SystemExit(f"未知阶段: {unknown}（可选 {STAGES}）")
    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="PLINK 全流程分析（科研级）")
    parser.add_argument("--version", action="version",
                        version=f"plink-analysis {PIPELINE_VERSION}")
    parser.add_argument("--plink-exe", default=None, help="plink 可执行文件路径")
    sub = parser.add_subparsers(dest="command", required=True)

    auto = sub.add_parser("auto", help="新建运行并执行完整流程")
    auto.add_argument("--input", required=True, help="基因型数据：目录/前缀/单文件")
    auto.add_argument("--pheno", default=None)
    auto.add_argument("--covar", default=None)
    auto.add_argument("--group", default=None)
    auto.add_argument("--species", default="auto",
                      help="auto|pig|cattle|sheep|goat|chicken|horse|dog|human|other")
    auto.add_argument("--run-id", default=None)
    auto.add_argument("--out-root", default=None)
    auto.add_argument("--threads", type=int, default=None)
    auto.add_argument("--memory", type=int, default=None, help="PLINK 内存上限 (MB)")
    auto.add_argument("--stages", default="all", help="逗号分隔，如 intake,qc,gwas")
    auto.add_argument("--title", default=None)
    auto.add_argument("--no-calibration", action="store_true")
    auto.add_argument("--dry-run", action="store_true")
    auto.add_argument("--force", action="store_true")
    auto.set_defaults(func=cmd_auto)

    run = sub.add_parser("run", help="在已有运行目录上续跑/重跑阶段")
    run.add_argument("--run-dir", required=True)
    run.add_argument("--stage", default="all")
    run.add_argument("--force", action="store_true")
    run.set_defaults(func=cmd_run)

    status = sub.add_parser("status", help="查看阶段状态")
    status.add_argument("--run-dir", required=True)
    status.set_defaults(func=cmd_status)

    report = sub.add_parser("report", help="重新生成汇总报告")
    report.add_argument("--run-dir", required=True)
    report.set_defaults(func=cmd_report)

    selftest = sub.add_parser("selftest", help="运行端到端自检")
    selftest.add_argument("rest", nargs=argparse.REMAINDER)
    selftest.set_defaults(func=cmd_selftest)
    return parser


def main(argv: list[str] | None = None) -> int:
    utils.ensure_utf8_stdout()
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
