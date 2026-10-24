"""Stage 4: simulation-based calibration + independent validation.

This runner is intentionally thin: the heavy lifting lives in
``pipeline/simulate.py`` (null / known-QTL calibration) and
``pipeline/validate.py`` (independent recomputation checks).  It exists so the
main CLI does not fail when an optional module is unavailable, and so the two
reports are merged into one stage metrics file.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from . import utils


def run(ctx) -> dict[str, Any]:
    cfg = ctx.config
    stage = ctx.stage_dir
    run_dir = ctx.run_dir
    warnings: list[str] = []
    errors: list[str] = []

    calibration_enabled = bool(cfg.get("calibration", {}).get("enabled", True))
    phenotype_present = bool(cfg.get("_derived", {}).get("phenotype_present"))
    simulation_metrics: dict[str, Any] | None = None

    if calibration_enabled and phenotype_present:
        try:
            from .simulate import run_calibration
            simulation_metrics = run_calibration(ctx)
        except ImportError as exc:
            warnings.append(f"仿真校准模块不可用，已跳过校准: {exc}")
        except Exception as exc:  # noqa: BLE001 - calibration failure must not abort reports
            warnings.append(f"仿真校准失败，已跳过: {type(exc).__name__}: {exc}")
    elif not phenotype_present:
        warnings.append("未提供表型，无法做已知 QTL 校准；仅执行独立重算校验")
    elif not calibration_enabled:
        warnings.append("配置关闭了仿真校准（calibration.enabled=false）")

    validation: dict[str, Any] = {}
    try:
        from .validate import run_validation
        validation = run_validation(ctx)
    except ImportError as exc:
        warnings.append(f"独立校验模块不可用，已跳过: {exc}")
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"独立校验失败: {type(exc).__name__}: {exc}")
        errors.append(str(exc))

    validation_summary = {
        "n_checks": len(validation.get("checks", [])),
        "n_pass": sum(1 for c in validation.get("checks", []) if c.get("status") == "pass"),
        "n_warn": sum(1 for c in validation.get("checks", []) if c.get("status") == "warn"),
        "n_fail": sum(1 for c in validation.get("checks", []) if c.get("status") == "fail"),
    }
    validation_status = validation.get("overall") or (
        "fail" if validation_summary["n_fail"] else
        ("warn" if validation_summary["n_warn"] else "pass"))

    if simulation_metrics is None:
        metrics: dict[str, Any] = {
            "stage": "calibration",
            "status": "warn" if phenotype_present or validation_summary["n_checks"] == 0 else "pass",
            "summary": "本次运行未完成仿真校准。" + ("独立校验已执行。" if validation_summary["n_checks"] else ""),
            "null": None, "power": None, "effect_bias": None, "concordance": None,
        }
    else:
        metrics = dict(simulation_metrics)

    metrics["stage"] = "calibration"
    metrics["validation"] = validation_summary
    metrics.setdefault("tables", {})
    if (stage / "tables").exists():
        for path in sorted((stage / "tables").glob("*.csv")):
            metrics["tables"].setdefault(path.stem, utils.rel_to_run(run_dir, path))
    metrics.setdefault("figures", [])
    figures_dir = stage / "plots"
    if figures_dir.exists():
        for path in sorted(figures_dir.glob("*.png")):
            if utils.rel_to_run(run_dir, path) not in metrics["figures"]:
                metrics["figures"].append(utils.rel_to_run(run_dir, path))

    statuses = [str(metrics.get("status", "pass")), str(validation_status)]
    if "fail" in statuses:
        status = "fail"
    elif "warn" in statuses:
        status = "warn"
    else:
        status = "pass"
    metrics["status"] = status
    metrics["warnings"] = list(dict.fromkeys(list(metrics.get("warnings", [])) + warnings))
    metrics["errors"] = list(dict.fromkeys(list(metrics.get("errors", [])) + errors))
    metrics["summary"] = _summary(metrics, validation_summary, calibration_enabled, phenotype_present)
    utils.write_json(stage / "metrics.json", metrics)
    ctx.logger.info("校准与校验完成，状态 %s", status)
    return metrics


def _summary(metrics: dict[str, Any], validation: dict[str, int], enabled: bool,
             pheno: bool) -> str:
    def num(value: Any, digits: int = 3) -> str:
        try:
            number = float(value)
            return f"{number:.{digits}f}" if number == number else "—"
        except (TypeError, ValueError):
            return "—"

    lines: list[str] = []
    if metrics.get("null"):
        null = metrics["null"]
        type1 = null.get("type1", {})
        lines.append(f"零模型校准：{null.get('n_replicates')} 次重复，λ 均值 "
                     f"{num(null.get('lambda_mean'))}（SD {num(null.get('lambda_sd'))}）；"
                     f"α=0.05 经验假阳性率 {num(type1.get('alpha_0.05'), 4)}。")
    elif enabled and pheno:
        lines.append("零模型校准：未获得结果。")
    if metrics.get("power"):
        power = metrics["power"]
        alpha5 = (power.get("alpha_5e-8") or {}).get("mean")
        alpha4 = (power.get("alpha_1e-4") or {}).get("mean")
        lines.append(f"已知 QTL 检出力：h²={power.get('h2')}、{power.get('n_causal')} 个因果位点，"
                     f"5e-8 阈值平均 power={num(alpha5)}，1e-4 阈值平均 power={num(alpha4)}。")
    if metrics.get("concordance"):
        conc = metrics["concordance"]
        lines.append(f"因果位点排名一致性：AUC={num(conc.get('rank_auc'))}，"
                     f"top100 命中率={num(conc.get('top_causal_in_top100'))}。")
    if validation.get("n_checks"):
        lines.append(f"独立重算校验：{validation['n_checks']} 项，通过 {validation['n_pass']}、"
                     f"警告 {validation['n_warn']}、失败 {validation['n_fail']}。")
    if metrics.get("status") == "fail":
        lines.append("警告：校准/校验发现失败项，结果需谨慎使用或修复后重跑。")
    return "\n".join(lines)
