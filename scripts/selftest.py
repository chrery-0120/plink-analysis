#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""plink-analysis 一键自检（DESIGN 第 14 节验收）。

做什么
------
1. 生成（或校验复用）演示数据集 ``E:\\deepseek\\demo_data``；
2. 用主 CLI 跑完整流程：``plink_pipeline.py auto --input <demo/geno>
   --pheno pheno.txt --covar covar.txt --group group.txt --run-id selftest ...``；
3. 断言各阶段状态、报告/metrics 产物、独立重算校验（validate）、仿真校准（λ/FPR/
   power/排序）、GWAS lead SNP、``status`` 子命令与断点续跑；
4. 打印 PASS/FAIL 清单（含定位信息：文件路径、命令、日志尾部）与总耗时。

用法（在 skill 目录下）::

    python scripts/selftest.py                       # 数量性状全流程（默认）
    python scripts/selftest.py --with-binary         # 追加二元性状流程
    python scripts/selftest.py --check-resume        # 追加断点续跑检查
    python scripts/selftest.py --keep-run            # 保留本次 run 目录
    python scripts/selftest.py --run-id mytest --out-root E:\\deepseek\\plink_runs

退出码：0 = 全部 PASS；1 = 有 FAIL；2 = 自检脚本自身异常/环境缺失。
目标耗时 ≤ 15 分钟（默认 demo 规模约 1~2 分钟）。
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = SKILL_DIR / "scripts"
CLI = SCRIPTS_DIR / "plink_pipeline.py"
SIMULATE_DEMO = SCRIPTS_DIR / "simulate_demo.py"
def _workspace() -> Path:
    """<workspace> when the skill lives in <workspace>/.agents/skills/plink-analysis."""
    return Path(__file__).resolve().parents[4]


def _default_plink() -> Path:
    for env_name in ("PLINK_EXE", "PLINK"):
        env = os.environ.get(env_name)
        if env:
            return Path(env)
    exe = "plink.exe" if os.name == "nt" else "plink"
    for candidate in (_workspace() / "tools" / "plink" / exe,
                      Path.cwd() / "tools" / "plink" / exe,
                      Path(exe)):
        if candidate.exists():
            return candidate
    return Path(exe)


DEFAULT_DEMO_DIR = _workspace() / "demo_data"
DEFAULT_OUT_ROOT = _workspace() / "plink_runs"
DEFAULT_PLINK = _default_plink()
DEMO_SEED = 20261023

# 验收阈值
LAMBDA_RANGE = (0.95, 1.05)        # 零模型 λ 均值（Lead/methods 共识的宽容差）
FPR05_RANGE = (0.02, 0.08)         # α=0.05 经验假阳性率
RANK_AUC_MIN = 0.65
TOP100_MIN = 0.10
MIN_VALIDATE_CHECKS = 12
MAX_VALIDATE_WARN = 2              # 已知可容忍项（如相关对计数口径差异）


# ---------------------------------------------------------------------------
# 检查器
# ---------------------------------------------------------------------------
@dataclass
class Check:
    name: str
    status: str                      # PASS / FAIL / SKIP
    detail: str
    hint: str = ""


@dataclass
class Report:
    checks: list[Check] = field(default_factory=list)

    def add(self, name: str, ok: bool | None, detail: str, hint: str = "") -> None:
        if ok is None:
            status = "SKIP"
        else:
            status = "PASS" if ok else "FAIL"
        self.checks.append(Check(name, status, detail, hint))

    @property
    def n_fail(self) -> int:
        return sum(1 for c in self.checks if c.status == "FAIL")

    @property
    def n_pass(self) -> int:
        return sum(1 for c in self.checks if c.status == "PASS")

    @property
    def n_skip(self) -> int:
        return sum(1 for c in self.checks if c.status == "SKIP")


def _load_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except Exception:
        return None


def _run(cmd: list[str], *, timeout: int = 1800, cwd: Path | None = None) -> dict[str, Any]:
    """运行子进程并返回 rc/stdout/stderr/elapsed（UTF-8）。"""
    started = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                              errors="replace", timeout=timeout,
                              cwd=str(cwd) if cwd else None)
        return {"rc": proc.returncode, "stdout": proc.stdout or "",
                "stderr": proc.stderr or "", "elapsed": time.time() - started,
                "cmd": " ".join(map(str, cmd)), "timeout": False}
    except subprocess.TimeoutExpired as exc:
        return {"rc": -1, "stdout": (exc.stdout or "") if isinstance(exc.stdout, str) else "",
                "stderr": f"超时 {timeout}s", "elapsed": time.time() - started,
                "cmd": " ".join(map(str, cmd)), "timeout": True}


def _tail(text: str, n: int = 12) -> str:
    lines = [ln for ln in (text or "").strip().splitlines() if ln.strip()]
    return "\n".join(lines[-n:])


# ---------------------------------------------------------------------------
# 1. 演示数据
# ---------------------------------------------------------------------------
def ensure_demo(demo_dir: Path, report: Report, *, force: bool = False) -> bool:
    required = ["geno.bed", "geno.bim", "geno.fam", "pheno.txt", "pheno_bin.txt",
                "covar.txt", "group.txt", "truth.csv", "truth_meta.json"]
    missing = [f for f in required if not (demo_dir / f).exists()]
    meta = _load_json(demo_dir / "truth_meta.json") or {}
    demo_ok = (not missing and int(meta.get("n_samples", 0)) >= 300
               and int(meta.get("n_snps", 0)) >= 5000
               and int(meta.get("seed", -1)) == DEMO_SEED)
    if demo_ok and not force:
        report.add("demo_dataset", True,
                   f"复用演示数据 {demo_dir}（{meta.get('n_samples')} 样本 x "
                   f"{meta.get('n_snps')} 位点，seed={meta.get('seed')}）")
        return True
    reason = ("缺少文件 " + ", ".join(missing[:4]) if missing
              else f"元数据不符（seed={meta.get('seed')}, n={meta.get('n_samples')}x"
                   f"{meta.get('n_snps')}）")
    print(f"[selftest] 需要重新生成演示数据：{reason}")
    if not SIMULATE_DEMO.exists():
        report.add("demo_dataset", False, f"缺少生成脚本 {SIMULATE_DEMO}",
                   hint="检查 skill 目录完整性")
        return False
    res = _run([sys.executable, str(SIMULATE_DEMO), "--out", str(demo_dir),
                "--seed", str(DEMO_SEED)], timeout=900)
    ok = res["rc"] == 0
    detail = (f"生成演示数据 {demo_dir}（{res['elapsed']:.1f}s）" if ok
              else f"演示数据生成失败 rc={res['rc']}")
    report.add("demo_dataset", ok, detail,
               hint="命令: " + res["cmd"] + "\n日志尾部:\n" + _tail(res["stderr"] or res["stdout"]))
    meta = _load_json(demo_dir / "truth_meta.json") or {}
    verify = meta.get("plink_verify") or {}
    report.add("demo_plink_verify", str(verify.get("status", "skipped")) == "pass",
               f"演示数据 PLINK 读取验证: {verify.get('status')}；"
               f"truth 因果位点 {meta.get('n_causal')} 个，h2={meta.get('h2_quantitative')}")
    return ok


# ---------------------------------------------------------------------------
# 2. 跑主流程
# ---------------------------------------------------------------------------
def run_pipeline(demo_dir: Path, out_root: Path, run_id: str, *, with_binary: bool,
                 timeout: int, report: Report) -> list[dict[str, Any]]:
    runs: list[dict[str, Any]] = []
    jobs = [("quantitative", "pheno.txt")]
    if with_binary:
        jobs.append(("binary", "pheno_bin.txt"))
    for label, pheno in jobs:
        rid = run_id if label == "quantitative" else f"{run_id}_bin"
        run_dir = out_root / rid
        if run_dir.exists():
            try:
                shutil.rmtree(run_dir, ignore_errors=True)
            except Exception as exc:  # pragma: no cover - Windows 文件锁
                print(f"[selftest] 警告：清理旧 run 目录失败（{exc}），改用新的 run id")
            if run_dir.exists():
                rid = f"{rid}_{_dt.datetime.now():%H%M%S}"
                run_dir = out_root / rid
        cmd = [sys.executable, str(CLI), "auto",
               "--input", str(demo_dir / "geno"),
               "--pheno", str(demo_dir / pheno),
               "--covar", str(demo_dir / "covar.txt"),
               "--group", str(demo_dir / "group.txt"),
               "--run-id", rid, "--out-root", str(out_root)]
        print(f"[selftest] 运行 {label} 流程: {run_id if label == 'quantitative' else rid}")
        res = _run(cmd, timeout=timeout)
        run_dir = out_root / rid
        res.update({"label": label, "run_dir": run_dir, "run_id": rid})
        runs.append(res)
        ok = res["rc"] in (0, 1) and (run_dir / "state.json").exists()
        report.add(f"pipeline_{label}", ok,
                   f"CLI auto 返回码 {res['rc']}（{res['elapsed']:.1f}s），run_dir={run_dir}",
                   hint="命令: " + res["cmd"] + "\n输出尾部:\n" + _tail(res["stdout"]))
        if not ok:
            report.add(f"pipeline_{label}_diagnosis", False,
                       "主流程未成功，后续阶段检查将基于缺失产物判定",
                       hint=_tail(res["stdout"] + "\n" + res["stderr"], 20))
    return runs


# ---------------------------------------------------------------------------
# 3. 逐项验收检查
# ---------------------------------------------------------------------------
STAGE_NAMES = ["intake", "qc", "popgen", "gwas", "calibration"]
STAGE_CN = {"intake": "数据接入", "qc": "质控", "popgen": "群体结构",
            "gwas": "GWAS", "calibration": "校准与校验"}


def check_run(run: dict[str, Any], demo_dir: Path, report: Report,
              prefix: str) -> None:
    run_dir: Path = run["run_dir"]
    label = run["label"]
    state = _load_json(run_dir / "state.json") or {}
    stages = state.get("stages") or {}

    # 阶段状态
    bad = []
    for stage in STAGE_NAMES:
        st = (stages.get(stage) or {}).get("status")
        if st is None:
            bad.append(f"{stage}=缺失")
        elif st in {"fail", "error"}:
            bad.append(f"{stage}={st}")
        elif st not in {"pass", "warn", "skipped"}:
            bad.append(f"{stage}={st}")
    report.add(f"{prefix}_stages_status", not bad,
               (f"阶段状态: " + ", ".join(f"{s}={(stages.get(s) or {}).get('status')}"
                                          for s in STAGE_NAMES))
               + ("" if not bad else f"；异常项 {bad}"),
               hint=f"state.json: {run_dir / 'state.json'}")

    # 阶段报告 / metrics / preflight
    required = {
        "00_intake/report.md": "数据接入报告",
        "00_intake/preflight.json": "接入预检",
        "00_intake/manifest.json": "数据清单",
        "01_qc/report.md": "质控报告",
        "02_popgen/report.md": "群体结构报告",
        "03_gwas/report.md": "GWAS 报告",
        "04_calibration/report.md": "校准报告",
        "04_calibration/validation.json": "独立校验结果",
        "05_reports/final_report.md": "最终报告(md)",
        "05_reports/final_report.html": "最终报告(html)",
        "05_reports/summary.json": "汇总 JSON",
    }
    missing = [name for name in required
               if not (run_dir / name).exists()
               or (run_dir / name).stat().st_size == 0]
    report.add(f"{prefix}_stage_reports", not missing,
               "所有阶段 report.md 与最终报告齐全" if not missing
               else f"缺失/空文件: {missing}",
               hint=f"run_dir={run_dir}")

    metrics_paths = {
        "intake": "00_intake/metrics.json",
        "qc": "01_qc/metrics.json",
        "popgen": "02_popgen/metrics.json",
        "gwas": "03_gwas/metrics.json",
        "calibration": "04_calibration/metrics.json",
    }
    bad_metrics = []
    loaded: dict[str, Any] = {}
    for stage, rel in metrics_paths.items():
        obj = _load_json(run_dir / rel)
        if obj is None:
            bad_metrics.append(rel)
            continue
        loaded[stage] = obj
        if obj.get("stage") not in (None, stage):
            bad_metrics.append(f"{rel}:stage={obj.get('stage')}")
    report.add(f"{prefix}_metrics_json", not bad_metrics,
               "各阶段 metrics.json 可解析" if not bad_metrics else f"异常: {bad_metrics}")

    # DESIGN 第 8 节关键字段
    qc = loaded.get("qc") or {}
    gwas = loaded.get("gwas") or {}
    cal = loaded.get("calibration") or {}
    contract_missing = []
    status_only: list[str] = []
    for stage, obj, keys in (
        ("qc", qc, ("n_samples_in", "n_samples_out", "n_variants_in", "n_variants_out",
                    "gates", "overall")),
        ("gwas", gwas, ("model", "n_tested", "lambda_gc", "lead_snps", "overall")),
        ("calibration", cal, ("null", "power", "validation", "overall")),
    ):
        for key in keys:
            if key in obj:
                continue
            # 兼容：stage 目前用 status 表达总体判定，overall 是其别名
            if key == "overall" and "status" in obj:
                status_only.append(f"{stage}.overall→status")
                continue
            contract_missing.append(f"{stage}.{key}")
    report.add(f"{prefix}_metrics_contract", not contract_missing,
               ("metrics 字段符合 DESIGN 第 8 节" if not contract_missing
                else f"缺少字段: {contract_missing}")
               + (f"（以 status 代替 overall: {status_only}）" if status_only else ""))

    # 独立重算校验
    validation = _load_json(run_dir / "04_calibration" / "validation.json") or {}
    checks = validation.get("checks") or []
    n_fail = int(validation.get("n_fail", sum(1 for c in checks if c.get("status") == "fail")))
    n_warn = int(validation.get("n_warn", sum(1 for c in checks if c.get("status") == "warn")))
    ok_val = (len(checks) >= MIN_VALIDATE_CHECKS and n_fail == 0 and n_warn <= MAX_VALIDATE_WARN)
    fail_names = [c.get("name") for c in checks if c.get("status") == "fail"]
    warn_names = [c.get("name") for c in checks if c.get("status") == "warn"]
    report.add(f"{prefix}_validation", ok_val,
               (f"独立重算校验 {len(checks)} 项：pass={len(checks) - n_fail - n_warn} "
                f"warn={n_warn} fail={n_fail}（overall={validation.get('overall')}）")
               + (f"；FAIL: {fail_names}" if fail_names else "")
               + (f"；WARN: {warn_names}" if warn_names else ""),
               hint=f"validation.json: {run_dir / '04_calibration' / 'validation.json'}；"
                    f"明细表: {run_dir / '04_calibration' / 'tables' / 'validation_checks.csv'}")
    # 关键容差项必须 pass
    key_checks = {"maf_recompute", "missingness_recompute", "gwas_beta_se_recompute",
                  "lambda_gc_recompute", "pca_recompute", "structure_counts"}
    status_of = {c.get("name"): c.get("status") for c in checks}
    bad_key = [k for k in key_checks
               if status_of.get(k) not in ("pass", None)]
    report.add(f"{prefix}_validation_key_checks", not bad_key,
               f"关键重算项 {sorted(key_checks)} 全部 pass" if not bad_key
               else f"未通过: {[(k, status_of.get(k)) for k in bad_key]}")

    # 校准数值验收
    null = cal.get("null") or {}
    lam_mean = null.get("lambda_mean")
    fpr05 = (null.get("type1") or {}).get("alpha_0.05")
    ok_lam = lam_mean is not None and LAMBDA_RANGE[0] <= float(lam_mean) <= LAMBDA_RANGE[1]
    ok_fpr = fpr05 is not None and FPR05_RANGE[0] <= float(fpr05) <= FPR05_RANGE[1]
    report.add(f"{prefix}_calibration_lambda", ok_lam,
               f"零模型 λ 均值 {float(lam_mean):.3f}（要求 ∈ {LAMBDA_RANGE}，"
               f"{null.get('n_replicates')} 次重复）" if lam_mean is not None else
               "零模型 λ 均值为空",
               hint=f"04_calibration/metrics.json → null.lambda_mean")
    report.add(f"{prefix}_calibration_fpr", ok_fpr,
               f"α=0.05 经验假阳性率 {float(fpr05):.4f}（要求 ∈ {FPR05_RANGE}）"
               if fpr05 is not None else "α=0.05 经验假阳性率为空")
    power = cal.get("power") or {}
    any_hit = ((power.get("alpha_5e-8") or {}).get("any_hit_rate"))
    rank_auc = (cal.get("concordance") or {}).get("rank_auc")
    top100 = (cal.get("concordance") or {}).get("top_causal_in_top100")
    model = str(cal.get("model") or "linear")
    if model == "logistic":
        # 二元性状：病例数少（~20%）+ log-OR 尺度，5e-8 的 any-hit 概率本就低，
        # 按 methods-writer 口径不作硬门槛（只记录）；以排序 AUC/λ/FPR 为验收。
        report.add(f"{prefix}_calibration_power", None,
                   f"二元性状（logistic，{cal.get('n_samples')} 样本）5e-8 至少检出 1 个因果 "
                   f"SNP 的重复比例 {any_hit}：不作硬性门槛（方法学上 n=400 时检出概率很低）")
    else:
        ok_power = (any_hit is not None and float(any_hit) >= 0.5)
        report.add(f"{prefix}_calibration_power", ok_power,
                   f"已知 QTL（h²={power.get('h2')}）5e-8 至少检出 1 个因果 SNP 的重复比例 "
                   f"{float(any_hit):.2f}（要求 ≥0.5）" if any_hit is not None else
                   "power.any_hit_rate 为空")
    ok_auc = rank_auc is not None and float(rank_auc) >= RANK_AUC_MIN
    report.add(f"{prefix}_calibration_concordance", ok_auc,
               f"因果 SNP 排序 AUC={float(rank_auc):.3f}（要求 ≥{RANK_AUC_MIN}），"
               f"top100 命中率={float(top100):.3f}" if rank_auc is not None
               else "因果 SNP 排序 AUC 为空")
    eb = cal.get("effect_bias") or {}
    if model == "logistic":
        sign_conc = eb.get("sign_concordance")
        report.add(f"{prefix}_calibration_effect_bias",
                   sign_conc is not None and float(sign_conc) >= 0.70,
                   f"logistic 效应方向一致率={sign_conc}（要求 ≥0.70；"
                   f"log-OR 与 liability β 尺度不同，只比方向）")
    else:
        ratio = eb.get("mean_beta_ratio")
        report.add(f"{prefix}_calibration_effect_bias",
                   ratio is not None and 0.7 <= float(ratio) <= 1.3,
                   f"效应量回归斜率 Σ(β_est·β_true)/Σ(β_true²)={ratio}（合理区间 0.7~1.3）；"
                   f"SE 标定比={eb.get('se_calibration')}")

    # GWAS lead SNP 与报告
    leads = gwas.get("lead_snps") or []
    lead_ok = bool(leads)
    gwas_report = (run_dir / "03_gwas" / "report.md")
    if lead_ok and gwas_report.exists():
        text = gwas_report.read_text(encoding="utf-8", errors="replace")
        lead_snp = str(leads[0].get("snp", ""))
        lead_ok = lead_snp in text
        lead_detail = (f"lead SNP {lead_snp}（P={leads[0].get('p')}）已出现在 GWAS 报告；"
                       f"λGC={float(gwas.get('lambda_gc')):.3f}"
                       if gwas.get("lambda_gc") is not None else
                       f"lead SNP {lead_snp} 已出现在 GWAS 报告")
    else:
        lead_detail = "lead_snps 为空或 GWAS 报告缺失"
    report.add(f"{prefix}_gwas_lead_snp", lead_ok, lead_detail,
               hint=f"03_gwas/metrics.json → lead_snps；报告 {gwas_report}")
    lam_gwas = gwas.get("lambda_gc")
    report.add(f"{prefix}_gwas_lambda",
               lam_gwas is not None and 0.85 <= float(lam_gwas) <= 1.25,
               f"主 GWAS λGC={float(lam_gwas):.3f}（合理区间 0.85~1.25；多基因性状可略>1）"
               if lam_gwas is not None else "主 GWAS λGC 为空")

    # 结果表可追溯
    res_csv = run_dir / "03_gwas" / "tables" / "gwas_results.csv"
    ok_res = res_csv.exists() and res_csv.stat().st_size > 0
    n_lines = 0
    if ok_res:
        try:
            n_lines = sum(1 for _ in res_csv.open(encoding="utf-8-sig")) - 1
        except Exception:
            n_lines = 0
    report.add(f"{prefix}_gwas_results_csv", ok_res and n_lines > 0,
               f"gwas_results.csv 存在，{n_lines} 行", hint=str(res_csv))


def check_status_cmd(run_dir: Path, report: Report) -> None:
    res = _run([sys.executable, str(CLI), "status", "--run-dir", str(run_dir)], timeout=120)
    out = (res["stdout"] + "\n" + res["stderr"])
    has_stages = sum(1 for s in STAGE_NAMES if s in out.lower()) >= 4
    ok = res["rc"] in (0, 1, 2) and has_stages
    report.add("status_command", ok,
               f"status 命令返回码 {res['rc']}，输出含阶段状态"
               if ok else f"status 命令异常 rc={res['rc']}",
               hint="命令: " + res["cmd"] + "\n输出尾部:\n" + _tail(out))


def check_resume(out_root: Path, run_id: str, report: Report, timeout: int) -> None:
    """断点续跑：删除 03_gwas 阶段目录后用 run --stage gwas 重跑。"""
    run_dir = out_root / run_id
    gwas_dir = run_dir / "03_gwas"
    if not gwas_dir.exists():
        report.add("resume_stage", False, f"03_gwas 不存在，无法测试续跑", hint=str(gwas_dir))
        return
    backup = run_dir / "03_gwas__bak"
    if backup.exists():
        shutil.rmtree(backup)
    shutil.move(str(gwas_dir), str(backup))
    try:
        res = _run([sys.executable, str(CLI), "run", "--run-dir", str(run_dir),
                    "--stage", "gwas", "--force"], timeout=timeout)
        ok = (res["rc"] in (0, 1) and gwas_dir.exists()
              and (gwas_dir / "metrics.json").exists()
              and (gwas_dir / "report.md").exists())
        report.add("resume_stage", ok,
                   f"删除 03_gwas 后 run --stage gwas 重跑：rc={res['rc']}，"
                   f"metrics/report 已重建" if ok
                   else f"断点续跑失败 rc={res['rc']}",
                   hint="命令: " + res["cmd"] + "\n输出尾部:\n" + _tail(res["stdout"]))
    finally:
        if backup.exists() and not gwas_dir.exists():
            shutil.move(str(backup), str(gwas_dir))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="plink-analysis 一键自检",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--demo-dir", default=str(DEFAULT_DEMO_DIR))
    ap.add_argument("--out-root", default=str(DEFAULT_OUT_ROOT))
    ap.add_argument("--run-id", default=None,
                    help="run id（默认 selftest_<YYYYmmdd_HHMMSS>，唯一，避免并发冲突）")
    ap.add_argument("--with-binary", action="store_true", help="额外跑二元性状全流程")
    ap.add_argument("--check-resume", action="store_true", help="额外验证断点续跑")
    ap.add_argument("--force-demo", action="store_true", help="强制重新生成演示数据")
    ap.add_argument("--keep-run", action="store_true", help="保留本次 run 目录")
    ap.add_argument("--timeout", type=int, default=1500, help="单次流程超时（秒）")
    ap.add_argument("--json", action="store_true", help="打印 JSON 摘要")
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    args = parse_args(argv)
    if not args.run_id:
        args.run_id = f"selftest_{_dt.datetime.now():%Y%m%d_%H%M%S}"
    t0 = time.time()
    demo_dir = Path(args.demo_dir)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    print("=" * 78)
    print("plink-analysis 一键自检")
    print(f"  skill:   {SKILL_DIR}")
    print(f"  demo:    {demo_dir}")
    print(f"  runs:    {out_root}")
    print(f"  python:  {sys.version.split()[0]}")
    print("=" * 78)

    report = Report()
    # 环境
    report.add("environment", CLI.exists() and SIMULATE_DEMO.exists(),
               f"CLI={CLI.name} 存在，PLINK={DEFAULT_PLINK}"
               + ("（存在）" if DEFAULT_PLINK.exists() else "（缺失！）"),
               hint="检查 .agents/skills/plink-analysis 目录与 tools/plink 安装")
    try:
        import numpy, pandas, scipy, statsmodels  # noqa: F401
        report.add("python_deps", True, "numpy/pandas/scipy/statsmodels 可用")
    except Exception as exc:  # noqa: BLE001
        report.add("python_deps", False, f"Python 依赖缺失: {exc}")
        return 2

    # 1) 演示数据
    ensure_demo(demo_dir, report, force=args.force_demo)

    # 2) 主流程（数量性状，可选二元）
    runs = run_pipeline(demo_dir, out_root, args.run_id, with_binary=args.with_binary,
                        timeout=args.timeout, report=report)

    # 3) 逐项验收
    for run in runs:
        if run["run_dir"].exists():
            check_run(run, demo_dir, report, run["label"])

    # 4) status 命令
    if runs and runs[0]["run_dir"].exists():
        check_status_cmd(runs[0]["run_dir"], report)

    # 5) 断点续跑（可选）
    if args.check_resume and runs and runs[0]["run_dir"].exists():
        check_resume(out_root, runs[0]["run_id"], report, args.timeout)

    # 6) 清理
    if not args.keep_run:
        for run in runs:
            run_dir = run["run_dir"]
            if run_dir.exists():
                shutil.rmtree(run_dir, ignore_errors=True)
        print(f"[selftest] 已清理 run 目录（--keep-run 可保留）："
              f"{[str(r['run_dir']) for r in runs]}")

    # 汇总
    elapsed = time.time() - t0
    print("=" * 78)
    print("自检结果")
    print("=" * 78)
    for c in report.checks:
        mark = {"PASS": "PASS", "FAIL": "FAIL", "SKIP": "SKIP"}[c.status]
        print(f"[{mark}] {c.name}: {c.detail}")
        if c.status == "FAIL" and c.hint:
            for line in c.hint.splitlines():
                print(f"        {line}")
    print("-" * 78)
    print(f"合计 {len(report.checks)} 项：PASS {report.n_pass}，"
          f"FAIL {report.n_fail}，SKIP {report.n_skip}；总耗时 {elapsed / 60:.1f} 分钟")
    if report.n_fail == 0:
        print("结论：PASS —— plink-analysis 全流程可用（阶段产物、独立校验、仿真校准均达标）")
    else:
        print("结论：FAIL —— 请按上面 FAIL 项的定位信息排查后重跑")
    if args.json:
        print(json.dumps({"n_pass": report.n_pass, "n_fail": report.n_fail,
                          "n_skip": report.n_skip, "elapsed_s": round(elapsed, 1),
                          "checks": [c.__dict__ for c in report.checks]},
                         ensure_ascii=False, indent=2))
    # 自检报告落盘（保留在 out_root 下，便于交付证据）
    try:
        report_path = out_root / f"selftest_{_dt.datetime.now():%Y%m%d_%H%M%S}.md"
        lines = ["# plink-analysis 自检报告", "",
                 f"- 时间：{_dt.datetime.now().astimezone().isoformat(timespec='seconds')}",
                 f"- skill：`{SKILL_DIR}`",
                 f"- demo：`{demo_dir}`",
                 f"- 合计：{len(report.checks)} 项，PASS {report.n_pass}，FAIL {report.n_fail}",
                 f"- 总耗时：{elapsed / 60:.1f} 分钟", "",
                 "| 检查 | 状态 | 详情 |", "|---|---|---|"]
        for c in report.checks:
            lines.append(f"| {c.name} | {c.status} | {c.detail.replace('|', '/')} |")
        if report.n_fail:
            lines += ["", "## FAIL 定位", ""]
            for c in report.checks:
                if c.status == "FAIL":
                    lines += [f"### {c.name}", "", c.detail, "", "```", c.hint, "```", ""]
        report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"自检报告：{report_path}")
    except Exception as exc:  # noqa: BLE001
        print(f"[selftest] 写自检报告失败: {exc}")
    return 0 if report.n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
