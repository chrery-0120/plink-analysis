# plink-analysis

[![CI](https://github.com/chrery-0120/plink-analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/chrery-0120/plink-analysis/actions/workflows/ci.yml)

**科研级 PLINK 1.9 全流程自动化分析 Skill** —— 从原始基因型数据到可写进论文的结果与中文报告。

A research-grade, end-to-end PLINK 1.9 pipeline packaged as an agent **skill**: data intake → strict QC → population structure / relatedness / LD / ROH / FST → GWAS → simulation-based calibration + independent recomputation → staged Chinese reports. Built for animal breeding / population genetics (pig, cattle, sheep, goat, chicken, …) and general cohorts.

---

## 特性

- **一条命令跑完 5 个阶段**：`intake → qc → popgen → gwas → calibration → report`
- **数据格式自动识别**：PLINK bed/ped/tped、VCF(.gz)、GEN/SAMPLE、23andMe；染色体命名与物种自动推断（猪/牛/羊/鸡等 `--chr-set` 规则自动处理）
- **严格质控**：先位点（缺失率 / MAF / HWE）后样本（缺失率 / 杂合度 ±3SD / 性别 / 重复与近亲 PI_HAT / PCA 离群），LD 剪枝，全部过滤可追溯
- **群体遗传**：PCA、亲缘矩阵、IBD、LD 衰减、ROH/FROH、FST、分组杂合度
- **GWAS**：数量性状 `--linear` / 二元性状 `--logistic`，PCA + 用户协变量，λGC、QQ/Manhattan、Bonferroni/FDR/建议阈值、lead SNP、LD clump
- **准确性校准**：零模型（λ、经验假阳性率）+ 已知 QTL（power、效应量、排名 AUC）仿真校准；独立重算校验（MAF/缺失率/β/SE/λGC 与 Python/statsmodels 对照）
- **自动查错**：预检 + PLINK 日志解析 + 中文诊断与可复制修复命令；失败可断点续跑
- **可复现**：阈值冻结在 `run.json`；每次 PLINK 调用记录在 `logs/commands.jsonl`；随机种子固定
- **中文报告**：每阶段 `report.md` + 自包含 `final_report.html`（图片内嵌）

## 目录结构

```
plink-analysis/
├── SKILL.md                 # Agent 主指令：前置检查、分段执行、分段汇报协议、失败处理
├── DESIGN.md                # 设计文档：阶段、模块契约、阈值与验收标准
├── LICENSE                  # MIT
├── .github/workflows/ci.yml # 每次 push/PR 跑端到端自检
├── examples/demo/           # 演示数据跑出的示例报告（final_report.html/.md）
├── references/              # 中文方法学文档（按需加载）
│   ├── 01-pipeline-stages.md
│   ├── 02-qc-thresholds.md
│   ├── 03-gwas.md
│   ├── 04-popgen.md
│   ├── 05-error-playbook.md
│   ├── 06-calibration.md
│   └── 07-report-template.md
├── scripts/
│   ├── plink_pipeline.py    # CLI：auto / run / status / report / selftest
│   ├── install_plink.py     # 自动下载安装 PLINK 1.9（含 DNS 污染绕过）
│   ├── simulate_demo.py     # 生成可复现演示数据（含 truth）
│   ├── selftest.py          # 一键端到端验收
│   └── pipeline/            # 5 个阶段 + 独立校验 + 仿真校准 + 报告/绘图
└── assets/report.css
```

## 安装

### 1) 作为 Agent Skill 使用

把本仓库放到 DSH / Claude 等支持 skill 的工作区扫描目录下（目录名即 skill 名）：

```bash
git clone <this-repo> <workspace>/.agents/skills/plink-analysis
# 或
git clone <this-repo> <workspace>/.dsh/skills/plink-analysis
```

之后在对话中直接说"用 plink-analysis 分析 `<数据路径>`"即可，Agent 会按 `SKILL.md` 的分段协议执行并逐段汇报。

### 2) 依赖

- Python ≥ 3.10（已验证 3.12；`requirements.txt` 已按验证通过的版本锁定）
- `pip install -r requirements.txt`
- PLINK 1.9（本机没有时）：

```bash
python scripts/install_plink.py            # 按平台自动下载（Windows/Linux/macOS）到 <workspace>/tools/plink/
# 或指定已有二进制
python scripts/plink_pipeline.py --plink-exe /path/to/plink ...
```

> 安装脚本会优先使用官方 S3；若 DNS 被污染，会用 DNS-over-HTTPS + 并行 Range 下载。

## 快速开始

```bash
python scripts/plink_pipeline.py auto \
  --input /data/pig/geno \
  --pheno /data/pig/pheno.txt \
  --covar /data/pig/covar.txt \
  --group /data/pig/group.txt \
  --out-root /data/plink_runs --run-id pig_2026
```

没有表型/协变量/分组就删掉对应参数（流程会自动跳过 GWAS / FST / 分组分析）。
只跑部分阶段：`--stages intake,qc`；续跑：`run --run-dir <dir> --stage qc`。

### 输入格式

| 文件 | 格式 |
|---|---|
| 基因型 | `*.bed + .bim + .fam`、`*.ped + .map`、`*.vcf[.gz]`、`*.tped + .tfam`、`*.gen + .sample`、23andMe txt |
| 表型 | `FID IID PHENO`（二元：1=对照、2=病例；也接受 0/1 自动转换；缺失 `-9/NA`） |
| 协变量 | `FID IID 数值列...`（文本列自动哑变量编码） |
| 分组 | `FID IID GROUP`（品种/家系/群体） |

## 输出

```
<out-root>/<run_id>/
├── run.json                  # 冻结的输入、物种、阈值、版本
├── state.json                # 各阶段状态、耗时、产物
├── logs/                     # pipeline.log + commands.jsonl（完整命令可复现）
├── 00_intake/                # manifest / preflight / 清洗后的表型·协变量·分组
├── 01_qc/                    # clean bfile、统计表、图、report.md
├── 02_popgen/                # PCA / 亲缘 / LD / ROH / FST
├── 03_gwas/                  # gwas_results.csv、clumps.csv、Manhattan/QQ
├── 04_calibration/           # 仿真校准 + validation.json
└── 05_reports/               # final_report.md / final_report.html / summary.json
```

## 准确性与验证

演示数据由 `scripts/simulate_demo.py` 按固定种子生成（400 样本 × 8000 SNP，2 个亚群、10 对重复/亲缘、2% 缺失、30 个已知因果 SNP），一键自检：

```bash
python scripts/selftest.py --with-binary
```

本机验收结果（2026-10-23）：

- 自检 **34 项：33 PASS / 0 FAIL / 1 SKIP**（二元性状 5e-8 power 按方法学原因不作硬门槛），约 2 分钟
- 独立重算：MAF 用整数等位基因计数逐位点一致；缺失率 N_MISS 精确一致；PCA 方向相关 0.999998；线性 β/SE 与 statsmodels 最大差 4.8e-5；λGC 独立复算一致
- 零模型校准：λ = 0.998（SD 0.035），α=0.05 经验假阳性率 0.0488
- 已知 QTL（h²=0.5）：5e-8 any-hit = 1.00，rank AUC = 0.773，效应量回归斜率 1.035
- 主 GWAS：λGC = 1.063，Bonferroni 显著 3 个，lead SNP 出现在报告中

2026-11-07 复核：按 `requirements.txt` 锁定版本（Python 3.12.3、PLINK 1.9.0 2026-09-27 构建）复跑，数量性状自检 20/20 PASS，关键数字与上表一致。

## 示例报告

`examples/demo/` 是演示数据跑出的完整报告：`final_report.html`（自包含，下载后用浏览器打开即可）与 `final_report.md`。不用跑流程就能看到输出长什么样。

## 常用命令

```bash
python scripts/plink_pipeline.py status --run-dir <run_dir>   # 阶段状态
python scripts/plink_pipeline.py report --run-dir <run_dir>   # 重新生成汇总报告
python scripts/selftest.py --with-binary                      # 端到端验收
```

## 局限

- PLINK 1.9 没有混合模型（LMM）：强亲缘/强分层数据建议再用 GCTA / SAIGE / GMMAT 复核，本流程用 PCA 协变量 + 亲缘诊断提示风险
- 芯片密度、品种混合、小样本会影响 LD/ROH/FST/λGC 的解释
- 显著 ≠ 因果；结果解释需结合系谱、品种背景与生物学知识
- 参考基因组版本与 SNP 坐标需一致；基因注释需自备注释文件

## License

MIT License，见 [LICENSE](LICENSE)。
