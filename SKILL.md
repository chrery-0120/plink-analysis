---
name: plink-analysis
description: 用 PLINK 1.9 完成科研级全流程遗传分析：数据接入与预检、严格质控、群体结构/亲缘/LD/ROH/FST、GWAS（数量/二元性状）、仿真校准与独立重算校验，并按阶段输出中文报告。适用于动物（猪/牛/羊/鸡等）或一般群体的芯片/重测序基因型数据；用户提供 bed/ped/VCF 数据（可附带表型、协变量、分组）时使用。
---

# PLINK 全流程分析（科研级）

本 skill 把"拿到基因型数据 → 可写进论文的结果"拆成 5 个可验证阶段。所有阈值来自 `run.json`，所有 PLINK 调用都有命令日志，每阶段产出 `metrics.json` + 中文 `report.md`，最后汇总为 `05_reports/final_report.md/.html`。

## 何时使用

用户要做 PLINK 相关分析（GWAS、群体结构、亲缘、LD、ROH、FST、质控），或给了一批基因型数据不知道怎么开始。数据没有表型时 flow 自动跳过 GWAS；没有分组时跳过 FST。

## 每次运行的前置检查（按顺序）

1. **定位本 skill**：`<skill_base>/scripts/`、`<skill_base>/references/`（skill_base 即加载本 skill 时给出的目录）。
2. **检查 PLINK**：`python scripts/install_plink.py --check-only`；若缺失执行 `python scripts/install_plink.py`（自动处理官方 S3 域名被 DNS 污染的情况），或用 `--plink-exe`/`PLINK_EXE` 指定已有 `plink.exe`。
3. **确认输入**（缺哪类就跳过对应模块，不要臆造）：
   - 基因型：目录/前缀/单文件（`*.bed+bim+fam`、`*.ped+map`、`*.vcf[.gz]`、`*.tped+tfam`、`*.gen+sample`、23andMe txt）
   - 表型：`FID IID PHENO`（二元 1=对照/2=病例，或 0/1 自动转换；缺失 `-9/NA`）
   - 协变量：`FID IID 数值列...`（分类列自动哑变量编码）
   - 分组：`FID IID GROUP`（用于 FST、分组着色、分组杂合度）
   - 物种：默认 `auto`（从染色体命名推断常染色体数与 XY/ZW）；不确定就保持 auto
   - 数据规模：`--threads`/`--memory` 按机器调整（内存参数是纯 MB 整数）
4. **固定参数再跑**：默认阈值写在 `scripts/pipeline/config.py`，会冻结进 `run.json`。运行开始后不得根据结果反复改阈值；确需修改时先改 `run.json` 并在汇报中说明原因。

## 一键运行

```powershell
python scripts/plink_pipeline.py auto `
  --input "<基因型目录或前缀>" `
  --pheno "<表型文件>" --covar "<协变量文件>" --group "<分组文件>" `
  --out-root "<workspace>/plink_runs" --run-id <项目名> --force
```

没有表型/协变量/分组就删掉对应参数。只跑部分阶段：`--stages intake,qc`。
运行目录：`<out-root>/<run_id>/`：

| 目录 | 内容 |
|---|---|
| `run.json` | 冻结的输入、物种、阈值、版本 |
| `state.json` | 各阶段状态、耗时、产物、错误 |
| `logs/pipeline.log`、`logs/commands.jsonl` | 流程日志与每次 PLINK 调用（可复现） |
| `00_intake/` | `manifest.json`、`preflight.json`、`base_prefix.txt`、清洗后的 `pheno_clean.txt`/`covar_clean.txt`/`group_clean.txt` |
| `01_qc/` | 清洗后 `data/clean.*`、LD 剪枝集、PCA、表格、图、`report.md` |
| `02_popgen/` | PCA、亲缘矩阵、IBD、LD 衰减、ROH/FROH、FST、`report.md` |
| `03_gwas/` | `tables/gwas_results.csv`、`clumps.csv`、Manhattan/QQ、`report.md` |
| `04_calibration/` | `simulate` 校准结果、`validation.json` 独立校验、`report.md` |
| `05_reports/` | `final_report.md`、`final_report.html`、`summary.json` |

## 分段汇报协议（必须执行）

阶段：`intake → qc → popgen → gwas → calibration`。每阶段结束后：

1. 读取该阶段 `metrics.json`（`summary`、`status`、`gates`、`warnings`、`errors`）和 `report.md`。
2. 用中文向用户汇报 5–10 行，包含：**状态**（PASS/WARN/FAIL/ERROR）、**关键数字**（样本与位点保留率、λGC、显著位点数、校准 λ/功率）、**异常与处理**、**对下一步和结论的影响**、**产物路径**。先结论后数字，不贴大段日志。
3. 遇到 `FAIL`/`ERROR`：
   - 读 `state.json` 中该阶段的 `error`、错误里的 `diagnosis`，并查 `references/05-error-playbook.md`；
   - 修复输入/参数后 `python scripts/plink_pipeline.py run --run-dir <run_dir> --stage <stage>` 续跑；
   - 禁止"带病继续"、禁止手工删除样本/位点掩盖问题；任何阈值放宽必须记录并在汇报中说明。
4. 全部完成后，向用户交付：`final_report.md/.html`、GWAS 显著位点表、QC 前后统计、PCA/LD/ROH 图，并明确说明局限（样本量、群体结构、PLINK 1.9 无混合模型等）。

## 各阶段要点与判定

- **intake**：格式识别、染色体/物种推断、ID 与表型/协变量/分组匹配、缺失与编码检查、文件校验和。`preflight` 出现 `fail` 必须先修复。
- **qc**：先位点过滤（缺失率 `--geno`、MAF、HWE）再样本过滤（`--mind`、杂合度 ±3SD、性别、重复 PI_HAT≥0.9、近亲 PI_HAT≥0.25、PCA 离群）。门限：样本保留率 ≥90% PASS、<60% FAIL；位点保留率 ≥80% PASS、<30% FAIL；平均缺失率、HWE 失败比例、性别/重复/离群比例均有独立门限。任何被删样本/位点都在 `tables/` 可追溯。
- **popgen**：PCA（LD 剪枝后）、亲缘矩阵、IBD、LD 衰减、ROH/FROH、FST、分组观察杂合度。用于判断群体分层与是否需要把 PC 放进 GWAS 协变量。
- **gwas**：数量性状 `--linear`、二元性状 `--logistic`，协变量 = 用户协变量 + PC1..k；输出 λGC、QQ/Manhattan、Bonferroni/FDR/建议阈值、lead SNP、LD clump。λGC 0.9–1.1 PASS，0.8–1.3 WARN，之外 FAIL（结构/亲缘/表型问题）。
- **calibration**：①零模型（表型与基因型无关）重复模拟，验收 λ∈[0.97,1.03] 且 α=0.05 经验假阳性率≈0.05；②已知 QTL 模拟（h²=0.3/0.5）估计 power、效应量偏差、排名 AUC；③独立重算校验（MAF/缺失率与 PLINK 对照、GWAS β/SE 与 statsmodels 对照、λ 独立计算）。任一项失败都要在报告中显著标注。

## 常用命令

```powershell
python scripts/plink_pipeline.py status --run-dir <run_dir>     # 查看阶段状态
python scripts/plink_pipeline.py report --run-dir <run_dir>     # 重新生成汇总报告
python scripts/plink_pipeline.py selftest                       # 演示数据端到端自检
```

## 参考文档（按需加载，不要一次全读）

- `references/01-pipeline-stages.md` — 阶段方法、PLINK 命令清单、执行顺序
- `references/02-qc-thresholds.md` — 阈值表与场景化调整
- `references/03-gwas.md` — 模型/协变量/多重检验/λGC/结果可信度
- `references/04-popgen.md` — PCA/亲缘/LD/ROH/FST 的解释
- `references/05-error-playbook.md` — 报错诊断与修复命令
- `references/06-calibration.md` — 仿真校准与独立校验设计
- `references/07-report-template.md` — 报告模板与汇报话术

## 硬性约束

- 不修改原始输入文件；所有产物写入 `plink_runs/<run_id>/`。
- 所有 PLINK 调用必须经 `pipeline/utils.py` 的 `Plink` 包装（自动记录命令、解析日志、失败给诊断）。
- 阈值、随机种子、软件版本必须落盘（`run.json`、`state.json`）后才用于结论。
- 报告中的数字必须来自产物文件，不得手写估计值；不确定就标注"需复核"。
- 数据涉及动物育种/保种时，注意：小群体近交、品种混合导致的结构、芯片密度对 LD/ROH 的影响；解释时引用 `references/04-popgen.md`。
