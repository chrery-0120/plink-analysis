# plink-analysis 设计文档（DESIGN.md）

本文档是代码中 `DESIGN.md 第 N 节` 引用的设计说明：流水线结构、模块契约、阈值与验收标准。
代码实现以本文件为准；方法学细节（阈值依据、模型选择、报错处理）见 `references/`。

版本：pipeline `1.0.0`，`run.json` 的 `schema_version` 为 `1`。

---

## 1. 目标与范围

把"拿到基因型数据 → 可写进论文的结果"拆成 5 个可验证阶段（`intake → qc → popgen →
gwas → calibration`），最后汇总报告（`report`）。适用对象：动物育种/群体遗传
（猪、牛、羊、鸡等）或一般群体的芯片/重测序基因型数据；分析引擎固定为 PLINK 1.9。

不在范围内的事：混合模型（LMM）关联分析（见 README"局限"）、基因注释（需自备注释文件）、
表型本身的采集与清洗。

## 2. 阶段图与依赖

```
intake → qc → popgen → gwas → calibration → report（目录 00–05）
数据接入 → 两遍法质控 → 结构/亲缘/LD/ROH/FST → 关联 + clump → 仿真 + 独立校验 → 分段 + 总报告
```

依赖语义（详见 `references/01-pipeline-stages.md`）：

- QC 产出的 LD 剪枝集（`01_qc/data/clean_pruned.prune.in`）是后面所有群体结构计算的共同输入
  （PCA、亲缘、LD 衰减、FST 全局均值、杂合度）。
- GWAS 协变量里的 PC（PC1–PC5）必须在 QC 冻结后重算一次，不要沿用含离群样本的 PCA。
- GWAS 命令一旦冻结，校准阶段必须复用同一条命令（同过滤、同协变量、同模型）。
- ROH/FST 不进入 GWAS 模型，只用于解释、育种决策与保种评估。

各阶段模块（`scripts/plink_pipeline.py` 的 `STAGE_MODULES`）：

| 阶段 | 模块 | 输入 | 关键产物 |
|---|---|---|---|
| intake | `pipeline.intake` | 原始基因型/表型/协变量/分组 | `00_intake/manifest.json`、`preflight.json`、`data/base.*` |
| qc | `pipeline.stage_qc` | `00_intake` 的规范化数据 | `01_qc/data/clean.*`、LD 剪枝集、统计表、图 |
| popgen | `pipeline.stage_popgen` | QC 后数据 + 剪枝集 | `02_popgen/` PCA/亲缘/LD/ROH/FST |
| gwas | `pipeline.stage_gwas` | QC 后数据 + PC 协变量 | `03_gwas/tables/gwas_results.csv`、`clumps.csv`、Manhattan/QQ |
| calibration | `pipeline.calibration_runner` | GWAS 命令 + 产物 | `04_calibration/` 仿真校准 + `validation.json` |
| report | `pipeline.report` | 全部阶段的 `metrics.json` | `05_reports/final_report.md/.html`、`summary.json` |

## 3. 三条硬规则

1. **不静默丢数据**：任何过滤都要同时记录"输入数 → 输出数、原因、阈值"，写入
   `metrics.json` 与阶段 `report.md`。
2. **顺序即结果**：单次 PLINK 调用内部按固定顺序执行；跨步骤必须显式落盘中间文件
   （`--make-bed`），不能依赖一条命令里"隐含的先后关系"。
3. **可复现**：每次调用记录完整 argv、耗时、返回码（`logs/commands.jsonl`）；输入文件记
   sha256；随机过程（仿真、PCA 近似、bootstrap）固定种子并写入 metrics。

## 4. 配置与阈值

全部参数由 `scripts/pipeline/config.py` 的 `default_config()` 给出，运行开始时冻结进
`run.json`；阈值绝不写在阶段模块里。运行开始后不得根据结果反复改阈值；确需修改时先改
`run.json` 并在汇报中说明原因。

- **物种预设**（`SPECIES_PRESETS`）：人/猪/牛/羊/山羊/马/狗/小鼠/鸡/鸭的常染色体数与性染色
  体系统（XY/ZW）；`infer_species()` 从染色体编码推断物种、`--chr-set` 与 `--allow-extra-chr`。
- **QC 默认值**：`mind 0.10`、`geno 0.10`、`maf 0.01`、`hwe 1e-6`、杂合度 ±3SD、重复
  PI_HAT ≥ 0.90、近亲 PI_HAT ≥ 0.25、PCA 离群 6SD、LD 剪枝窗口 50/步长 5/r² 0.2。
- **GWAS 默认值**：协变量 = 用户协变量 + PC1–5；clump `p1 5e-8 / p2 1e-4 / r² 0.2 / kb 1000`。
- **校准默认值**：开启，10 次重复，20 个因果 SNP，h² 0.5，种子 `20261023`。
- **运行默认值**：4 线程、4000 MB、单步超时 7200 秒。
- `validate_config()` 在运行前检查阈值合法性（范围、重复/近亲阈值逻辑、线程与内存下限），
  返回 `(errors, warnings)`；有 error 时拒绝运行。

## 5. 运行目录布局

```
<out-root>/<run_id>/
├── run.json        # 冻结的配置（含 pipeline_version、schema_version、物种、阈值、版本）
├── state.json      # 各阶段状态、耗时、产物、错误
├── logs/           # pipeline.log + commands.jsonl（每次 PLINK 调用的完整 argv）
├── 00_intake/      # manifest / preflight / 清洗后的表型·协变量·分组
├── 01_qc/          # clean bfile、统计表、图、report.md
├── 02_popgen/      # PCA / 亲缘 / LD / ROH / FST
├── 03_gwas/        # gwas_results.csv、clumps.csv、Manhattan/QQ
├── 04_calibration/ # 仿真校准 + validation.json
└── 05_reports/     # final_report.md / final_report.html / summary.json
```

## 6. CLI 入口

`scripts/plink_pipeline.py`：`auto`（新建并跑全流程）、`run`（在已有 run 目录上续跑/重跑阶段）、
`status`（查看阶段状态）、`report`（重新生成汇总报告）、`selftest`（转调 `scripts/selftest.py`）。
全局选项 `--plink-exe` 指定 PLINK 可执行文件，`--version` 打印流水线版本。

执行模型：`build_context()` 组装 `PipelineContext`（run 目录、配置、`Plink` 包装、日志、状态），
`run_stage()` 按 `STAGE_MODULES` 动态导入阶段模块并调用 `module.run(ctx)`；阶段返回的字典写入
该阶段的 `metrics.json`，状态写入 `state.json`。

## 7. utils 与 Plink 包装契约

`scripts/pipeline/utils.py` 是唯一允许接触文件系统、PLINK 可执行文件与 PLINK 日志格式的模块。
阶段模块不得自行调 `subprocess` 跑 PLINK，必须经过 `Plink` 包装：

- `Plink(exe, run_dir, threads, memory_mb, logger, timeout_s)`：构造时用 `--version` 探测并记录
  PLINK 版本；每次调用记录完整 argv、耗时、返回码到 `logs/commands.jsonl`。
- `PlinkResult`：`returncode`、`elapsed_s`、`log_path`、`log_text`、`counts`、`warnings`、`errors`；
  `ok` 当且仅当返回码为 0 且无 errors。
- `PlinkError`：携带 `command`、`returncode`、`log_path`、`log_tail`、`diagnosis`（可复制的修复建议）；
  阶段失败时以诊断信息呈现，不贴原始大段日志。
- 通用工具：`write_json`（NaN/Inf 转 null，保证严格 JSON 可序列化）、`read_bim`、路径与
  `stage_dir` 映射（`STAGE_DIRS`：intake→`00_intake`，qc→`01_qc`，popgen→`02_popgen`，
  gwas→`03_gwas`，calibration→`04_calibration`，reports→`05_reports`）。

## 8. 阶段 metrics 契约与独立校验

每个阶段结束时产出 `metrics.json`，字段约定：`summary`（关键数字）、`status`
（`pass/warn/fail/error`）、`gates`（各门限的判定）、`warnings`、`errors`。
SKILL.md 的分段汇报协议只读 `metrics.json` 与 `report.md`，不读原始日志。

**门限示例**（完整表见 `references/02-qc-thresholds.md` 与 SKILL.md）：样本保留率 ≥90% PASS、
<60% FAIL；位点保留率 ≥80% PASS、<30% FAIL；λGC 0.9–1.1 PASS、0.8–1.3 WARN，之外 FAIL。

**独立校验**（`scripts/pipeline/validate.py`，由 calibration 阶段调用，产物 `validation.json`）：
不看阶段的自述数字，而是从产物文件独立重算再比对——`.fam`/`.bim`/`.bed` 字节数 vs PLINK 日志
vs `metrics.json`；`--freq`/`--missing` vs 从 `--recode A` 剂量逐位点重算；PCA 用 `.eigenvec`/
`.eigenval` 与剪枝后基因型独立重做并比对方向与特征值尺度（特征值之和 = 样本数 N）；亲缘
`.genome` PI_HAT 范围与 `.rel` 对称性；GWAS 随机抽 20 个 SNP，用 statsmodels OLS/Logit 在同样
本、同协变量下重算 β/SE，与 PLINK 结果比对（容差 1e-4）；从 P 值独立重算 λGC；HWE 比例检查。
任一项失败都要在报告中显著标注。

## 9. 报告生成

`scripts/pipeline/report.py`：

- `stage_report(run_dir, stage, metrics, config, summary_override=None)`：写
  `<run_dir>/<stage_dir>/report.md`（中文、固定章节结构）。
- `final_report(run_dir, config, state)`：汇总全部阶段的 `metrics.json`/`report.md`，产出
  `05_reports/final_report.md`、自包含的 `final_report.html`（内联 CSS + base64 图片，可直接
  用浏览器打开或打印为 PDF）与机器可读的 `summary.json`。
- `render_markdown_to_html(md_text, out_path, ...)`：小型 Markdown→HTML 渲染器（表格、代码、
  列表、引用、链接、图片）。
- 防御性原则：缺失的 metrics 字段、空表、缺失图片与跳过的阶段绝不抛异常；报告照常产出，
  并把该阶段明确标为未执行。报告中的数字必须来自产物文件，不得手写估计值；不确定就标注"需复核"。
- 报告第九章"复现命令"根据 `run.json` 重建 `auto` 命令（含输入、表型、协变量、分组、物种、
  run-id 与 out-root），并给出 `run`/`report`/`status` 的续跑命令。

## 10. 绘图

`scripts/pipeline/plots.py`：

- 接口契约：每个公开函数形如 `f(data, out_path, **kwargs) -> Path`（`plot_pca_scatter` 与
  `plot_kinship_heatmap` 的额外数据作位置参数），用 matplotlib `Agg` 后端写 PNG（DPI=150）。
- 图内文字刻意用英文：Windows 服务器上 CJK 字体不可靠，缺字会静默变成方框；中文解释放在
  Markdown/HTML 报告里。
- 函数对缺失列、空 Series、全 NaN 输入绝不抛异常，而是出一张写明原因的占位图，让阶段能继续、
  报告能照常产出。
- 坐标轴总有标题、轴标签，含义明确时加图例。

## 11. 仿真校准

`scripts/pipeline/simulate.py`（重活）+ `scripts/pipeline/calibration_runner.py`（薄 runner，
合并仿真与独立校验的报告）：

- **零模型**：表型与基因型无关（随机正态 / 随机 1/2 病例标签），与主 GWAS 完全相同的 PLINK
  命令与协变量重复 R 次，度量 λGC 与 α=0.05/1e-3/5e-8 的经验假阳性率。验收：λ 均值 ∈ [0.95, 1.05]，
  α=0.05 的 FPR ∈ [0.02, 0.08]。
- **已知 QTL**：用 truth 中的因果 SNP 按 h²=0.3/0.5 生成表型（因果效应 + 多基因背景 + 环境噪声），
  统计 5e-8 与 1e-4 阈值的检出率、效应量偏差（β_est/β_true）、SE 标定、因果 SNP 排序 AUC 与
  top100 命中率。
- **可复现**：种子来自 `config.calibration.seed`（默认 `20261023`）。

## 12. 错误处理与断点续跑

- 运行前 `preflight` 出现 `fail` 必须先修复，不许带病继续。
- PLINK 失败时 `PlinkError` 给出中文诊断与可复制修复命令（完整手册见
  `references/05-error-playbook.md`，含日志解析规则）。
- 修复输入/参数后用 `run --run-dir <dir> --stage <stage>` 续跑；禁止手工删除样本/位点掩盖问题；
  任何阈值放宽必须记录并在汇报中说明。
- 阶段状态机：`state.json` 记录每阶段 `running/pass/warn/fail/error`、开始时间与错误信息。

## 13. 可复现性

- 阈值、随机种子、软件版本（PLINK 版本、Python 版本、流水线版本）落盘后才用于结论。
- `run.json` 冻结全部配置；`logs/commands.jsonl` 记录每次 PLINK 调用的完整命令；输入文件记 sha256。
- 演示数据由 `scripts/simulate_demo.py` 按固定种子（`20261023`）生成，含 truth（因果 SNP 与效应量），
  任何人可复现同一份数据与同一组验收数字。

## 14. 验收

`scripts/selftest.py` 是一键端到端验收（默认数量性状；`--with-binary` 追加二元性状；`--check-resume`
追加断点续跑检查；`--keep-run` 保留 run 目录）：

1. 生成（或校验复用）演示数据（400 样本 × 8000 SNP，2 个亚群、10 对重复/亲缘、2% 缺失、30 个已知因果 SNP）；
2. 用主 CLI 跑完整流程；
3. 断言各阶段状态、报告/metrics 产物、独立重算校验、仿真校准（λ/FPR/power/排序）、GWAS lead SNP、
   `status` 子命令与断点续跑；
4. 打印 PASS/FAIL 清单（含定位信息：文件路径、命令、日志尾部）与总耗时。

验收阈值：零模型 λ 均值 ∈ (0.95, 1.05)；α=0.05 经验假阳性率 ∈ (0.02, 0.08)；因果 SNP 排序 AUC ≥ 0.65；
top100 命中率 ≥ 0.10；独立校验至少 12 项、最多 2 个 warn；二元性状 5e-8 power 按方法学原因不作硬门槛。
退出码：0 = 全部 PASS；1 = 有 FAIL；2 = 自检脚本自身异常/环境缺失。目标耗时 ≤ 15 分钟
（默认 demo 规模约 1–2 分钟）。

`scripts/pipeline/_report_selftest.py` 是报告层（plots.py / report.py / report.css）的独立自检：
用全伪造的 run 目录（5 个阶段的假 metrics、真实 CSV 表格与 plots.py 生成的真实 PNG）演练
`plots.*`（含退化输入）、`stage_report`、`final_report` 与 `render_markdown_to_html`。
