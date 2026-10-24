# 全流程阶段参考（01-pipeline-stages）

> **适用场景**：拿到 bed/ped/vcf/gen/23andMe 基因型数据，需要在猪、牛、羊、鸡等畜禽的遗传育种或保种项目中跑完整 PLINK 流程；本文件给出阶段图、可直接复制的命令序列、PLINK 执行顺序语义、失败判定与回退策略。所有阈值来自 `run.json`（DESIGN 第 4 节），下列默认值均可在 config 中修改。

相关文档：`02-qc-thresholds.md`（阈值依据）、`03-gwas.md`（关联模型）、`04-popgen.md`（PCA/亲缘/LD/ROH/FST）、`05-error-playbook.md`（报错处理）、`06-calibration.md`（仿真校准）、`07-report-template.md`（报告模板）。

## 0. 三条硬规则

1. **不静默丢数据**：任何过滤都要同时记录「输入数 → 输出数、原因、阈值」，写入 `metrics.json` 与阶段 `report.md`。只写"已质控"而不写删了多少样本/位点的报告不可接受。
2. **顺序即结果**：单次 PLINK 调用内部按固定顺序执行（见第 8 节），与命令行里 flag 的书写顺序无关（PLINK 1.9 按字典序解析 flag）。跨步骤必须显式落盘中间文件（`--make-bed`），不能依赖一条命令里"隐含的先后关系"。
3. **可复现**：每次调用记录完整 argv、耗时、返回码（`logs/commands.jsonl`）；输入文件记 sha256；随机过程（仿真、PCA 近似、bootstrap）固定种子并写入 metrics。

## 1. 阶段图与依赖

```
intake → qc → popgen → gwas → calibration → report（00–05）
  数据接入 → 两遍法质控 → 结构/亲缘/LD/ROH/FST → 关联+clump → 仿真+独立校验 → 分段+总报告
  产物：00_intake/geno.*；01_qc/data/qc_final.* + prune.in；02_popgen/tables；03_gwas/tables；04_calibration/tables + validation.json；05_reports/final_report.*
```

依赖语义：

- **QC 产出的 LD 剪枝集（`prune.in`）是后面所有群体结构计算的共同输入**（PCA、亲缘、LD 衰减、FST 全局均值、杂合度）。先剪枝再算，是 PLINK 官方对 `--pca`/`--genome`/`--het`/`--ibc`/`--fst` 的一致建议（LD 会放大少数位点的权重）。
- **PCA 的协变量（PC1–PC5/10）必须在 QC 冻结后重算一次**，不要沿用包含离群样本的 PCA：离群个体会改变主成分方向。
- **GWAS 命令一旦冻结，校准阶段必须复用同一条命令**（同过滤、同协变量、同模型），否则校准测的不是实际使用的模型。
- ROH/FST 不进入 GWAS 模型，只用于解释、育种决策与保种评估；不要为了"好看的 FST"而改 QC 阈值。

## 2. 数据接入（intake）

**输入**：`.bed/.bim/.fam`、`.ped/.map`、`.vcf/.vcf.gz`（可 bcf）、Oxford `.gen/.bgen`、23andMe 文本。
**输出**：统一二进制 `00_intake/geno.*`、`manifest.json`（文件/大小/sha256）、`preflight.json`、`metrics.json`（第 8 节 intake 契约）。

### 2.1 转换命令

```bat
rem bed/bim/fam 只校验；ped/map 与 VCF 转换时统一为二进制
plink --bfile raw --freq --out 00_intake\preflight
plink --file raw --make-bed --out 00_intake\geno
rem VCF：转换时按 2.2 加 --chr-set；有非标准 contig 再加 --allow-extra-chr 0
plink --vcf raw.vcf.gz --double-id --make-bed --out 00_intake\geno
plink --vcf raw.vcf.gz --double-id --allow-extra-chr 0 --make-bed --out 00_intake\geno
rem 23andMe（缺失 FID 用 --double-id；文本先清洗空格）
plink --23file raw.txt --double-id --make-bed --out 00_intake\geno
```

### 2.2 物种染色体集（必须先定）

| 物种 | 推荐 flag | 等价含义 | X 的数值编码 |
|---|---|---|---|
| 人 | 默认 | 22 常染色体 + X/Y/XY/MT | 23 |
| 猪（Sus scrofa） | `--chr-set 18` | 18 常染色体 | 19 |
| 牛（Bos taurus） | `--cow` | `--chr-set 29 no-xy` | 30 |
| 绵羊（Ovis aries） | `--sheep` | `--chr-set 26 no-xy no-mt` | 27 |
| 山羊（Capra hircus） | `--chr-set 29 no-xy` | 29 常染色体 | 30 |
| 马 | `--horse` | `--chr-set 31 no-xy no-mt` | 32 |
| 鸡（Gallus gallus） | `--chr-set 38 no-xy`（见下） | 38 常染色体，ZW 系统需特殊处理 | 39（Z 按 X 处理时） |

来源：PLINK 1.9 输入页 “Nonstandard chromosome IDs”（https://www.cog-genomics.org/plink/1.9/input#allow_extra_chr）。规则是：**n 对常染色体时，X=n+1、Y=n+2、XY=n+3、MT=n+4**；`no-x/no-y/no-xy/no-mt` 是"该物种没有这条染色体"的合法性检查。

- `--allow-extra-chr` 与 `--chr-set` **不是二选一**：前者是"允许任意未识别编码"（如 `chrUn`、`scaffold12`），后者是"重新定义标准染色体集"。混用会把物种的 X 染色体当成普通编号，性别核查和 X 关联分析随之失效。正确顺序：先 `--chr-set`，确实还有 contig 才加 `--allow-extra-chr`。
- **家禽 ZW**：PLINK 1.9 没有原生 ZW 支持。实践做法是把 Z 编码为 X（`--chr-set 38 no-xy` 后 X=39）并把 W 作为额外染色体，但 PLINK 的 X 模型假设"雄性半合子"，与 ZW（雌性半合子）相反，因此 **`--check-sex` 与 X 模型关联在家禽中不可直接解释**；主分析建议只用常染色体 1–38，Z 单独用 PLINK 2 或专门的 ZW 工具处理，并在报告中写明局限。

### 2.3 preflight 检查（失败即停，不进入 QC）

1. 文件齐全：三件套存在且非空；`.fam` 行数=样本数、`.bim` 行数=位点数；`.bed` 大小 = 3 + ceil(样本数/4)×位点数。
2. 编码合法：染色体集合属于当前 `--chr-set`（与 config `species.chr_set` 比对）；等位基因只含 A/C/G/T、0/1 或已声明 indel（`<DEL>` 等要先转换）。
3. ID 可对齐：样本 ID 唯一、位点 ID 唯一或可用 `--set-missing-var-ids` 补齐；表型/协变量 ID 100% 对上基因型样本（见 `05-error-playbook.md`）。

### 2.4 失败与回退

| 判定 | 触发 | 回退 |
|---|---|---|
| FAIL | 文件缺失/为 0 字节/.bed 魔数错误 | 停止，输出诊断；提示重新下载或转换，不写入 qc |
| FAIL | 染色体编码非法且无法用 `--chr-set` 解释 | 生成"编码清单"表，交用户确认物种；必要时 `--allow-extra-chr` 仅保留常染色体分析 |
| WARN | 样本 <200 或位点 <5000 | 继续，但在报告中标明"小样本/低密度芯片，PCA、ROH、λ 的稳定性下降" |
| WARN | 表型缺失 >10% | 继续，报告缺失数与后续有效样本数；不自动插补 |

## 3. 质控（qc）

原则：**两遍法**。先描述统计（不删），再按 PLINK 顺序过滤。原因见第 8 节：同一条命令里 `--mind` 在 `--geno` 之前执行，样本缺失率会把"稍后会被删掉的坏位点"算进去，导致过度剔除样本。

### 3.1 第一遍：描述统计（只读）

```bat
plink --bfile 00_intake\geno ^
  --missing --freq --hardy --het ^
  --out 01_qc\tables\raw_stats
rem 输出 .imiss/.lmiss（缺失）、.frq（MAF）、.hwe（ALL/AFF/UNAFF）、.het（F 值）
rem --hardy 在 case/control 下会给出 ALL/AFF/UNAFF 三套检验；--het 默认只算创始人
```

### 3.2 位点过滤 → 样本过滤 → 杂合度（分步落盘）

```bat
rem (1) 位点：缺失率 → HWE → MAF（PLINK 内部顺序，可在一条命令里）
plink --bfile 00_intake\geno --geno 0.10 --hwe 1e-6 --maf 0.01 ^
  --make-bed --out 01_qc\data\step1_snp
rem (2) 样本缺失率（在已过滤位点集上重算）
plink --bfile 01_qc\data\step1_snp --mind 0.10 ^
  --make-bed --out 01_qc\data\step2_sample
rem (3) 常染色体 LD 剪枝（供亲缘/PCA/杂合度/HWE 复核共用）
plink --bfile 01_qc\data\step2_sample --chr 1-18 ^
  --indep-pairwise 50 5 0.2 --out 01_qc\data\prune
rem (4) 杂合度离群（±3 SD，在剪枝集上；Python 算 mean±3SD 后 --remove）
plink --bfile 01_qc\data\step2_sample --extract 01_qc\data\prune.prune.in ^
  --het --out 01_qc\tables\het
rem (5) 亲缘/重复：PI_HAT>0.90 重复；0.25–0.90 近亲（是否剔除由 config.remove_related 决定）
plink --bfile 01_qc\data\step2_sample --extract 01_qc\data\prune.prune.in ^
  --genome --min 0.2 --out 01_qc\tables\ibd
rem (6) 性别核查：先剪枝 X，再检查
plink --bfile 01_qc\data\step2_sample --chr 19 --indep-pairwise 50 5 0.2 ^
  --out 01_qc\data\prune_x
plink --bfile 01_qc\data\step2_sample --extract 01_qc\data\prune_x.prune.in ^
  --check-sex 0.2 0.8 --out 01_qc\tables\sexcheck
rem (7) PCA（用 header 让 .eigenvec 可直接当协变量）；(8) 应用剔除清单冻结 QC
plink --bfile 01_qc\data\step2_sample --extract 01_qc\data\prune.prune.in ^
  --pca 20 header --out 01_qc\tables\pca
plink --bfile 01_qc\data\step2_sample --remove 01_qc\data\remove_samples.txt ^
  --make-bed --out 01_qc\data\qc_final
```

要点与依据：

- `--check-sex` 默认女性 F<0.2、男性 F>0.8，**要求 X 的 PAR 已用 `--split-x` 拆开**；但 `--split-x` 只内置人类 build 代码（b36/b37/b38），畜禽必须给自定义 PAR 区间，或接受 PAR 位点带来的少量偏差。官方文档明确警告：不剪枝时 `--check-sex` 会给出错误性别判定（https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex）。
- 重复样本处理：**保留缺失率更低、表型/协变量更完整的一条**；若二者数据都有价值，可先合并再重新基因型，而不是简单删一条。
- `remove_related` 与 `remove_pca_outliers` 在 config 中默认 `false`：科研惯例是"亲属关系用混合模型/GRM 处理，而不是直接删"；育种群体本身就由家系构成，盲目剔除会引入选择偏倚。仅当用户明确要求"独立样本 GWAS"时才打开。
- QC 冻结后要记录：`n_samples_in/out`、`n_variants_in/out`、各原因剔除数、缺失率均值、MAF 谱、HWE 极值数、PI_HAT 分布、PCA 离群数（第 8 节 qc metrics）。

### 3.3 判定（建议阈值，可在 config 调整）

| gate | PASS | WARN | FAIL |
|---|---|---|---|
| 样本保留率 | ≥90% | 70–90% | <70% 或剩余 <50 |
| 位点保留率 | ≥70% | 50–70% | <50% 或剩余 <1000（芯片项目） |
| 单个样本缺失率最大值 | ≤mind | — | >2×mind |
| HWE 极值比例（p<1e-6） | <1% | 1–5% | >5%（提示分型/混合群体问题） |
| 性别不一致 | 0 | 1–3 | >3 或比例 >2% |
| PI_HAT 最大 | ≤0.90 | >0.90 但已处理 | 未处理 |

## 4. 群体结构 / 亲缘 / LD / ROH / FST（popgen）

**输入**：`01_qc/data/qc_final`、`prune.prune.in`、分组文件（`--within`）。
**输出**：第 8 节 popgen metrics：PCA、relatedness、LD bins、ROH/FROH、FST、diversity。

```bat
rem PCA（剪枝集；header 便于直接当协变量）
plink --bfile 01_qc\data\qc_final --extract 01_qc\data\prune.prune.in ^
  --pca 20 header --out 02_popgen\tables\pca
rem 亲缘 PI_HAT；GRM 供混合模型（GCTA/SAIGE）
plink --bfile 01_qc\data\qc_final --extract 01_qc\data\prune.prune.in ^
  --genome --min 0.2 --out 02_popgen\tables\ibd
plink --bfile 01_qc\data\qc_final --extract 01_qc\data\prune.prune.in ^
  --make-grm-bin --out 02_popgen\data\grm
rem LD 衰减：成对 r²，Python 按 50 kb 分箱（PLINK 1.9 --r2 无 bins 修饰符）
plink --bfile 01_qc\data\qc_final --extract 01_qc\data\prune.prune.in ^
  --r2 --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0 ^
  --out 02_popgen\tables\ld_pairs
rem ROH/FROH（config 默认 ≥50 位点、≥1000 kb、密度 1/50kb、gap 1000kb；随物种调整）
plink --bfile 01_qc\data\qc_final --chr 1-18 ^
  --homozyg --homozyg-snp 50 --homozyg-kb 1000 --homozyg-density 50 --homozyg-gap 1000 ^
  --out 02_popgen\tables\roh
rem FST（Weir & Cockerham 1984；需要 --within 分组文件 FID IID group）
plink --bfile 01_qc\data\qc_final --chr 1-18 ^
  --within groups.txt --fst --out 02_popgen\tables\fst
rem 杂合度 / 近交系数
plink --bfile 01_qc\data\qc_final --extract 01_qc\data\prune.prune.in ^
  --het --ibc --out 02_popgen\tables\inbreed
```

解释与坑见 `04-popgen.md`。三个必须写进报告的元信息：剪枝位点数、染色体范围（是否含性染色体）、物种基因组总长（FROH 的分母）。

## 5. 关联分析（gwas）

**输入**：`qc_final`、表型文件、协变量文件（PC1–PC5/10、性别、批次、品种/群体）。
**输出**：`.assoc.linear` / `.assoc.logistic`、`.adjusted`、clump 表、lead SNP 表、Manhattan/QQ 图。

```bat
rem 数量性状
plink --bfile 01_qc\data\qc_final ^
  --pheno pheno.txt --pheno-name trait ^
  --covar covar.txt --covar-name PC1-PC5,sex,batch ^
  --linear hide-covar ^
  --adjust --out 03_gwas\tables\gwas
rem 二元性状（--ci 只对 --logistic 输出 OR 的置信区间）
plink --bfile 01_qc\data\qc_final ^
  --pheno pheno.txt --pheno-name status ^
  --covar covar.txt --covar-name PC1-PC5,sex ^
  --logistic hide-covar beta --ci 0.95 ^
  --adjust --out 03_gwas\tables\gwas
rem LD clump（p1 为 index 阈值，p2 为从属位点阈值，r2/kb 定义同一信号）
plink --bfile 01_qc\data\qc_final ^
  --clump 03_gwas\tables\gwas.assoc.linear ^
  --clump-p1 5e-8 --clump-p2 1e-4 --clump-r2 0.2 --clump-kb 1000 ^
  --out 03_gwas\tables\clump
```

注意：

- 表型文件缺失值编码必须与 config 一致（默认 `-9/NA/0`）；`NA` 会被 PLINK 视为缺失，但 `0` 在 case/control 里是缺失、在数量性状里是合法值，需在 intake 阶段判定表型类型。
- PLINK 1.9 **没有混合模型**（无 LMM/GRM 关联）。家系结构明显的育种群体，`--linear/--logistic` 会低估 SE、抬高 λ；此时报告必须写明"未校正亲缘"，并给 GCTA/SAIGE/GMMAT/regenie 的替代建议（`03-gwas.md`）。
- `--adjust` 同时给出 Bonferroni、Holm、Sidak、FDR_BH、FDR_BY 与 GC 校正列（加 `gc` 修饰符）。默认流程只报告 λ 与原始 p，不自动做 GC 校正（`gwas.genomic_control = "none"`）。
- 若用 `.eigenvec` 当协变量：**必须带 header**（`--pca ... header`），否则 PLINK 读不到列名；更稳妥的做法是脚本生成 `FID IID PC1...` 的协变量文件。

## 6. 仿真校准（calibration）

**输入**：QC 冻结数据（建议取 20k 位点子集，加速）、GWAS 冻结命令。
**输出**：零模型 λ 分布与经验假阳性率、已知 QTL 的 power/效应偏差/排名 AUC、`validation.json`。

流程要点（详细设计见 `06-calibration.md`）：

1. 零模型：表型与基因型独立（随机正态或置换真实表型），重复 R 次跑**同一条 GWAS 命令**，记录每次 λ 与最小 p，统计 α=0.05/1e-3/5e-8 的经验假阳性率。
2. 已知 QTL：h2=0.3/0.5、20 个因果 SNP（MAF>0.05、彼此近似独立），按 `y = Σβg + ε` 生成表型，统计 power、`β_hat/β_true`、因果 SNP 排名 AUC。
3. 独立重算：MAF、缺失率、λ、20 个 SNP 的 beta/SE 与 statsmodels 对照（容差见验收标准）；验收：λ 均值 ∈[0.97,1.03]（DESIGN 第 14 节）、power 与解析预期相符、validation 全 PASS。

**校准不是可选项**：它证明"这条流水线在你的数据规模下没有系统性错误"，也是报告里唯一能支撑"结果可信"的定量证据。

## 7. 报告（report）

每阶段固定 7 段：目的 → 方法与参数 → 结果（表/图）→ 质控判定（PASS/WARN/FAIL）→ 异常与处理 → 对后续/结论的影响 → 产物清单。最终报告结构见 `07-report-template.md` 与 DESIGN 第 9 节。

CLI 对 agent 友好的输出约定（DESIGN 第 11 节）：每阶段末尾打印单行 `SUMMARY_JSON: {...}`，再打印 `SUMMARY_MD: ...`；退出码 0=全通过、1=有 WARN、2=有 FAIL/ERROR。

## 8. PLINK 1.9 顺序语义（写报告和排错都要用）

官方 order of operations（https://www.cog-genomics.org/plink/1.9/order）中与本流程相关的固定顺序（节选）：

| 阶段 | 固定顺序 |
|---|---|
| 初始化 | 定义染色体集（`--chr-set`/`--cow`…，未指定按人类）→ `--out` → 记录 `--chr/--autosome` 过滤器 |
| 读入与改名 | 读数据 → `--set-missing-var-ids` → 单区间选择 → `--snps-only` → 更新位点/等位基因/`--flip` |
| 位点选择 | `--extract/--exclude`（**在 `--keep/--remove` 之前**）→ 属性/质量分 → `--thin` |
| 样本选择 | `--update-ids/--update-sex` → `--keep/--remove` → `--filter` |
| 表型处理 | 缺失性别表型置缺（除非 `--allow-no-sex`）→ `--prune`（丢缺失表型）→ 性别/表型/创始人过滤 → `--thin-indiv` |
| 样本缺失 | **`--mind`**（此时位点尚未做 `--geno`）→ `--within` 定义分组 |
| 协变量 | **载入 `--covar`** → 去常量协变量 |
| 位点统计与过滤 | **`--missing` → `--geno` → `--hardy/--hwe` → `--maf/--max-maf/--mac` → `--bp-space`** |
| 近交/关系 | `--make-rel/--make-grm-*` → `--rel-cutoff` → `--ibc` |
| 主成分 | **`--pca`**（在 `--indep` 之前！） |
| 性别/输出 | `--check-sex/--impute-sex` → 写 .cov / 二进制 / 文本 |
| 群体参数 | **`--homozyg`（ROH）→ `--indep*`（LD 剪枝）→ `--flip-scan` → `--r/--r2`** |
| 关系与 F 统计 | `--distance → --read-dists` → **`--genome`（IBD）→ `--het` → `--fst`** |
| 关联 | **`--assoc/--linear/--logistic` → `--clump` → 退出** |

由此得到四条实用推论：

- **`--mind` 在 `--geno` 之前**：所以样本缺失率要在"位点过滤后的文件"上重算（第 3.2 节分步的原因）。
- **`--extract/--exclude` 在 `--keep/--remove` 之前**：可以同时指定，互不影响；但 `--extract` 列表里不存在的 ID 默认报错，脚本要先确认列表与数据一致（或容忍缺失）。
- **`--pca` 在 `--indep*` 之前**：不要写"`--indep-pairwise ... --pca`"企图一次完成；PCA 会用到未剪枝的位点。必须两步：先 `--indep-pairwise` 生成 `prune.in`，再 `--extract prune.in --pca`。
- **`--het`/`--fst` 在 `--genome` 之后、关联分析之前**：因此 QC 阶段的 `--het`/`--genome` 与 GWAS 是两次独立调用；同一次调用里再指定 `--maf` 等过滤器，会作用于该次调用的所有分析。

## 9. 最小可用流程 vs 完整流程

### 9.1 最小可用（30 分钟内出结果，适合快速评估）

```bat
plink --vcf raw.vcf.gz --double-id --make-bed --out geno
plink --bfile geno --geno 0.10 --hwe 1e-6 --maf 0.01 --mind 0.10 --make-bed --out qc
plink --bfile qc --indep-pairwise 50 5 0.2 --out prune
plink --bfile qc --extract prune.prune.in --pca 10 header --out pca
plink --bfile qc --pheno pheno.txt --covar pca.eigenvec --covar-name PC1-PC5 --linear hide-covar --out gwas
```

### 9.2 完整流程（科研交付）

在最小流程基础上补齐：preflight 全项、两遍法 QC、性别核查、亲缘与重复、PCA 离群、ROH/FROH、FST、LD 衰减、clump、λGC/FDR、仿真校准、独立重算校验、分段报告 + 最终报告。完整流程的每次过滤都必须可回溯到具体参数与计数。

## 参见

- PLINK 1.9 文档：https://www.cog-genomics.org/plink/1.9/ （Input / Filter / Basic stats / LD / IBD / Stratification / Assoc / Postprocessing / Order of operations）
- Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*. 2015;4:7.
- Purcell S, et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. *Am J Hum Genet*. 2007;81(3):559–575.
- Anderson CA, et al. Data quality control in genetic case-control association studies. *Nat Protoc*. 2010;5(9):1564–1573.

## 10. 常见误区 / 检查清单

- [ ] 有没有把 `--allow-extra-chr` 当成 `--chr-set` 用，导致 X 染色体被当普通染色体、`--check-sex` 结果无效？
- [ ] 有没有在一条命令里同时写 `--indep-pairwise` 和 `--pca`，以为 PCA 用了剪枝集？（PLINK 顺序里 `--pca` 在前）
- [ ] 样本缺失率是不是在 `--geno` 之后的文件上算的？（同一条命令里 `--mind` 先执行）
- [ ] 家禽数据是否误用 `--check-sex`/X 模型？报告中是否写明 ZW 局限？协变量 `.eigenvec` 是否带 header 且在最终 QC 集上重算？
- [ ] 重复/近亲样本是删除还是保留？删除理由是否与"家系群体/保种群体"的科研目标一致？
- [ ] 每次 PLINK 调用的过滤参数、输入输出计数、sha256 是否都进了 `commands.jsonl` 与 metrics？
- [ ] ROH 的 `--homozyg-snp/kb/density/gap` 是否随芯片密度和物种调整？FST 分组 ID 是否一致、全局均值是否在剪枝集上算？
- [ ] 校准是否复用了"同一条冻结 GWAS 命令"？有没有因为校准而偷偷改过过滤器？
