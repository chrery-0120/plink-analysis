# 关联分析参考（03-gwas）

> **适用场景**：在完成 QC、PCA 与亲缘评估后，对数量性状（体重、产奶量、日增重、胴体品质…）或二元性状（发病/抗性、性别、有无某表型）做全基因组关联分析，并判断结果是否可信。命令默认参数来自 `run.json` 的 `gwas` 段，可修改。

## 1. 模型选择

| 表型类型 | 推荐命令 | 输出的效应量 | 说明 |
|---|---|---|---|
| 连续数量性状 | `--linear` | 回归系数 BETA（每增加 1 个 A1 等位基因） | 加性模型；`standard-beta` 标准化表型与预测子；`genotypic`/`dominant`/`recessive` 为其他遗传模型 |
| 二元性状（0/1 或 1/2） | `--logistic` | 默认 OR；加 `beta` 输出回归系数 | 需表型编码正确；事件数少时慎用（无 Firth，见第 9 节） |
| 计数/有序/生存时间 | 不建议 PLINK 1.9 | — | 先转换（如 log、二分类）或用专门软件；转换后的解释要在报告中写明 |
| 多个相关表型 | 不建议 | — | 多变量模型（如 MANOVA/多性状混合模型）在 PLINK 之外完成 |

命令核对：

- `--linear` 与 `--logistic` 都要求表型值与协变量文件 ID 对得上；缺失协变量默认会把样本表型置缺（`--covar` 的 `keep-pheno-on-missing-cov` 可关闭该行为，但会产生协变量缺失的模型估计，**不推荐**）。
- PLINK 1.9 的 `--linear`/`--logistic` 默认只对加性剂量（A1 计数）建模；等位基因方向由 A1 决定。默认 A2 是高频等位、A1 是次等位（除非 `--keep-allele-order`）。**跨数据集比较效应量时必须统一 A1**，否则符号会反。
- 性染色体：`--xchr-model` 的 0/1/2/3 分别表示"跳过性染色体与单倍体"、"仅在 X 上加性别协变量（默认 1）"、"男性编码 0/2"、"基因型×性别交互"。`sex` 修饰符把性别加到所有染色体，`no-x-sex` 去掉 X 上的性别协变量。Y/MT 位点应单独报告，不混入常染色体多重检验。

## 2. 协变量策略

**必选**：PC1–PC5（或 PC1–PC10）、性别（若有）、明显的技术批次/芯片版本/场次。
**可选**：年龄、胎次、体重、季节、品种/群体（虚拟变量）。

原则：

1. **PC 数量**：先看 PCA 方差解释率与 λ 的变化。常用 5 或 10；多品种混合数据可用 10–20。PC 太少 → 残余结构抬高 λ；PC 太多 → 吃掉真实信号（尤其 PC 与品种/家系高度相关时）。
2. **品种/群体效应**：如果 PCA 的 PC1 已经把品种分开，再同时放品种哑变量与 PC 会共线（PLINK 会给出 VIF 警告，`--vif` 默认阈值 50）。做法二选一并在报告中说明：用 PC 代表结构，或分层分析（每个品种单独 GWAS + meta 分析）。
3. **性别**：PLINK 1.9 默认只在 X 染色体上自动加性别协变量；常染色体要显式加 `sex` 修饰符或把 sex 写进 `--covar` 文件。家禽 ZW 的性别编码与哺乳动物不同，需先确认 1/2 的含义。
4. **批次**：芯片批次与 DNA 提取批次是动物数据里最常见的隐形混杂；如果批次与群体/表型完全混杂，任何模型都无法救，只能在局限中说明。
5. **缺失值**：协变量缺失的样本默认从该次分析中剔除（表型置缺）；不要用均值插补后不记录——插补会低估 SE。
6. **共线性**：PLINK 1.9 在 `--linear` 中做 VIF 检查（默认 50），常数协变量会被自动剔除（`--no-const-covar` 是默认行为的一部分）。若某个协变量被剔除，报告中要写。

## 3. 标准命令与输出

```bat
rem 数量性状（协变量文件：FID IID PC1 ... PC5 sex batch，带表头）
plink --bfile 01_qc\data\qc_final ^
  --pheno pheno.txt --pheno-name trait ^
  --covar covar.txt --covar-name PC1-PC5,sex,batch ^
  --linear hide-covar ^
  --adjust --out 03_gwas\tables\gwas_qt

rem 二元性状（--ci 输出 OR 的 95% CI；beta 输出回归系数）
plink --bfile 01_qc\data\qc_final ^
  --pheno pheno.txt --pheno-name status ^
  --covar covar.txt --covar-name PC1-PC5,sex ^
  --logistic hide-covar beta --ci 0.95 ^
  --adjust --out 03_gwas\tables\gwas_cc

rem 条件分析（控制 lead SNP 后再找次级信号）
plink --bfile 01_qc\data\qc_final --pheno pheno.txt --covar covar.txt ^
  --covar-name PC1-PC5,sex --condition rs123456 ^
  --linear hide-covar --out 03_gwas\tables\gwas_cond
```

输出列（`--linear`）：

| 列 | 含义 |
|---|---|
| CHR / SNP / BP | 位点信息 |
| A1 | 被计数的等位基因（默认次等位） |
| TEST | `ADD`（加性）；协变量行默认隐藏（`hide-covar`），X 上会有性别项 |
| NMISS | 该位点参与回归的样本数（缺失基因型导致逐位点变化） |
| BETA | 每增加 1 个 A1 的效应（表型单位）；SE 可由 `BETA/STAT` 反推 |
| STAT | t 统计量（`--linear`）或 z 统计量（`--logistic`） |
| P | 该检验的 p 值 |

`--logistic` 默认输出 `OR`，加 `beta` 后输出 `BETA`；加 `--ci 0.95` 会输出 `SE` 与 `L95/U95`（OR 的置信区间）。`--adjust` 生成 `.adjusted` 文件，包含 `UNADJ / GC / BONF / HOLM / SIDAK_SS / SIDAK_SD / FDR_BH / FDR_BY` 列。

## 4. 多重检验

设进入检验的位点数为 M（以 `.assoc.*` 中实际检验数或 QC 后位点数为准，报告里要写清楚）：

| 阈值 | 公式 | 含义 | 适用 |
|---|---|---|---|
| Bonferroni | 0.05/M | 控制族错误率（FWER） | 保守、可复现的显著性声明 |
| 建议阈值（suggestive） | 1/M | 每个检验期望 1 个假阳性 | 筛选候选、后续验证 |
| FDR 0.05 | Benjamini–Hochberg | 控制错误发现比例 | 位点较多的探索性研究，需报告 q 值 |
| 5e-8 | 人类常用 | 约 100 万独立检验的量级 | 仅作参照；畜禽芯片 M 更小，5e-8 往往过严 |

- 5e-8 的来源：Risch & Merikangas（*Science* 1996;273:1516–1517）提出约 5e-8 的基因组显著水平；Dudbridge & Gusnanto（*Genet Epidemiol* 2008;32:227–234）用实际 LD 结构估计量级一致。它是**人类**经验值，不是自然常数。
- 畜禽芯片（50K）M≈40k–50k：Bonferroni ≈1e-6，1/M≈2e-5。用 5e-8 会显著降低 power。**推荐同时报告 Bonferroni、1/M 与 FDR 三种结果**，不要只报一个"显著"。
- 位点间存在 LD，Bonferroni 过保守；可用"有效独立检验数"（Li MX, et al. *Genet Epidemiol* 2012;36:214–224 的方法）或 LD 剪枝后的 M 做校正，但要在方法中说明。
- 多重检验的 M 不要用"全基因组位点数"和"通过 MAF 的位点数"混着写；GWAS 实际检验的是后者。

## 5. λGC（基因组控制膨胀因子）

定义：λGC = median(χ²) / 0.4549（χ²₁ 的中位数），等价于用 p 值中位数换算。PLINK 1.9 在 `--adjust` 中给出 λ（`.log` 里也会打印），`gc` 修饰符输出 GC 校正后的 p；`--lambda <值>` 可手动指定 λ 做校正。

判读（人类经验，畜禽可参考）：

| λ 范围 | 解释 | 处理 |
|---|---|---|
| 0.90–1.10 | 可接受 | 默认不做 GC 校正，报告中记录 λ |
| 1.10–1.20 | 轻度膨胀 | 检查 PC 数、亲缘、批次、表型变换；必要时加 GC 或换 LMM |
| >1.20 | 明显膨胀 | 优先怀疑群体结构/亲缘/表型-批次混杂；**不要靠 GC 掩盖问题** |
| <0.90 | 过度校正/检验统计量有偏 | 检查是否误用了 GC、表型是否离散化过度、位点是否过少 |

要点：

- λ 随样本量增大而趋近 1 的抽样波动变小：零模型下 SD(λ) ≈ 2.33/√M（M 为检验位点数；M=6000 时约 0.03）。所以"λ=1.08"在 M=6000、N=400 时可能只是噪声；在 M=500k 时则说明问题。
- case/control 的 λ 受病例/对照比例影响，不同研究之间不可直接比较。
- λ 不能区分"结构膨胀"与"多基因真实信号"：多基因性状也会有 λ>1。若要区分，用 LDSC 截距（Bulik-Sullivan BK, et al. *Nat Genet* 2015;47:291–295）或混合模型。
- 本流程默认 `gwas.genomic_control = "none"`：**只报告 λ，不自动改 p 值**；如果用户要求 GC 校正，用 `--adjust gc` 并在报告中写明"p 值经 GC 校正，λ=…"。

## 6. QQ 图判读

QQ 图 = 观测 −log10(p) 对期望 −log10(p)（均匀分布）。典型形态：

| 形态 | 可能原因 |
|---|---|
| 早期（p>1e-3）就整体上翘 | 检验统计量有偏、表型未正态化、协变量不足、GC 过度使用 |
| 大部分贴合、尾部单点/少数点上翘 | 真实关联信号（尤其与 Manhattan 峰对应） |
| 整体均匀上移 | 群体结构/亲缘/批次混杂（λ 大） |
| 整体下移（低于对角线） | 过度校正、离散表型当连续、位点极少、统计量计算异常 |
| 阶梯状/离散点 | 小样本、低 MAF、case/control 不平衡导致 p 值离散 |

报告要求：QQ 图必须标 λ 与 M，并与 Manhattan 图、clump 结果一起解释；**不能只看 λ 判定"质控通过"**。

## 7. LD clump 与 lead SNP

```bat
plink --bfile 01_qc\data\qc_final --clump 03_gwas\tables\gwas_qt.assoc.linear ^
  --clump-p1 5e-8 --clump-p2 1e-4 --clump-r2 0.2 --clump-kb 1000 ^
  --out 03_gwas\tables\clump
```

- 参数语义：`--clump-p1` 是 index 位点的 p 上限；`--clump-p2` 是允许并入 clump 的从属位点 p 上限；`--clump-r2` 是同 clump 的 LD 上限；`--clump-kb` 是物理距离半径。PLINK 1.9 的默认值较宽松（p1=1e-4、p2=0.01、r2=0.5、250 kb），本流程按 config 收紧为 5e-8 / 1e-4 / 0.2 / 1000 kb。
- clump 需要**同一套基因型数据**（`--bfile`）来计算 LD，报告文件里必须有 `SNP` 与 `P` 列（默认列名，可用 `--clump-snp-field`/`--clump-field` 指定）。
- **index SNP** 是 clump 中 p 最小的位点，**lead SNP** 通常指该 clump 的代表位点（可能与 index 相同）。报告时要区分"信号区域"与"代表位点"，不要把 clump 内所有位点都称为 lead SNP。
- 低密度芯片（如 10K–50K）里 1000 kb 半径可能把两个独立信号并成一个；小样本里 r² 估计噪声大。建议：至少用 r²=0.1 与 0.2 两档做敏感性分析；对每个 clump 报告 index SNP、位点数、p 范围、物理区间。
- `--clump-verbose` 给出 clump 内所有位点，`--clump-best` 给出最佳代理位点，写报告时用得上。

## 8. 效应量的解释与陷阱

- `--linear` 的 BETA 单位是"表型单位 / 每拷贝 A1"；若表型做了标准化（`standard-beta`），BETA 是标准差单位，可直接比较不同性状。
- `--logistic` 默认 OR：OR>1 表示 A1 增加风险（或病例概率）；case/control 的 OR 依赖基线患病率，不能当相对风险解释。
- 报告必须给出：A1/A2、A1 频率（AF）、BETA 或 OR、SE 或 95% CI、p、N、模型、协变量列表。
- **Winner's curse**：在发现集里挑出来的 top SNP，其效应量系统性高估；小样本尤其严重。没有独立验证集时，报告要写"效应量为发现集估计，存在高估可能"；有条件时用同一群体的另一半样本或独立群体验证。
- 显著 ≠ 因果：芯片位点可能只是因果变异的 LD 代理；需要注释（基因位置、eQTL、保守性）与功能证据。
- 跨数据集比较前，先统一基因组版本（bp）、等位基因方向（A1）、链方向（strand）；动物芯片不同版本间常有坐标差异，需 liftover 或按 rsID/位置匹配。

## 9. PLINK 1.9 没有混合模型：局限与替代

PLINK 1.9 的 `--linear/--logistic` 假设残差独立，**无法处理亲缘/家系/群体结构**。PC 只能吸收粗尺度结构，对"家系内相关"几乎无效。后果：SE 低估、p 值偏小、λ 偏高，家系/保种群体尤其明显。

判断是否需要混合模型：

- PI_HAT>0.25 的配对数多（例如 >5% 的样本对）、有系谱、已知半同胞家系、同一场/同一父本的大批后代；
- 加了 PC 后 λ 仍 >1.1；
- 表型在家族/群体间聚集。

替代方案（按场景选择，报告中写明"因 PLINK 1.9 无 LMM，改用 X 软件"）：

| 工具 | 方法 | 适用 |
|---|---|---|
| GCTA `--mlma` / `--fastGWA` | GRM 混合模型 / fastGWA | 中大规模、已有 PLINK 数据；`plink --make-grm-bin` 可直接产 GRM |
| BOLT-LMM | 贝叶斯混合模型 | 大样本人类/动物，速度快 |
| GEMMA | LMM（含多性状） | 小到中等样本、家系结构明确 |
| SAIGE / GMMAT | 广义混合模型（logistic） | 病例-对照、相关样本、不平衡设计 |
| regenie | 两阶段混合模型 | 生物样本库规模、多性状 |
| ssGBLUP / BLUPF90 | 单步法（系谱+基因组） | 畜禽育种值估计、已知系谱的育种群体 |

引用：Yang J, et al. *Nat Genet* 2011;43:519–525（GCTA）；Loh PR, et al. *Nat Genet* 2015;47:284–290（BOLT-LMM）；Zhou X, Stephens M. *Nat Genet* 2012;44:821–824（GEMMA）；Chen H, et al. *Am J Hum Genet* 2016;98:653–665（GMMAT）；Zhou W, et al. *Nat Genet* 2018;50:1335–1345（SAIGE）；Mbatchou J, et al. *Nat Genet* 2021;53:1097–1103（regenie）；Legarra A, et al. *J Dairy Sci* 2014;97:7998–8010（ssGBLUP 综述）。

## 参见

- PLINK 1.9 关联分析：https://www.cog-genomics.org/plink/1.9/assoc ；后处理 clump：https://www.cog-genomics.org/plink/1.9/postproc#clump ；顺序：https://www.cog-genomics.org/plink/1.9/order
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B*. 1995;57:289–300.
- Storey JD, Tibshirani R. Statistical significance for genomewide studies. *PNAS*. 2003;100:9440–9445.
- Devlin B, Roeder K. Genomic control for association studies. *Biometrics*. 1999;55:997–1004.
- Price AL, et al. Principal components analysis corrects for stratification in genome-wide association studies. *Nat Genet*. 2006;38:904–909.
- Bulik-Sullivan BK, et al. LD score regression distinguishes confounding from polygenicity. *Nat Genet*. 2015;47:291–295.
- Chang CC, et al. Second-generation PLINK. *GigaScience*. 2015;4:7.

## 10. 常见误区 / 结果可信度清单

- [ ] 表型与协变量编码正确（case/control 的 1/2 或 0/1、缺失值、数量性状分布）？PC 是否在最终 QC 集上重算，PC 数多/少各试一次？
- [ ] 性别是否加入模型？X 染色体用了哪个 `--xchr-model`？Y/MT 是否混进了常染色体多重检验？
- [ ] 是否报告了 M（实际检验位点数）、Bonferroni、1/M、FDR 三种阈值？5e-8 是否被误当成畜禽的"金标准"？
- [ ] λ 是否报告并解释？>1.1 时是否先排查结构/亲缘/批次，而不是直接 GC 校正？显著信号是否做了 clump 并写明参数与 r² 敏感性？
- [ ] 效应量是否带 A1、AF、SE/CI、N？是否提示 winner's curse？
- [ ] 家系/近亲结构是否明确标注"未用混合模型"？是否给替代软件建议？
- [ ] 是否有独立验证或至少"发现集/验证集"的设计说明？
- [ ] p 值是否在报告中保留足够位数（如 1.2e-9），而不是只写"显著"？
