# 群体遗传与亲缘参考（04-popgen）

> **适用场景**：在 QC 之后描述群体结构、亲缘关系、连锁不平衡、近交与分化，用于 GWAS 协变量、育种群体评估、保种群体多样性监测。命令默认参数来自 `run.json` 的 `popgen` 段，可修改。所有指标都应在 **LD 剪枝后的常染色体位点集**上计算（除非明确说明）。

## 1. PCA（主成分分析）

### 1.1 原理与命令

PLINK 1.9 的 `--pca` 先计算方差标准化的关系矩阵（variance-standardized relationship matrix），再取特征向量作为样本得分（样本数 × PC）。官方建议在近似连锁平衡的位点集上做。

```bat
rem 第一步：剪枝（常染色体）
plink --bfile 01_qc\data\qc_final --chr 1-18 ^
  --indep-pairwise 50 5 0.2 --out 02_popgen\data\prune

rem 第二步：在剪枝集上提取主成分（header 便于直接当协变量）
plink --bfile 01_qc\data\qc_final --extract 02_popgen\data\prune.prune.in ^
  --pca 20 header --out 02_popgen\tables\pca
```

输出：`.eigenvec`（FID IID PC1…PCk）、`.eigenval`（特征值，只写出请求的 k 个）。

### 1.2 参数

- **剪枝参数**：`50 5 0.2` 是 PLINK 文档的经典起点（窗口 50 个位点、步长 5、r²>0.2 剔除）。畜禽芯片 LD 跨度长，建议用 `100 10 0.2` 或 `50 5 0.1` 复核；高密度芯片（>500K）用 `--indep-pairwise 50kb 5 0.2`（窗口按物理距离）。
- **PC 个数**：`--pca 20` 足够做筛查；GWAS 协变量取前 5 或 10。PC 个数的选择依据是"再加 PC 后 λ 不再下降、且不吞掉候选信号"。
- **方差解释率**：λ_i / Σλ。PLINK 只输出请求的 k 个特征值，Σλ 用样本数 n 近似（方差标准化关系矩阵对角线≈1，迹≈n），因此报告写"前 k 个 PC 解释率（近似）"并注明计算方式。

### 1.3 离群判定

- 简单规则：样本在 PC1–PC10 空间到群体中心的欧氏距离（按各 PC 标准差标准化）> `qc.pca_outlier_sd`（默认 6 SD）。
- 稳健规则：用中位数与 MAD（中位数绝对偏差）代替均值与 SD；或 Mahalanobis 距离 + 卡方阈值。
- **不要只看 PC1–PC2 的图**：有些离群只在 PC3–PC10 显现；建议同时输出 PC1–2、PC3–4、PC5–6 三张图。
- 群体结构不是错误：多品种/多场次数据在 PC 上的分离是生物学信号；只有当离群样本**同时**具有异常缺失率、杂合度、性别或表型时才按数据错误处理。

### 1.4 解释与坑

- PC1 常对应品种/地理/场次；在动物数据中也可能对应芯片批次（看 PC 与批次是否完全重合）。
- 小样本（N<100）的 PC 不稳定：重抽样（bootstrap）后 PC 方向可能翻转；报告结论要保守。
- 芯片 ascertainment 偏差：商用芯片位点是按商业品种设计的，地方品种的 PC 可能被人为夸大或压缩。
- 新样本投影：PLINK 1.9 可用 `--pca-clusters`/`--pca-cluster-names` 定义参考群体后再投影外部样本（help 中有说明）；跨芯片/跨版本投影前必须统一位点与等位基因方向，否则 PC 无意义。

## 2. IBD / 亲缘（PI_HAT）

### 2.1 命令与输出

```bat
plink --bfile 01_qc\data\qc_final --extract 02_popgen\data\prune.prune.in ^
  --genome --min 0.2 --out 02_popgen\tables\ibd

rem 需要 GCTA 格式 GRM（供 GCTA/混合模型使用）时：
plink --bfile 01_qc\data\qc_final --extract 02_popgen\data\prune.prune.in ^
  --make-grm-bin --out 02_popgen\data\grm
```

`.genome` 关键列：`PI_HAT`（= Z2 + 0.5·Z1）、`Z0/Z1/Z2`（P(IBD=0/1/2)）、`IBS0/IBS1/IBS2`、`SNP_COUNT`。`--min` 控制报告下限（常用 0.2），`--max` 可限制上限；`--rel-check` 只保留同 FID 的配对（家系内验证）。

### 2.2 PI_HAT 判读

| PI_HAT | 关系（理想期望） | 处理 |
|---|---|---|
| >0.90 | 同一个体/单卵双生/重复样本 | 删除一条，保留质量更高者 |
| 0.40–0.60 | 亲子、全同胞 | 保留（用混合模型处理）；若要求独立样本则每对删一条 |
| 0.20–0.30 | 半同胞、祖孙、叔侄 | 同上；阈值边界附近不要机械判定 |
| 0.05–0.20 | 较远亲属/同群体背景 | 保留，作为群体结构信息 |
| <0.05 | 无关（大群体） | — |

注意：

- PI_HAT 是**矩估计**，受 MAF、LD、样本量影响；位点数 <1000 时方差很大，阈值附近的判定不可靠。先在剪枝集上算，并检查 `SNP_COUNT` 是否一致。
- 封闭保种群体即使"无关"个体也可能有 0.05–0.15 的背景 PI_HAT；阈值应按本群体 PI_HAT 分布的低谷位置定，而不是照搬人类 0.25。
- 育种群体中 PI_HAT 不是系谱亲缘系数：半同胞家系在芯片上可能落在 0.2–0.35，全同胞可到 0.5–0.6；有系谱时应同时报告系谱亲缘并检查一致性（不一致是系谱错误/采样错误的线索）。
- `--make-grm-bin` 输出的 GRM 与 PI_HAT 尺度不同（GRM 非对角元 ≈ 2×亲缘系数），不要混用阈值。KING-robust 亲缘估计在 PLINK 2 中提供（`--make-king`），PLINK 1.9 没有。

### 2.3 家系/育种中的用途

- 重复样本识别（同个体多次采样、样本标签错误）；
- 亲子关系验证（`--genome --rel-check` 配合 FID 分组）；
- 选配设计：避免高亲缘个体交配；
- 近交监测：与 FROH、`--het` 联合解释；
- 混合模型：把 GRM 交给 GCTA/SAIGE 等软件。

## 3. LD 衰减

### 3.1 命令

PLINK 1.9 的 `--r2` **没有按距离分箱的 `bins` 修饰符**（只有 matrix/inter-chr/形态修饰符），因此做法是输出成对 r²，再由脚本按距离分箱：

```bat
plink --bfile 01_qc\data\qc_final --chr 1-18 ^
  --r2 --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0 ^
  --out 02_popgen\tables\ld_pairs
rem 由 Python 读取 .ld，按 50 kb 分箱计算 mean/median r² 与配对数
```

参数：`--ld-window`（按位点计数，默认约 10）、`--ld-window-kb`（默认 1000）、`--ld-window-r2`（默认 0.2，报告阈值）。`--r2` 与 `--r` 的区别是 r² 与 r；`dprime` 可加 D′。

### 3.2 解释

- LD 衰减曲线回答：**多远的距离内位点还高度相关**。畜禽中常见量级：商业猪/牛品种 r²≈0.2 的距离在几十到几百 kb；地方品种通常更短；鸡的 LD 衰减更快（基因组小、有效群体历史不同）。
- 用途：决定 GWAS 分辨率、标记密度、基因组选择（GS）的标记数、LD 剪枝窗口、ROH 阈值。
- 小样本偏差：无 LD 位点的 r² 基线约 1/N 量级（N 为个体数），N=100 时基线可达 0.01 左右；小样本的 LD 曲线整体抬高，比较不同研究前要核对 N。
- 芯片 ascertainment：商用芯片偏向常见变异与商业品种，LD 曲线不能直接外推到地方品种或全基因组。
- 用 r² 还是 D′：r² 更适合"标记能否替代因果位点"（GWAS/GS），D′ 更适合单倍型/重组历史；报告里要写清用哪个。
- 与剪枝的关系：LD 衰减的窗口（如 1000 kb）与剪枝窗口不是一回事；剪枝是为了去掉冗余，衰减是为了描述结构。

## 4. ROH 与 FROH

### 4.1 原理与参数

PLINK 1.9 `--homozyg` 用**扫描窗口**算法（不是 HMM）：窗口内纯合比例达到阈值即"命中"，再合并成 ROH。官方默认：

| 参数 | 默认 | 含义 |
|---|---|---|
| `--homozyg-snp` | 100 | ROH 至少包含的位点数 |
| `--homozyg-kb` | 1000 | ROH 最小长度（kb） |
| `--homozyg-density` | 50 | 平均每 50 kb 至少 1 个位点 |
| `--homozyg-gap` | 1000 | 相邻位点间隔超过 1000 kb 则不连成同一 ROH |
| `--homozyg-het` | 无限制 | ROH 内允许的杂合数（可设 0） |
| `--homozyg-window-snp` | 50 | 扫描窗口位点数 |
| `--homozyg-window-het` / `-missing` | 1 / 5 | 窗口内允许的杂合/缺失数 |
| `--homozyg-window-threshold` | 0.05 | 窗口命中率下限 |

本流程 config 默认：`snp=50, kb=1000, density=50, gap=1000, het_missing=0.0`（比 PLINK 默认更适合中低密度芯片）。

### 4.2 畜禽芯片的阈值调整

| 场景 | 建议 |
|---|---|
| 50K 芯片（猪/牛/羊） | `--homozyg-snp 30–50 --homozyg-kb 500–1000 --homozyg-density 100 --homozyg-gap 1000` |
| 高密度芯片（>500K）或 WGS | `--homozyg-snp 50–100 --homozyg-kb 1000–2000`；WGS 可提高 density 要求 |
| 低密度芯片（<20K） | 降低 snp 到 20–30、kb 到 300–500，并**先 imputation** 再算；否则 FROH 系统性低估 |
| 小家系/近交群体 | 允许 ROH 内少量杂合（`--homozyg-het 1–2`），或检查是否把真实重组片段切碎 |
| 家禽（38 常染色体、基因组小） | 按物理距离调整 kb 阈值；微染色体（<10 Mb）上的 ROH 要单独看，容易因位点稀少产生假 ROH |

### 4.3 FROH 的计算与解释

- 定义：**FROH = Σ ROH 长度 / 常染色体基因组总长**。分母必须明确：是用芯片覆盖的染色体跨度求和、还是用参考基因组组装长度（如猪 Sscrofa11.1、牛 ARS-UCD1.2、绵羊 Oar_rambouillet_v1.0、鸡 GRCg6a）。**报告里必须写分母来源**，否则 FROH 不可比。
- 与系谱近交系数比较：FROH 只捕捉"足够长、可被当前芯片检出"的片段，通常与系谱 F 相关性中等（文献常见 r≈0.4–0.7），且反映的是较久远的近交（取决于 ROH 长度下限：1000 kb ≈ 数十代以前的共同祖先）。
- 与 `--het` 的 F 比较：F 是全基因组纯合度的矩估计，受 MAF 影响；FROH 对片段的分布更敏感。两者结合看更可靠。
- 拆分：可按 ROH 长度分箱（如 1–2 Mb、2–4 Mb、4–8 Mb、>8 Mb），长片段代表近期近交，短片段代表古老近交。report 中给出各长度段的占比。
- 选择信号：ROH 频率在群体中异常高的区域（ROH island）可能受选择，但也可能是低重组/低基因密度/组装缺口造成的假象，需与 FST、iHS 等独立证据结合。

## 5. FST（群体分化）

### 5.1 命令与输出

```bat
rem groups.txt: FID IID group（空格或 tab 分隔，无表头）
plink --bfile 01_qc\data\qc_final --chr 1-18 ^
  --within groups.txt --fst --out 02_popgen\tables\fst
rem 只有两个群体时，也可用表型当分组：
plink --bfile 01_qc\data\qc_final --fst case-control --out 02_popgen\tables\fst2
```

PLINK 1.9 用 **Weir & Cockerham (1984)** 方法逐位点估计（仅常染色体二倍体位点），日志中给出**原始均值与加权均值**（加权均值按位点方差/样本量加权，通常更稳健）。输出 `.fst` 文件含位点级 FST。

### 5.2 解释

- Wright 的经验分级（*Evolution and the Genetics of Populations*, 1978）：<0.05 分化很小，0.05–0.15 中等，0.15–0.25 大，>0.25 很大。畜禽品种间 FST 常见 0.05–0.30；同一品种内不同场/家系通常 <0.05。
- 全局均值要在**剪枝集**上算（官方建议）：LD 会让少数区域重复计权。
- 逐位点 FST 噪声大（尤其小样本/低 MAF），用于扫描选择信号时要设最小样本量与 MAF 下限，并用窗口平均（如 100 kb 窗口）平滑；PLINK 1.9 没有窗口 FST，需脚本或 VCFtools。
- 多群体时：先算总体 FST，再做两两 FST 矩阵（循环 `--within` 两两分组）；报告矩阵比单个数字更有信息量。
- 偏差来源：样本量不等（W&C 已部分校正）、芯片 ascertainment（位点按商业品种设计会低估地方品种分化）、群体混合（FST 对混合比例敏感）、小样本（偏大）。
- 近交与分化的区别：FST 描述群体间等位基因频率差异；FROH/`--het` 描述个体纯合度。保种项目两者都要报。

## 6. 杂合度与近交系数

```bat
plink --bfile 01_qc\data\qc_final --extract 02_popgen\data\prune.prune.in ^
  --het --ibc --out 02_popgen\tables\inbreed
```

- `--het` 的 `F = (O(HOM) − E(HOM)) / (N − E(HOM))`，method-of-moments；F>0 纯合过多（近交），F<0 杂合过多（污染/混合群体/重复样本）。
- `--ibc` 输出 GCTA 的三个估计（Yang et al. 2011）：Fhat1 = 方差标准化关系减 1，Fhat2 类似 `--het`，Fhat3 基于配子间相关。三者接近时结论更稳；差异大时先查 MAF 与群体结构。
- 观测杂合度 Ho、期望杂合度 He、FIS = 1 − Ho/He 可按群体分层计算：`--freq --within groups.txt` 得到群体内频率，再手算 He=2pq；PLINK 1.9 没有直接的 FIS 报告。
- 小样本/近交群体的 F 分布偏斜，用分群体 + 稳健统计（中位数/MAD）解释；不要用全群体均值掩盖品种差异。
- 位点必须剪枝；样本 <50 时用 `--read-freq` 载入可靠频率，否则 F 的期望值有偏。

## 7. 动物育种 / 保种中的解释要点

1. **小群体与品种混合**：所有指标（PCA、PI_HAT、FROH、FST）抽样方差都大，报告要给 n、位点数与不确定性说明；单年数据不足以判断近交趋势，最好有 3–5 个世代/年度的重复。PCA/FST 会强烈反映品种差异，保种分析先按品种分层；跨品种 FST 不能当"保种价值"。
2. **芯片密度**：低密度芯片系统性低估 FROH、LD 与 PI_HAT 的可靠性；做多样性/近交结论前尽量 imputation 到高密度参考面板。
3. **性染色体**：ROH/PCA/亲缘默认只用常染色体；性染色体上另有剂量与单倍型问题（家禽 ZW 尤其），不要混入常染色体指标。
4. **Ne 与保种决策**：Ne 可由 LD 或系谱估计（PLINK 1.9 不直接提供），引用时要写方法与来源；He 高、FROH 低、Ne 大、特有等位基因多是正向指标，但不要用单一指标做淘汰决策，需结合系谱、表型与群体历史。

## 参见

- PLINK 1.9：分层/PCA https://www.cog-genomics.org/plink/1.9/strat ；IBD/ROH https://www.cog-genomics.org/plink/1.9/ibd ；LD https://www.cog-genomics.org/plink/1.9/ld ；基础统计（--het/--ibc/--fst）https://www.cog-genomics.org/plink/1.9/basic_stats
- Weir BS, Cockerham CC. Estimating F-statistics for the analysis of population structure. *Evolution*. 1984;38:1358–1370.
- Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. *Am J Hum Genet*. 2011;88:76–82.
- Purfield DC, Berry DP, McParland S, Bradley DG. Runs of homozygosity and population history in cattle. *BMC Genet*. 2012;13:70.
- Price AL, et al. Principal components analysis corrects for stratification in genome-wide association studies. *Nat Genet*. 2006;38:904–909.
- Patterson N, Price AL, Reich D. Population structure and eigenanalysis. *PLoS Genet*. 2006;2:e190.

## 8. 常见误区 / 检查清单

- [ ] PCA/亲缘/LD/ROH/FST 是否都在剪枝后的常染色体位点上算？剪枝参数是否随物种/芯片密度调整并写入报告？
- [ ] PCA 离群是否与品种/场次分离混为一谈？删除离群是否有 ID 与理由？
- [ ] PI_HAT 阈值是否用本群体分布校准？小样本（位点 <1000）的判定是否标注不可靠？家系群体是否误删真实亲属？
- [ ] LD 衰减是否报告了 r² 类型（r² vs D′）、窗口、N、位点数？ROH 阈值是否随芯片密度/物种调整、FROH 分母是否写明？
- [ ] FST 全局均值是否在剪枝集上算？多群体是否给了两两矩阵？样本量/MAF 下限是否说明？
- [ ] `--het`/`--ibc` 是否在剪枝集上算？小样本是否用 `--read-freq`？家禽 ZW 是否把 Z/W 混入常染色体指标？
- [ ] 所有指标是否与样本量、位点数、软件版本、参数一起报告，保证可复现？
