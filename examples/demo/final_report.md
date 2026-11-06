# verify2_20261107 — PLINK 全流程分析最终报告

> **运行编号**：`verify2_20261107` ｜ **生成时间**：2026-11-07 00:52:36 +0800 ｜ **综合判定**：⚠️ WARN ｜ **物种**：猪
> **运行创建时间**：2026-11-07T00:48:37+08:00 ｜ **报告版本**：1（plink-analysis skill）

## 一、执行摘要

- 数据概况：400 个样本 × 8,000 个位点，格式 bed，物种 猪。
- 质量控制：样本保留率 99.2%（397/400），位点保留率 80.8%（6,464/8,000），总体判定 ✅ PASS。
- 群体结构：PC1 解释 1.60% 方差，PCA 离群样本 2 个，PI_HAT ≥ 0.25 的样本对 8 对。
- GWAS：linear 模型，检验 6,464 个位点，λGC = 1.063，Bonferroni 显著 3 个，最显著位点 snp00324（P = 1.273e-12）。
- 校准与校验：零模型 λGC 均值 0.998，独立校验通过 16/16 项，总体判定 ✅ PASS。

**综合判定：⚠️ WARN**

流程整体可用，但存在警告项；建议按各阶段报告中的处理建议复核后，再用于正式结论。

**建议关注**

- 结果解读时应结合表型测量误差、芯片覆盖度与群体历史，避免将统计关联直接解释为因果。

## 二、阶段状态总览

| 阶段 | 目录 | 状态 | 耗时(s) | 关键结果 | 阶段报告 |
| :--- | :--- | :---: | ---: | :--- | :--- |
| 数据接入与格式识别 | `00_intake/` | ⚠️ WARN | 0.841 | 样本 400，位点 8,000，格式 bed | [00_intake/report.md](../00_intake/report.md) |
| 质量控制 | `01_qc/` | ✅ PASS | 11.83 | 样本 397/400，位点 6,464/8,000 | [01_qc/report.md](../01_qc/report.md) |
| 群体结构与亲缘分析 | `02_popgen/` | ✅ PASS | 7.688 | PC1 方差 1.60%，PI_HAT≥0.25 对 8 | [02_popgen/report.md](../02_popgen/report.md) |
| 全基因组关联分析 | `03_gwas/` | ✅ PASS | 8.317 | λGC 1.063，显著 3，检验位点 6,464，最小 P 1.273e-12 | [03_gwas/report.md](../03_gwas/report.md) |
| 仿真校准与独立校验 | `04_calibration/` | ✅ PASS | 209.16 | 零模型 λ 均值 0.998，校验通过 16/16 | [04_calibration/report.md](../04_calibration/report.md) |

## 三、数据接入与数据概况

> **阶段结论**：数据接入完成：400 个样本 × 8000 个位点。 物种推断：猪（常染色体数 18，性染色体系统 XY，依据：常染色体数 18 与 猪 匹配（无性染色体））。 表型：quantitative，匹配 400/400，编码 continuous numeric。 协变量：2 列（含哑变量）。 分组：2 组。 需注意：ambiguous_strand_snps（详见 00_intake/preflight.json）。

| 项目 | 取值 |
| :--- | :--- |
| 基因型格式 | bed |
| 样本数 | 400 |
| 位点数 | 8,000 |
| 染色体 | 1, 10, 11, 12, 13, 14, 15, 16, 17, 18, …（共 18 项） |
| 推断物种 | 猪 |
| 表型 | 数量性状 |
| 协变量 | sex、age |
| 分组 | A（200）、B（200） |
| 预检状态 | ⚠️ WARN |

**物种推断明细**

| 项目 | 取值 |
| :--- | :--- |
| 推断物种名称 | pig |
| species_inferred.label | 猪 |
| species_inferred.autosome_count | 18 |
| 性别决定系统 | XY |
| 染色体组设置 | — |
| species_inferred.allow_extra_chr | 否 |
| species_inferred.evidence | 常染色体数 18 与 猪 匹配（无性染色体） |
| species_inferred.chromosome_info.n_codes | 18 |
| species_inferred.chromosome_info.numeric_max | 18 |
| species_inferred.chromosome_info.numeric_min | 1 |
| species_inferred.chromosome_info.n_numeric | 18 |
| species_inferred.chromosome_info.has_x | 否 |
| species_inferred.chromosome_info.has_y | 否 |
| species_inferred.chromosome_info.has_z | 否 |
| species_inferred.chromosome_info.has_w | 否 |
| species_inferred.chromosome_info.has_mt | 否 |
| species_inferred.chromosome_info.nonstandard_codes | — |
| species_inferred.chromosome_info.has_nonstandard | 否 |
| species_inferred.chromosome_info.has_chr_prefix | 否 |

**输入文件预检**

| 检查项 | 状态 | 说明 |
| :--- | :--- | :--- |
| genotype_files | ✅ PASS | PLINK 前缀 /home/hatch/plink-runs/verify2_20261107/00_intake/data/base；样本 400，位点 8000 |
| sample_ids_unique | ✅ PASS | 重复 IID 数量 0 |
| variant_ids_unique | ✅ PASS | 重复位点 ID 数量 0（normalize 已自动加后缀） |
| allele_codes | ✅ PASS | 非 A/C/G/T/N 等位基因记录 0 条 |
| ambiguous_strand_snps | ⚠️ WARN | AT/GC 模糊位点 2500 个（链方向无法从数据判断，跨平台合并时需注意） |
| variant_positions | ✅ PASS | BP=0 的位点 0 个 |
| sex_codes | ✅ PASS | 性别编码分布 {1: 211, 2: 189}（0=未知/缺失）；未知性别样本过多时性别检查关闭 |
| chromosome_codes | ✅ PASS | 染色体 ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '2', '3', '4', '5', '6']；推断 猪，chr-set=None |
| phenotype_file | ✅ PASS | 匹配 400/400 个样本（0 行未匹配）；类型=quantitative，编码=continuous numeric |
| phenotype_coding | ✅ PASS | 数量性状 n=400，均值=0.047，标准差=1.462，缺失=0 |
| covariate_file | ✅ PASS | 匹配 400/400；数值列 ['sex', 'age']；分类列自动哑变量编码 [] |
| group_file | ✅ PASS | 匹配 400/400；分组数 2：[{'name': 'A', 'n': 200}, {'name': 'B', 'n': 200}] |

**输入文件清单**

| 文件 | 大小 | SHA-256（前 24 位） |
| :--- | :--- | :--- |
| /home/hatch/plink-demo-data/geno.bed | 781.3 KB | 5a49f43078d3152a75c8609… |
| /home/hatch/plink-demo-data/geno.bim | 206.1 KB | 3eb4b6037046babb0dba729… |
| /home/hatch/plink-demo-data/geno.fam | 7.0 KB | b3c1d99b5edde20565ec24e… |
| /home/hatch/plink-demo-data/pheno.txt | 7.2 KB | cbfde0558e01633af7244eb… |
| /home/hatch/plink-demo-data/covar.txt | 5.9 KB | fb5805b4715dacb4448d150… |
| /home/hatch/plink-demo-data/group.txt | 4.3 KB | 4d6641cf75c20bcffa4de64… |

## 四、质量控制（QC）

> **阶段结论**：质控完成：样本 400→397（保留 99.2%），常染色体位点 8000→6464（保留 80.8%）。 清理后平均 MAF=0.186，平均位点缺失率=0.0202，HWE p<1e-6 比例=0.0000。 样本剔除原因：duplicate 2、heterozygosity_outlier 1。 高亲缘对（PI_HAT≥0.25）：10 对；PCA 离群样本：2 个。 所有质控门限均通过。

### 4.1 质控前后对比

| 项目 | 质控前 | 质控后 | 保留率 |
| :--- | ---: | ---: | :--- |
| 样本数 | 400 | 397 | 99.2% |
| 位点数 | 8,000 | 6,464 | 80.8% |

### 4.2 剔除明细

| 层面 | 剔除原因 | 数量 |
| :--- | :--- | ---: |
| 样本 | 重复样本 | 2 |
| 样本 | heterozygosity_outlier | 1 |
| 位点 | 缺失率 | 0 |
| 位点 | MAF 过低 | 1,144 |
| 位点 | HWE 偏分离 | 392 |
| 位点 | total | 1,536 |

### 4.3 质量指标

| 指标 | 取值 |
| :--- | :--- |
| 样本平均缺失率 | 0.02015 |
| 位点平均缺失率 | 0.02015 |
| 平均 MAF | 0.1863 |
| MAF < 0.01 位点数 | 5 |
| HWE p < 1e-6 位点数 | 0 |
| 杂合度均值 F | 0.00597 |
| 杂合度标准差 | 0.01962 |
| 性别检查 | 未执行 |
| 性别不一致样本数 | 0 |
| 重复样本对数 | 2 |
| 最大 PI_HAT | 1 |
| 平均 PI_HAT | 0.0108 |
| PCA 离群样本数 | 2 |
| LD 剪枝后位点数 | 6,463 |

### 4.4 质控判定

共 7 项判定：✅ PASS 7 项，⚠️ WARN 0 项，❌ FAIL 0 项。

| 判定项 | 状态 | 实测值 | 阈值 | 说明 |
| :--- | :--- | ---: | :--- | :--- |
| sample_retention | ✅ PASS | 0.9925 | {"warn": 0.9, "fail": 0.6} | 保留 397/400 个样本 |
| variant_retention | ✅ PASS | 0.808 | {"warn": 0.8, "fail": 0.3} | 保留 6464/8000 个位点（常染色体） |
| genotype_call_rate | ✅ PASS | 0.9798 | {"warn": 0.97, "fail": 0.95} | 清理后平均位点缺失率 0.0202 |
| hwe_fail_fraction | ✅ PASS | 0 | {"warn": 0.01, "fail": 0.05} | HWE p<1e-6 的位点比例 0.0000 |
| het_outlier_fraction | ✅ PASS | 0.0025 | {"warn": 0.02, "fail": 0.05} | 杂合度离群样本 1 |
| duplicate_samples | ✅ PASS | 2 | {"duplicate_pihat": 0.9} | 已识别并剔除 2 个重复样本（PI_HAT≥0.9；数据质量问题，已处理） |
| pca_outliers | ✅ PASS | 0.00504 | {"warn": 0.02, "fail": 0.05} | PCA 离群样本 2（默认仅标记不删除） |

**图 1**　杂合度（F）分布

![杂合度（F）分布](figures/qc_heterozygosity.png)

**图 2**　HWE p 值 QQ 图（质控诊断）

![HWE p 值 QQ 图（质控诊断）](figures/qc_hwe_qq.png)

**图 3**　MAF 频谱

![MAF 频谱](figures/qc_maf_spectrum.png)

**图 4**　样本/位点缺失率分布

![样本/位点缺失率分布](figures/qc_missingness.png)

**图 5**　pca

![pca](figures/qc_pca.png)

**图 6**　PCA 碎石图

![PCA 碎石图](figures/qc_pca_scree.png)

**图 7**　亲缘关系（PI_HAT）分布

![亲缘关系（PI_HAT）分布](figures/qc_relatedness.png)

## 五、群体结构与亲缘分析

> **阶段结论**：群体结构：PCA 前 5 个主成分解释 3.5% 方差；离群样本 2 个。 亲缘：PI_HAT 均值 0.01065，最大值 0.6635，≥0.25 的对数 8。 LD 衰减：50kb 处平均 r²=0.00281（550836 对）。 ROH：平均每个个体 0.0 段，平均 FROH=0.0（最大值 0.0）。 FST：全局均值 0.01226（加权 0.01517）。

### 5.1 主成分分析（PCA）

| 项目 | 取值 |
| :--- | :--- |
| 主成分数 | 20 |
| 离群样本数 | 2 |
| PC1 方差解释 | 1.603% |
| PC2 方差解释 | 0.489% |
| PC3 方差解释 | 0.468% |
| PC4 方差解释 | 0.466% |
| PC5 方差解释 | 0.455% |
| PC6 方差解释 | 0.452% |
| PC7 方差解释 | 0.441% |
| PC8 方差解释 | 0.430% |
| PC9 方差解释 | 0.424% |
| PC10 方差解释 | 0.385% |
| PC1–PC5 累计方差（占总方差） | 3.48% |

**离群样本**：ind086, ind171

### 5.2 亲缘关系

| 项目 | 取值 |
| :--- | ---: |
| PI_HAT ≥ 0.25 的样本对数 | 8 |
| 最大 PI_HAT | — |

### 5.3 连锁不平衡（LD）

共统计 550,836 个位点对；按物理距离分箱的平均 r² 如下表。

### 5.4 纯合片段（ROH）与近交

| 指标 | 取值 |
| :--- | ---: |
| ROH 数量 | — |
| 平均每样本 ROH 数 | 0 |
| 平均 FROH | 0 |

| 样本 IID | FROH | ROH 数量 | ROH 总长度（kb） |
| :--- | ---: | ---: | ---: |
| ind001 | 0 | 0 | 0 |
| ind002 | 0 | 0 | 0 |
| ind003 | 0 | 0 | 0 |
| ind004 | 0 | 0 | 0 |
| ind005 | 0 | 0 | 0 |
| ind006 | 0 | 0 | 0 |
| ind007 | 0 | 0 | 0 |
| ind008 | 0 | 0 | 0 |
| ind009 | 0 | 0 | 0 |
| ind010 | 0 | 0 | 0 |

### 5.5 群体分化（FST）

| 指标 | 取值 |
| :--- | ---: |
| 全局 FST | — |
| 参与位点数 | 6,464 |

| 位点 | 染色体 | 物理位置 | FST |
| :--- | ---: | ---: | ---: |
| snp05692 | 13 | 53,455,000 | 0.4857 |
| snp02332 | 6 | 22,655,000 | 0.4843 |
| snp01599 | 4 | 14,990,000 | 0.4794 |
| snp01761 | 4 | 15,800,000 | 0.4471 |
| snp03230 | 8 | 31,145,000 | 0.4387 |
| snp02103 | 5 | 19,510,000 | 0.4354 |
| snp03983 | 9 | 36,910,000 | 0.4182 |
| snp04454 | 11 | 43,265,000 | 0.4174 |
| snp02277 | 6 | 22,380,000 | 0.4055 |
| snp05530 | 13 | 52,645,000 | 0.4036 |

### 5.6 遗传多样性

| 指标 | 取值 |
| :--- | ---: |
| 平均杂合度 | 0.2588 |
| 分组数 | — |

**图 8**　FST 分布

![FST 分布](figures/popgen_fst.png)

**图 9**　LD 衰减曲线

![LD 衰减曲线](figures/popgen_ld_decay.png)

**图 10**　pca

![pca](figures/popgen_pca.png)

**图 11**　PCA 碎石图

![PCA 碎石图](figures/popgen_pca_scree.png)

**图 12**　亲缘关系（PI_HAT）分布

![亲缘关系（PI_HAT）分布](figures/popgen_relatedness.png)

**图 13**　样本 FROH 分布

![样本 FROH 分布](figures/popgen_roh.png)

## 六、全基因组关联分析（GWAS）

> **阶段结论**：GWAS（linear）完成：397 个样本、6464 个位点，协变量 ['sex', 'age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']。 多重检验：Bonferroni 阈值 7.74e-06，显著位点 3 个，FDR<0.05 位点 4 个，建议阈值 1.55e-04 下 4 个。 基因组膨胀因子 λGC=1.0627。 最强信号：snp00324（1:2615000，P=1.273e-12，效应=-0.7071）。 LD clump：3 个独立信号区间。

### 6.1 模型与参数

| 项目 | 取值 |
| :--- | :--- |
| 关联模型 | linear |
| 表型类型 | quantitative |
| 参与样本数 | 397 |
| 检验位点数 | 6,464 |
| 协变量 | sex, age, PC1, PC2, PC3, PC4, PC5 |
| λGC | 1.063 |
| λGC 来源 | p-values (median chi-square) |
| 总体判定 | ✅ PASS |

### 6.2 显著性阈值与结果计数

| 阈值类型 | 数值 |
| :--- | ---: |
| bonferroni | 7.735e-06 |
| suggestive | 1.547e-04 |
| fdr_0.05 | 0.05 |

| 统计口径 | 数量 |
| :--- | ---: |
| Bonferroni 显著 | 3 |
| FDR < 0.05 | 4 |
| 建议性阈值 | 4 |

**P 值分布**

| 项目 | 取值 |
| :--- | ---: |
| 最小 P 值 | 1.273e-12 |
| P 值中位数 | 0.4869 |
| P < 5e-8 的位点数 | 2 |

### 6.3 Lead SNP

| 位点 | 染色体 | 物理位置 | A1 | 等位基因频率 | 效应量 beta | 标准误 | P 值 | n |
| :--- | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: |
| snp00324 | 1 | 2,615,000 | T | 0.3615 | -0.7071 | 0.09629 | 1.273e-12 | 390 |
| snp00095 | 1 | 1,470,000 | G | 0.7238 | 0.7255 | 0.1113 | 2.293e-10 | 391 |
| snp00331 | 1 | 2,650,000 | G | 0.394 | 0.5113 | 0.1048 | 1.590e-06 | 382 |
| snp04973 | 12 | 47,860,000 | A | 0.4158 | 0.4547 | 0.1052 | 1.989e-05 | 392 |
| snp02761 | 7 | 26,800,000 | A | 0.6851 | -0.4189 | 0.1113 | 1.949e-04 | 389 |
| snp07359 | 17 | 69,790,000 | G | 0.8453 | 0.5194 | 0.1438 | 3.455e-04 | 391 |
| snp06748 | 16 | 64,735,000 | T | 0.05711 | 0.7861 | 0.2179 | 3.499e-04 | 394 |
| snp01240 | 3 | 11,195,000 | T | 0.8338 | 0.4906 | 0.1364 | 3.651e-04 | 394 |
| snp00994 | 3 | 9,965,000 | C | 0.2258 | -0.4432 | 0.1282 | 6.097e-04 | 383 |
| snp02998 | 7 | 27,985,000 | C | 0.3432 | 0.3905 | 0.1143 | 7.013e-04 | 389 |

### 6.4 结果解读

λGC = 1.063，p 值分布与零假设基本一致，未见明显系统性膨胀。

**图 14**　全基因组曼哈顿图（-log10 P，含显著性阈值）

![全基因组曼哈顿图（-log10 P，含显著性阈值）](figures/gwas_manhattan.png)

**图 15**　GWAS p 值 QQ 图（含 λGC 标注）

![GWAS p 值 QQ 图（含 λGC 标注）](figures/gwas_qq.png)

## 七、仿真校准与独立校验

> **阶段结论**：零模型校准：10 次重复，λ 均值 0.998（SD 0.035）；α=0.05 经验假阳性率 0.0488。 已知 QTL 检出力：h²=0.5、20 个因果位点，5e-8 阈值平均 power=0.145，1e-4 阈值平均 power=0.150。 因果位点排名一致性：AUC=0.773，top100 命中率=0.255。 独立重算校验：16 项，通过 16、警告 0、失败 0。

### 7.1 零模型校准

| 指标 | 取值 |
| :--- | ---: |
| 重复次数 | 10 |
| λGC 均值 | 0.998 |
| λGC 标准差 | 0.035 |
| α = 0.05 经验假阳性率 | 0.04878 |
| α = 1e-3 经验假阳性率 | 0.001037 |
| α = 5e-8 经验假阳性率 | 0 |

零模型 λGC 均值 0.998 落在可接受区间 [0.97, 1.03]，提示检验统计量校准良好。

### 7.2 已知 QTL 功效

| 指标 | 取值 |
| :--- | ---: |
| 设定遗传力 h² | 0.5 |
| 因果位点数 | 20 |
| α = 5e-8 功效均值 | 0.145 |
| α = 5e-8 功效标准差 | 0.01581 |
| α = 1e-4 功效均值 | 0.15 |
| α = 1e-4 功效标准差 | 2.926e-17 |

### 7.3 效应量与排名一致性

| 指标 | 取值 |
| :--- | ---: |
| 平均 beta 比值 | 1.035 |
| SE 校准系数 | 1.35 |
| 因果位点进入前 100 的比例 | 0.255 |
| 因果位点排名 AUC | 0.7728 |

### 7.4 独立重算校验

| 指标 | 取值 |
| :--- | :--- |
| 校验项数 | 16 |
| 通过 | 16 |
| 警告 | 0 |
| 失败 | 0 |
| 总体判定 | ✅ PASS |

| 校验项 | 状态 | 实测值 | 阈值 | 说明 |
| :--- | :--- | :--- | :--- | :--- |
| structure_counts | ✅ PASS | {"n_samples": 397, "n_variants": 6464, "n_chromosomes": 18, "bed_bytes": 646403} | .bed 字节数 = 3 + n_v*ceil(n_s/4)，且与 metrics/PLINK 日志一致 | 397 样本 x 6464 位点（18 条染色体），.bed 646403 字节，魔数 6c1b01 |
| duplicate_ids | ✅ PASS | {"duplicate_samples": 0, "duplicate_snps": 0} | 0 | 重复样本 ID 0 个，重复 SNP ID 0 个 |
| allele_legality | ✅ PASS | {"invalid": 0, "identical": 0} | A1/A2 ∈ {A,C,G,T} 且 A1≠A2 | 6464 个位点等位基因合法 |
| qc_filter_arithmetic | ✅ PASS | {"n_samples_out": 397, "n_variants_out": 6464, "samples_removed": 3, "variants_… | n_in - removed == 文件实际计数 | 样本/位点移除计数与 QC 后文件一致 |
| ld_prune_counts | ✅ PASS | {"n_pruned": 6463, "n_removed": 1} | prune.in 计数 == metrics.pruned.n_variants | LD 剪枝保留 6463 个位点，剔除 1 个 |
| maf_recompute | ✅ PASS | {"max_abs_maf_diff": 0.0, "allele_count_mismatch": 0, "n_snps": 100, "n_allele_… | 等位基因计数精确一致；max\\|ΔMAF\\| <= 1e-05 | 随机 100 个位点：A1 等位基因计数/缺失数逐位点与 --freqx 精确比对，不一致 0 项；重算 MAF 最大绝对差 0.00e+00（容差 1e-05） |
| missingness_recompute | ✅ PASS | {"max_fraction_diff": 0.0, "variant_mismatch": 0, "sample_mismatch": 0, "n_snps… | 整数 N_MISS 精确一致；max\\|ΔF_MISS\\| <= 1e-06 | 2000 个位点、397 个样本的缺失数分别与 --missing 的 N_MISS 精确比对（不一致 0 项）；子集总体缺失率 0.0204 |
| hwe_extremes | ✅ PASS | 0 | {"warn": 0.005, "fail": 0.02} | QC 后 6464 个位点（全部样本(ALL)）HWE(midp) p<1e-6 比例 0.0000%（warn>0.5%，fail>2.0%） |
| phenotype_coding | ✅ PASS | {"type": "quantitative", "n_rows": 400, "n_in_fam": 397, "n_usable": 397, "n_un… | 二元 1/2 编码且患病率合理；数量性状取值多样、无缺失异常 | 表型类型=quantitative，400 行，匹配 QC 后样本 397，可用 397 |
| pca_recompute | ✅ PASS | {"n_samples": 397, "n_pcs": 20, "k_checked": 5, "eigenvector_source": "01_qc/ta… | {"var_explained": "eigval / N_samples", "pca_corr": 0.95, "eigen_scale": 0.15} | 397 样本 x 20 个 PC；PC1 特征值 6.3628，方差解释 1.603%；独立重算方向相关 PC1=0.999998（6463 位点） |
| relatedness_range | ✅ PASS | {"n_pairs": 78606, "min": 0.0, "max": 0.6635, "mean": 0.010647744446988778, "n_… | PI_HAT ∈ [-0.05, 1.05]，计数与 metrics 一致 | 78606 个样本对：PI_HAT ∈ [0.000, 0.663]，>0.25 有 8 对，>0.90 有 0 对（近重复） |
| kinship_matrix | ✅ PASS | {"n": 397, "asymmetry": 0.0, "diag_mean": 1.0048829949622164, "diag_max_dev": 0… | 对称（<1e-6）、对角≈1、非对角 ∈ [-2, 2] | 397x397 关系矩阵：对称性 0.0e+00，对角均值 1.0049，非对角范围 [-0.069, 0.720] |
| gwas_assoc_table | ✅ PASS | {"n_rows": 6464, "n_significant": 3, "alpha_used": 7.735148514851485e-06} | 行数 == metrics.n_tested，NMISS/显著性计数一致 | 6464 个位点的加性检验结果 |
| gwas_beta_se_recompute | ✅ PASS | {"model": "linear", "n_compared": 20, "max_abs_beta_diff": 4.770030011341797e-0… | max\\|Δ\\| <= 0.0001 | linear 模型，20 个随机 SNP：max\\|Δbeta\\|=4.77e-05，max\\|ΔSE\\|=4.83e-05（容差 1e-04） |
| lambda_gc_recompute | ✅ PASS | {"lambda_independent": 1.0627431181240865, "comparisons": {"metrics.lambda_gc":… | 相对容差 0.005 | 独立重算 λGC=1.06274，metrics.lambda_gc=1.06270，metrics.lambda_plink=1.06233，plink.log.lambda=1.06233，adjusted.recompute=1.0… |
| gwas_results_table | ✅ PASS | {"n_rows": 6464, "n_lead_snps": 10} | 列齐全、行数一致、P 与 .assoc 一致、lead SNP 可追溯 | gwas_results.csv 6464 行，含 10 个 lead SNP |

**图 16**　效应量校准（β_est / β_true）

![效应量校准（β_est / β_true）](figures/calibration_effect_bias.png)

**图 17**　零模型 λGC 分布（重复模拟）

![零模型 λGC 分布（重复模拟）](figures/calibration_null_lambda.png)

**图 18**　已知 QTL 检出力（power）曲线

![已知 QTL 检出力（power）曲线](figures/calibration_power.png)

## 八、方法学与参数

### 8.1 方法学依据

- PLINK 1.9 是全基因组关联与群体遗传分析的标准工具链（Chang et al., 2015, *GigaScience* 4:7；https://www.cog-genomics.org/plink/1.9/）。
- 质控阈值的推荐范围参考 Anderson et al., 2010, *PLoS One* 5:e11342（样本/位点过滤与 HWE 处理）。
- PCA 作为群体分层校正的标准做法（Price et al., 2006, *Nature Genetics* 38:904；Patterson et al., 2006, *PLoS Genetics* 2:e190）。
- 亲缘/IBD 估计与重复样本识别使用 PI_HAT（Purcell et al., 2007, *AJHG* 81:559）。
- 多重检验采用 Bonferroni 与 Benjamini–Hochberg FDR；建议性阈值仅供探索。
- ROH/FROH 作为近交代理指标，需注意芯片密度与群体历史对估计的影响（Browning & Browning, 2011）。
- PLINK 1.9 不含线性混合模型（LMM）；强亲缘/强分层数据建议使用 GCTA、GEMMA、BOLT-LMM 等复核（Yang et al., 2011, *AJHG* 88:76）。
- 仿真校准与独立重算用于验证流程实现，不替代独立队列的重复验证。

### 8.2 关键参数

**质控参数**

| 参数 | 取值 | 说明 |
| :--- | :--- | :--- |
| qc.mind | 0.1 | 样本缺失率上限（PLINK --mind） |
| qc.geno | 0.1 | 位点缺失率上限（PLINK --geno） |
| qc.maf | 0.01 | 最小等位基因频率下限（PLINK --maf） |
| qc.hwe | 1.000e-06 | HWE 检验 p 值下限（PLINK --hwe） |
| qc.het_sd | 3 | 杂合度离群判定标准差倍数（--het） |
| qc.duplicate_pihat | 0.9 | 重复样本 PI_HAT 阈值 |
| qc.related_pihat | 0.25 | 近亲样本 PI_HAT 阈值 |
| qc.sex_check | auto | 性别一致性检查模式（--check-sex） |
| qc.pca_outlier_sd | 6 | PCA 离群判定标准差倍数 |
| qc.ld_prune.window | 50 | LD 剪枝窗口（SNP 数，--indep-pairwise） |
| qc.ld_prune.step | 5 | LD 剪枝步长（SNP 数） |
| qc.ld_prune.r2 | 0.2 | LD 剪枝 r² 阈值 |
| qc.remove_related | 否 | 是否剔除近亲样本 |
| qc.remove_pca_outliers | 否 | 是否剔除 PCA 离群样本 |

**群体结构参数**

| 参数 | 取值 | 说明 |
| :--- | ---: | :--- |
| popgen.pca | 20 | 输出主成分数（--pca） |
| popgen.ld_window_kb | 1,000 | LD 计算窗口（kb，--r2） |
| popgen.ld_bin_kb | 50 | LD 距离分箱宽度（kb） |
| popgen.roh.snp | 50 | ROH 检测最少 SNP 数（--homozyg） |
| popgen.roh.kb | 1,000 | ROH 检测最小长度（kb） |
| popgen.roh.density | 50 | ROH 区域 SNP 密度阈值（每 kb） |
| popgen.roh.gap | 1,000 | ROH 允许的最大间隔（kb） |
| popgen.roh.het_missing | 0 | ROH 内允许的杂合/缺失比例 |

**GWAS参数**

| 参数 | 取值 | 说明 |
| :--- | :--- | :--- |
| gwas.model | auto | 关联模型：auto/linear/logistic |
| gwas.covariates | auto | 协变量设置（auto 表示自动选择） |
| gwas.n_pcs | 5 | 纳入 GWAS 的 PCA 主成分数 |
| gwas.ci | 0.95 | 效应量置信区间水平 |
| gwas.maf | 0.01 | 最小等位基因频率下限（PLINK --maf） |
| gwas.geno | 0.1 | 位点缺失率上限（PLINK --geno） |
| gwas.hwe | 1.000e-06 | HWE 检验 p 值下限（PLINK --hwe） |
| gwas.clump.p1 | 5.000e-08 | clump 显著阈值 p1 |
| gwas.clump.p2 | 1.000e-04 | clump 次级阈值 p2 |
| gwas.clump.r2 | 0.2 | clump 的 r² 聚类阈值 |
| gwas.clump.kb | 1,000 | clump 的物理距离窗口（kb） |
| gwas.genomic_control | none | 基因组控制（none / beta / 1） |

**校准参数**

| 参数 | 取值 | 说明 |
| :--- | :--- | :--- |
| calibration.enabled | 是 | 是否启用仿真校准 |
| calibration.n_replicates | 10 | 仿真重复次数 |
| calibration.n_causal | 20 | 因果位点数 |
| calibration.h2 | 0.5 | 设定遗传力 h² |
| calibration.n_snps_subset | 20,000 | 仿真使用的 SNP 子集大小 |
| calibration.alpha | 5.000e-08 | 显著性阈值 α |
| calibration.seed | 20,261,023 | — |
| calibration.truth_path | — | — |
| calibration.model | auto | 关联模型：auto/linear/logistic |

**运行时参数**

| 参数 | 取值 | 说明 |
| :--- | :--- | :--- |
| runtime.threads | 4 | PLINK 线程数 |
| runtime.memory_mb | 4,000 | PLINK 内存上限（MB） |
| runtime.timeout_s | 7,200 | 单步命令超时（秒） |
| runtime.keep_intermediate | 是 | 是否保留中间文件 |

## 九、复现命令

在 skill 根目录（`plink-analysis/`）下执行：

```powershell
cd /home/hatch/workspace/plink-analysis
python scripts\plink_pipeline.py auto --input "/home/hatch/plink-demo-data/geno" --pheno "/home/hatch/plink-demo-data/pheno.txt" --covar "/home/hatch/plink-demo-data/covar.txt" --group "/home/hatch/plink-demo-data/group.txt" --species pig --run-id verify2_20261107 --out-root "/home/hatch/plink-runs"
```

对已有 run 目录重跑或补跑阶段：

```powershell
python scripts\plink_pipeline.py run --run-dir "/home/hatch/plink-runs/verify2_20261107" --stage all --force
python scripts\plink_pipeline.py report --run-dir "/home/hatch/plink-runs/verify2_20261107"
python scripts\plink_pipeline.py status --run-dir "/home/hatch/plink-runs/verify2_20261107"
```

> 复现需保持相同版本的 PLINK（1.9）、Python 依赖与 run.json 中的参数；随机种子与仿真设置记录在 `run.json` / `04_calibration/metrics.json` 中。

## 十、局限性

- 统计关联不等于因果：GWAS 显著位点仅提示候选区域，仍需独立群体重复、功能注释与实验验证。
- PLINK 1.9 不提供混合模型（LMM）；PCA 协变量只能部分校正群体分层与亲缘结构，强结构群体中可能仍有残余膨胀。
- 芯片标记的 ascertainment bias、未填补基因型与低 MAF 位点会限制稀有变异与因果位点的检出。
- 多重检验 Bonferroni 较为保守；FDR 与建议性阈值提高了发现数但增加了假阳性风险，需在结论中区分。
- 样本量决定检验功效；小样本下的效应量估计与显著性均不稳定。
- 质控阈值、LD 剪枝参数与主成分数量会影响结果；本报告固定记录参数，但未做全参数敏感性分析。
- 表型测量误差、缺失机制与协变量编码依赖上游数据质量，本流程无法纠正系统性表型偏差。
- 仿真校准验证的是流程实现与统计性质，不能替代独立队列或前瞻性验证；仿真中的遗传结构是对真实群体的简化。
- FROH、PI_HAT 等指标受芯片密度与群体历史影响，跨研究比较时应谨慎。

## 十一、文件与产物清单

| 文件 | 大小 | 说明 |
| :--- | :--- | :--- |
| run.json | 3.4 KB | 运行配置 / 状态 |
| state.json | 5.2 KB | 运行配置 / 状态 |
| 00_intake/base_prefix.txt | 59 B | 文件 |
| 00_intake/data/base.bed | 781.3 KB | 中间数据 |
| 00_intake/data/base.bim | 206.1 KB | 中间数据 |
| 00_intake/data/base.fam | 7.0 KB | 中间数据 |
| 00_intake/data/base.log | 1.0 KB | 中间数据 |
| 00_intake/data/covar_clean.txt | 5.9 KB | 中间数据 |
| 00_intake/data/group_clean.txt | 4.3 KB | 中间数据 |
| 00_intake/data/pheno_clean.txt | 7.2 KB | 中间数据 |
| 00_intake/data/within.txt | 4.3 KB | 中间数据 |
| 00_intake/manifest.json | 2.1 KB | 文件 |
| 00_intake/metrics.json | 5.7 KB | 阶段指标（机器可读） |
| 00_intake/preflight.json | 2.0 KB | 文件 |
| 00_intake/report.md | 10.7 KB | 阶段报告（Markdown） |
| 01_qc/data/clean.bed | 631.3 KB | 中间数据 |
| 01_qc/data/clean.bim | 166.5 KB | 中间数据 |
| 01_qc/data/clean.fam | 7.0 KB | 中间数据 |
| 01_qc/data/clean.log | 1.2 KB | 中间数据 |
| 01_qc/data/clean_pruned.log | 2.0 KB | 中间数据 |
| 01_qc/data/clean_pruned.prune.in | 56.8 KB | 中间数据 |
| 01_qc/data/clean_pruned.prune.out | 9 B | 中间数据 |
| 01_qc/data/exclude_variants.txt | 13.5 KB | 中间数据 |
| 01_qc/data/pca_covar.txt | 24.2 KB | 中间数据 |
| 01_qc/data/pruned.log | 2.0 KB | 中间数据 |
| 01_qc/data/pruned.prune.in | 56.8 KB | 中间数据 |
| 01_qc/data/pruned.prune.out | 9 B | 中间数据 |
| 01_qc/data/qc_variants.bed | 631.3 KB | 中间数据 |
| 01_qc/data/qc_variants.bim | 166.5 KB | 中间数据 |
| 01_qc/data/qc_variants.fam | 7.0 KB | 中间数据 |
| 01_qc/data/qc_variants.log | 1.2 KB | 中间数据 |
| 01_qc/data/remove_samples.txt | 70 B | 中间数据 |
| 01_qc/metrics.json | 4.2 KB | 阶段指标（机器可读） |
| 01_qc/plots/heterozygosity.png | 48.1 KB | 图（PNG） |
| 01_qc/plots/hwe_qq.png | 67.0 KB | 图（PNG） |
| 01_qc/plots/maf_spectrum.png | 48.2 KB | 图（PNG） |
| 01_qc/plots/missingness.png | 64.1 KB | 图（PNG） |
| 01_qc/plots/pca.png | 74.2 KB | 图（PNG） |
| 01_qc/plots/pca_scree.png | 65.6 KB | 图（PNG） |
| 01_qc/plots/relatedness.png | 54.6 KB | 图（PNG） |
| 01_qc/report.md | 14.5 KB | 阶段报告（Markdown） |
| 01_qc/tables/clean_pca.eigenval | 160 B | 数据表（CSV） |
| 01_qc/tables/clean_pca.eigenvec | 85.6 KB | 数据表（CSV） |
| 01_qc/tables/clean_pca.log | 1.4 KB | 数据表（CSV） |
| 01_qc/tables/final_stats.frq | 303.0 KB | 数据表（CSV） |
| 01_qc/tables/final_stats.het | 25.7 KB | 数据表（CSV） |
| 01_qc/tables/final_stats.hwe | 549.3 KB | 数据表（CSV） |
| 01_qc/tables/final_stats.imiss | 20.2 KB | 数据表（CSV） |
| 01_qc/tables/final_stats.lmiss | 271.5 KB | 数据表（CSV） |
| 01_qc/tables/final_stats.log | 1.5 KB | 数据表（CSV） |
| 01_qc/tables/genome.genome | 7.3 MB | 数据表（CSV） |
| 01_qc/tables/genome.log | 1.1 KB | 数据表（CSV） |
| 01_qc/tables/het_clean.csv | 12.7 KB | 数据表（CSV） |
| 01_qc/tables/pca.eigenval | 156 B | 数据表（CSV） |
| 01_qc/tables/pca.eigenval.csv | 716 B | 数据表（CSV） |
| 01_qc/tables/pca.eigenvec | 85.7 KB | 数据表（CSV） |
| 01_qc/tables/pca.log | 1.4 KB | 数据表（CSV） |
| 01_qc/tables/raw_stats.frq | 375.0 KB | 数据表（CSV） |
| 01_qc/tables/raw_stats.het | 25.8 KB | 数据表（CSV） |
| 01_qc/tables/raw_stats.hwe | 679.8 KB | 数据表（CSV） |
| 01_qc/tables/raw_stats.imiss | 20.4 KB | 数据表（CSV） |
| 01_qc/tables/raw_stats.lmiss | 336.0 KB | 数据表（CSV） |
| 01_qc/tables/raw_stats.log | 1.5 KB | 数据表（CSV） |
| 01_qc/tables/related_pairs.csv | 273 B | 数据表（CSV） |
| 01_qc/tables/removed_samples.csv | 89 B | 数据表（CSV） |
| 01_qc/tables/removed_variants.csv | 13.5 KB | 数据表（CSV） |
| 01_qc/tables/sample_qc_clean.csv | 10.9 KB | 数据表（CSV） |
| 01_qc/tables/sample_stats.frq | 303.0 KB | 数据表（CSV） |
| 01_qc/tables/sample_stats.het | 25.8 KB | 数据表（CSV） |
| 01_qc/tables/sample_stats.imiss | 20.4 KB | 数据表（CSV） |
| 01_qc/tables/sample_stats.lmiss | 271.5 KB | 数据表（CSV） |
| 01_qc/tables/sample_stats.log | 1.3 KB | 数据表（CSV） |
| 01_qc/tables/variant_qc_clean.csv | 238.9 KB | 数据表（CSV） |
| 01_qc/tables/variant_qc_raw.csv | 283.0 KB | 数据表（CSV） |
| 02_popgen/data/ld_subset.txt | 56.8 KB | 中间数据 |
| 02_popgen/metrics.json | 10.9 KB | 阶段指标（机器可读） |
| 02_popgen/plots/fst.png | 59.3 KB | 图（PNG） |
| 02_popgen/plots/ld_decay.png | 25.4 KB | 图（PNG） |
| 02_popgen/plots/pca.png | 69.7 KB | 图（PNG） |
| 02_popgen/plots/pca_scree.png | 65.6 KB | 图（PNG） |
> 另有 166 个文件未列出。

**阶段报告索引**

- [数据接入与格式识别报告](../00_intake/report.md)（⚠️ WARN）
- [质量控制报告](../01_qc/report.md)（✅ PASS）
- [群体结构与亲缘分析报告](../02_popgen/report.md)（✅ PASS）
- [全基因组关联分析报告](../03_gwas/report.md)（✅ PASS）
- [仿真校准与独立校验报告](../04_calibration/report.md)（✅ PASS）

---

> 本报告由 [plink-analysis skill](https://www.cog-genomics.org/plink/1.9/) 自动生成，所有数值以各阶段 `metrics.json` 与 `tables/` 中的原始结果为准。阶段报告的详细诊断信息请查阅对应 `report.md`。
