# 错误手册（05-error-playbook）

> **适用场景**：PLINK 运行失败或结果异常时，按"症状 → 原因 → 自动/手动修复 → 验证"定位问题。本文是 `utils.diagnose_log()` 中文建议的知识库：诊断器按日志子串匹配条目，返回下面"自动修复"与"手动命令"两栏。日志原文以 v1.9.0 实际输出为准，差几个字不影响匹配。

## 0. 修复分级（必须先明确）

| 级别 | 含义 | 允许自动执行 |
|---|---|---|
| A（安全、可逆） | 改文件编码/路径/换行/BOM、复制到 ASCII 短路径、释放占用、调低 `--memory/--threads` | 是 |
| B（语义明确、需记录） | 显式指定 `--chr-set`、`--double-id`、`--missing-phenotype`、`--biallelic-only` | 是，但写入 run.json 并告知用户 |
| C（改变数据） | 加 `--allow-extra-chr`、删样本/位点、放宽阈值、GC 校正 | **否**，必须用户确认并记录在案 |

底线：任何"删数据"的修复都要在报告里留下计数与原因；不允许为了跑通而静默放宽阈值。

## E01 Invalid chromosome code（染色体编码非法）

**症状**：`Error: Invalid chromosome code 'chr1' on line 1 of plink.bim.`（`chrZ`、`chrUn`、`scaffold12` 同理）
**原因**：VCF/注释带 `chr` 前缀；物种染色体数未声明（默认人类 22 对，猪的 23 号会被当 X）；contig/scaffold；全角字符或前后空格。
**自动修复**：扫描 `.bim`/VCF 头部的染色体集合，与 config `species.chr_set` 比对，给出 `--chr-set` 建议；**不自动加 `--allow-extra-chr`**。
**手动**：
```bat
rem 看点名的编码
Get-Content raw.bim | ForEach-Object { ($_ -split '\s+')[0] } | Sort-Object -Unique
rem 按物种声明（推荐转换时一次完成）：猪/牛/羊/鸡
plink --vcf raw.vcf.gz --chr-set 18 --double-id --make-bed --out geno
plink --vcf raw.vcf.gz --cow     --double-id --make-bed --out geno
plink --vcf raw.vcf.gz --sheep   --double-id --make-bed --out geno
plink --vcf raw.vcf.gz --chr-set 38 no-xy --double-id --make-bed --out geno
rem 仍有 contig：允许额外编码（不参与多数分析）；更推荐用 bcftools 先统一命名
plink --vcf raw.vcf.gz --chr-set 18 --allow-extra-chr 0 --double-id --make-bed --out geno
bcftools annotate --rename-chrs rename.txt raw.vcf.gz -Oz -o renamed.vcf.gz
```
`rename.txt` 每行两列：`chr1 1`、`chrX 19`（按物种）、`chrZ 39`。
**验证**：`.bim` 染色体集合 ⊆ 期望集合；X 编码 = 常染色体对数 + 1；样本/位点数与转换前一致；`--check-sex` 能读到 X 位点。

## E02 Duplicate ID（样本 ID 重复）

**症状**：`Error: Duplicate ID 'sampleA'.` / `Error: Duplicate ID 'FAM1 IND1'.`
**原因**：同一 FID+IID 出现两次；VCF 样本名重复；默认 VCF ID 解析是 `--double-id` + `--id-delim`（默认 `_`），含下划线的样本名会被拆成 FID/IID 而撞名；多芯片合并编号体系不同。
**自动修复**：定位重复 ID 并给出重命名表；不自动删样本。
**手动**：
```bat
rem 查重复
Get-Content geno.fam | ForEach-Object { ($_ -split '\s+')[0..1] -join ' ' } | Group-Object | Where-Object Count -gt 1
rem 更新 ID（rename.txt 四列：旧FID 旧IID 新FID 新IID）
plink --bfile geno --update-ids rename.txt --make-bed --out geno_fixed
rem VCF 场景明确解析规则：FID=IID=样本名 / FID=0 / 自定义分隔符
plink --vcf raw.vcf.gz --double-id --make-bed --out geno
plink --vcf raw.vcf.gz --const-fid 0 --make-bed --out geno
plink --vcf raw.vcf.gz --id-delim - --make-bed --out geno
```
**验证**：`.fam` ID 唯一；`plink --bfile geno_fixed --freq` 正常；样本数不变。
**位点 ID 重复**通常不致命，但会让 `--clump`/注释/跨数据集匹配产生歧义：用 `--list-duplicate-vars ids-only` 列出，再 `--set-missing-var-ids @:#` 或 `--update-name` 统一。

## E03 No variants / No people remaining

**症状**：`Error: No variants remaining after main filters.` / `Error: No people remaining after main filters.`
**原因**：阈值叠加过严；`--maf` 与极小样本冲突；case/control 的 `--hwe` 只在对照中检验而对照数为 0；`--chr` 范围与物种编码不符；表型/协变量全缺失导致 `--prune` 后无人；`--extract` 列表与数据不匹配。
**自动修复**：解析日志中每步 removed 数，复现"逐级计数"，指出哪一步删到 0 并给放宽建议。
**手动**：
```bat
rem 逐级跑，看每步剩余数（.log）
plink --bfile raw --geno 0.10 --make-bed --out t1
plink --bfile t1  --hwe 1e-6 --make-bed --out t2
plink --bfile t2  --maf 0.01 --make-bed --out t3
plink --bfile t3  --mind 0.10 --make-bed --out t4
Get-Content t1.log,t2.log,t3.log,t4.log | Select-String 'variants remaining|people remaining'
```
**验证**：每步剩余 >0；最终数 = 输入 − 各条 removed 之和。`--allow-no-vars/--allow-no-samples` 仅用于探查，不得作为交付结果。

## E04 .bed magic number（不是 PLINK1 二进制）

**症状**：`Error: .bed file has wrong magic number.`
**原因**：文件不是 PLINK 1 的 `.bed`（正确前两字节 `0x6C 0x1B`）：可能是 PLINK 2 的 `.pgen`、gzip 压缩文件、VCF 改后缀、UCSC BED（区间文本）、下载不完整或文本模式传输。
**自动修复**：读前 4 字节判断真实格式，给出转换命令；不自动重命名/覆盖原文件。
**手动**：
```powershell
Format-Hex -Path geno.bed -Count 4        rem 期望 6C 1B ...
```
```bat
plink2 --pgen geno.pgen --pvar geno.pvar --psam geno.psam --make-bed --out geno   rem pgen→bed
gzip -dc geno.bed.gz > geno.bed                                                     rem gzip
plink --vcf real.vcf.gz --double-id --make-bed --out geno                           rem 实为 VCF
```
**验证**：前两字节正确；`plink --bfile geno --freq` 成功；计数与来源一致。

## E05 Failed to open（打不开文件）

**症状**：`Error: Failed to open geno.bed.` / `Error: Failed to open --pheno file ...`
**原因**：路径不存在或大小写不符；工作目录不同；文件名含空格/中文；输出目录不存在；OneDrive/杀毒软件占用；权限不足；扩展名被隐藏（实为 `geno.bed.txt`）；路径超过 260 字符。
**自动修复**：`Resolve-Path` 校验 + 创建输出目录 + 复制到 ASCII 短路径；占用类提示关闭软件，不强制。
**手动**：
```powershell
Resolve-Path .\geno.bed
New-Item -ItemType Directory -Force E:\plink_work\run1 | Out-Null
```
```bat
plink --bfile E:\plink_work\geno --freq --out E:\plink_work\run1\check
```
**验证**：`check.log` 生成且含 "Analysis finished"；返回码 0。

## E06 Invalid phenotype（表型编码错误）

**症状**：`Error: Invalid phenotype value ...`，或**不报错但异常**：`NMISS=0`、p 全 NA、日志提示样本无表型。
**原因**：缺失值不是 `-9`（`NA` 会被当缺失，但 `0` 在 case/control 里是缺失、在数量性状里合法）；case/control 用 0/1 却未加 `--1`；表型列有小数/字符串；`--pheno` 第三列不是目标表型（多列需 `--pheno-name`/`--mpheno`）；表头行被当数据。
**自动修复**：识别表型列取值分布 → 判定类型（0/1/2、连续、全缺失）→ 给编码建议与命令；不自动改写原文件。
**手动**：
```bat
rem 默认 1=对照、2=病例、-9=缺失
plink --bfile qc --pheno pheno.txt --pheno-name status --logistic --out test
rem 0/1 编码：--1 使 0=对照、1=病例（并强制 case/control）
plink --bfile qc --pheno pheno.txt --pheno-name status --1 --logistic --out test
rem 其他缺失值编码 / 多列表型
plink --bfile qc --pheno pheno.txt --missing-phenotype -99 --linear --out test
plink --bfile qc --pheno pheno.txt --pheno-name milk_yield --covar covar.txt --linear --out test
```
**验证**：日志给出非缺失表型人数与病例/对照数；`NMISS>0`；数量性状先看分布（直方图/偏度）再建模。

## E07 表型/协变量 ID 对不上

**症状**：无硬报错但**样本全被排除**或协变量没起作用：`NMISS=0`、只输出截距、`Warning: ... people ... missing phenotype`、PC 列全 NA。
**原因**：FID/IID 顺序颠倒；分隔符不同（空格/制表符/逗号）；UTF-8 BOM 粘在第一个 FID；VCF 默认 ID 解析把 `A_01` 拆成 FID=A、IID=01；大小写/前导零差异；文件内重复 ID；样本不在基因型数据中。
**自动修复**：比较两个 ID 集合，输出"仅基因型有/仅表型有/重复"三张清单并给规范化建议；不自动覆盖用户文件。
**手动**：
```powershell
$g = Get-Content geno.fam | ForEach-Object { ($_ -split '\s+')[0..1] -join ' ' }
$p = Get-Content pheno.txt | Select-Object -Skip 1 | ForEach-Object { ($_ -split '\s+')[0..1] -join ' ' }
Compare-Object $g $p | Select-Object -First 20
```
```bat
plink --bfile geno --update-ids idmap.txt --make-bed --out geno_fixed   rem 四列 ID 对照
plink --bfile geno --pheno pheno.txt --pheno-merge --linear --out test  rem 保留原表型值
plink --bfile geno --pheno pheno.txt --allow-no-sex --linear --out test rem 性别缺失时保留表型（谨慎）
```
**验证**：`Compare-Object` 为空；`.assoc` 的 NMISS 接近 QC 后样本数；协变量出现在 `.log`（非 `hide-covar` 时可见）。

## E08 VCF 多等位 / 样本 ID 问题

**关键事实**：PLINK 1.9 读多等位 VCF **默认不报错**——只保留"参考等位 + 最常见替代等位"，其余替代等位被编码为缺失（官方 input 页）。这是最危险的静默丢数据。
**症状**：导入后缺失率异常升高；AF 分布异常；某些位点基因型全缺失；`--vcf-half-call` 报 Half-missing（默认 `error`）；样本名含空格报错。
**自动修复**：统计多等位位点数与半缺失调用数，建议先拆分/过滤再导入；不自动丢弃。
**手动**：
```bat
bcftools norm -m -any -f ref.fa raw.vcf.gz -Oz -o split.vcf.gz
plink --vcf split.vcf.gz --biallelic-only strict --vcf-half-call m ^
  --vcf-idspace-to _ --double-id --make-bed --out geno
plink --vcf split.vcf.gz --vcf-min-qual 20 --vcf-filter --make-bed --out geno   rem 可选
```
**验证**：`.bim` 等位基因只含 A/C/G/T（或声明 indel）；导入前后位点数差 = 多等位/低质量位点数（写入报告）；多等位计数：`bcftools query -f '%ALT\n' raw.vcf.gz | grep -c ','`。

## E09 `--allow-extra-chr` 与 `--chr-set` 的区别（最常混用）

| 维度 | `--chr-set <n> [no-x/no-y/no-xy/no-mt]` | `--allow-extra-chr ['0']` |
|---|---|---|
| 作用 | 重新定义标准染色体集（常染色体数、X/Y/XY/MT 编码） | 允许未识别的额外编码通过 |
| 典型用途 | 猪 18、牛 29、羊 26、鸡 38、马 31 | unplaced contig、scaffold |
| 对 X 的影响 | 正确设置 X 编码与性别/单倍体语义 | 额外编码不参与 X 语义，**不能替代 `--chr-set`** |
| 误用风险 | 数量设错 → X 与常染色体混淆 | 把 X 当 contig → 性别核查失效、X 关联丢失 |
| 正确顺序 | 先 `--chr-set`，仍有 contig 再加 `--allow-extra-chr` | 同左 |

**验证**：`--freq` 正常；X 编码 = 常染色体数 + 1；`--check-sex` 能读到 X；PCA/ROH 用 `--chr 1-N` 排除性染色体与 contig。

## E10 内存 / 线程

**症状**：`Error: Failed to allocate ...`、程序直接退出、日志截断、机器卡死、`--threads` 越多越慢。
**原因**：内存需求随"样本×位点"增长；32 位版本有 2 GB 限制；`--threads` 超过物理核心会争抢；网络盘/机械盘 I/O 瓶颈。
**自动修复**：按物理核心数与可用内存下调 `runtime.memory_mb/threads`；建议先 `--memory` 试跑；不自动换软件。
**手动**：
```bat
rem --memory 是初始工作区大小(MB)，不是硬上限；默认让 PLINK 自行决定
plink --bfile big --memory 8000 --threads 8 --freq --out check
plink --bfile big --chr 1-5 --make-bed --out part1   rem 分块减负
```
**验证**：4/8 线程耗时对比；内存峰值 < 物理内存 80%；必要时改用 PLINK 2。

## E11 中文路径 / 编码 / 换行 / BOM

**症状**：`Failed to open`；ID 比较失败但肉眼看不出差异；日志出现 `锟斤拷`/`��`；第一个样本 ID 前有不可见字符；Excel 另存的 CSV 表型失效。
**原因**：控制台 GBK；文件带 UTF-8 BOM；Excel 存成 GBK 或逗号分隔（PLINK 按空白分词）；中文路径在部分 API/杀毒软件下失败。
**自动修复（A 级）**：复制到 ASCII 短路径；BOM/GBK→UTF-8（无 BOM）转换副本；不改原文件与用户目录。
**手动**：
```powershell
Format-Hex -Path .\pheno.txt -Count 3          rem 期望开头不是 EF BB BF
$txt = Get-Content -Raw -Encoding Default .\pheno_gbk.txt
[System.IO.File]::WriteAllText('E:\plink_work\pheno.txt', $txt, [System.Text.UTF8Encoding]::new($false))
```
```bat
chcp 65001
plink --bfile E:\plink_work\geno --pheno E:\plink_work\pheno.txt --linear --out E:\plink_work\gwas
```
**验证**：无 BOM；ID 与基因型完全一致；输出文件名与日志可读。

## E12 文件权限 / 被占用

**症状**：写 `.log`/输出时 `Failed to open`、`Access is denied`；日志只有一半；同一命令昨天能跑今天不能。
**原因**：输出目录只读或被同步盘锁定；Excel 打开着正在写的 `.csv`；OneDrive 冲突；杀毒软件扫描；另一个 PLINK 进程持有文件。
**自动修复**：检测只读属性与占用进程，建议换到本地目录/新目录；不强制结束进程。
**手动**：
```powershell
Get-ItemProperty E:\plink_work\run1\gwas.log | Select-Object Attributes
Get-Process plink -ErrorAction SilentlyContinue
Copy-Item -Recurse E:\plink_work\run1 E:\plink_work\run1_copy
```
**验证**：新目录重跑成功；日志完整。

## E13 .bim/.fam 行格式错误 / 等位基因非法

**症状**：`Error: Line 10 of plink.bim has fewer tokens than expected.` / `.fam` 同理 / `Error: Invalid allele code ...`
**原因**：`.bim` 必须 6 列、`.fam` 必须 6 列（FID IID 父 母 性别 表型），空白分隔；多字符等位基因（`<DEL>`、`I/D`）或多等位分隔符；Excel 粘贴引入逗号/全角空格；行尾空行。
**自动修复**：统计每行 token 数并定位坏行，给修复建议；不自动重写。
**手动**：
```powershell
Get-Content geno.bim | ForEach-Object { if (($_ -split '\s+').Count -ne 6) { $_ } }
Get-Content geno.fam | ForEach-Object { if (($_ -split '\s+').Count -ne 6) { $_ } }
```
```bat
plink --bfile geno --snps-only just-acgt --make-bed --out geno_snp    rem 去多字符等位基因
plink --bfile geno --set-missing-var-ids @:# --make-bed --out geno_named
```
**验证**：坏行检查为空；`--freq` 正常；报告记录被丢弃的位点数。

## E14 性染色体 / 性别核查异常

**症状**：`--check-sex` 无位点可分析或大量"PROBLEM"；X 关联结果反常；家禽 F 值分布与哺乳动物相反。
**原因**：`--chr-set` 未设对（X 编码错误）；X 的 PAR 未拆（`--split-x` 只内置人类 build，畜禽需自定义坐标或接受偏差）；未做 LD 剪枝（官方警告会给出错误性别判定）；MAF 不可靠；家禽 ZW 与 PLINK 的 XY 模型相反。
**自动修复**：检查 X 编码与位点是否存在；给出剪枝 + check-sex 命令；家禽场景直接标注"不适用"。
**手动**：
```bat
rem X 先剪枝再核查（X 编码：猪 19、牛 30、羊 27、鸡 39）
plink --bfile 01_qc\data\step2_sample --chr 19 --indep-pairwise 50 5 0.2 --out prune_x
plink --bfile 01_qc\data\step2_sample --extract prune_x.prune.in ^
  --check-sex 0.2 0.8 --out 01_qc\tables\sexcheck
plink --bfile geno --split-x b38 no-fail --make-bed --out geno_split   rem 人类 PAR
```
**验证**：F 值双峰（哺乳动物雌低雄高）；不一致样本数写入 metrics；家禽报告写明"ZW 不支持 XY 式性别推断"。

## E15 clump / fst / within 相关错误

| 症状 | 原因 | 修复 |
|---|---|---|
| `--clump` 找不到基因型/无输出 | clump 需要 `--bfile` 算 LD；或报告无 `SNP`/`P` 列 | 加 `--bfile`；`--clump-snp-field`/`--clump-field` 指定列名 |
| `--fst` 报不足 2 个群体 | `--within` 文件格式不对或所有样本同组 | 检查 `FID IID group` 三列、无表头、组数 ≥2 |
| `--within` 报 ID 不存在 | 分组文件 ID 与基因型不一致/BOM/分隔符 | 按 E07 规范化 ID |
| `--extract` 报列表项不存在 | 列表含表头/逗号/重复，或位点已被删 | 一行一个 ID；先与当前数据对照 |
| `--clump` 结果为空 | p1 阈值过严或没有位点达到 p1 | 探索时可放宽 `--clump-p1`，报告中说明 |
| `--adjust` 列全 NA | p 值列名不匹配或报告为空 | 检查输入报告与字段解析 |

## E16 快速索引（日志子串 → 条目）

| 日志关键子串 | 条目 | 日志关键子串 | 条目 |
|---|---|---|---|
| `Invalid chromosome code` | E01 | `Failed to open` | E05/E11/E12 |
| `Duplicate ID` | E02 | `Invalid phenotype` / 无表型提示 | E06/E07 |
| `No variants/people remaining` | E03 | `multiallelic` / `Half-missing call` | E08 |
| `wrong magic number` | E04 | `Failed to allocate` / 卡死 | E10 |
| `fewer tokens than expected` / `Invalid allele code` | E13 | `--check-sex` / X 异常 | E14 |
| `--clump` / `--fst` / `--within` | E15 | `--allow-extra-chr` / `--chr-set` 提示 | E09 |

## 附录：通用诊断顺序（诊断器/人共用）

1. **读日志尾部**：致命错误以 `Error:` 开头或在 `.log` 末尾；先拿原始文本，不要只看返回码。
2. **找计数**：`N people/variants remaining`、`removed due to ...`，确认在哪一步归零。
3. **核对输入**：文件存在性、前 4 字节、行数、ID 集合、染色体编码、表型分布。80% 的"PLINK 报错"其实是输入问题。
4. **最小复现**：只留 `--bfile + --freq` 跑一次，再逐个加 flag 找出触发项。
5. **版本核对**：`plink --version`（本机应为 v1.9.0 64-bit）；PLINK 2 与 1.9 的 flag 语义不同，不要混用。
6. **修复后验证**：返回码 0 + 产物存在 + 计数符合预期 + 报告记录变更；仅"能跑通"不算修复成功。

## 参见

- PLINK 1.9 错误与警告：https://www.cog-genomics.org/plink/1.9/errors
- 输入/染色体集/VCF/表型/协变量：https://www.cog-genomics.org/plink/1.9/input
- 过滤：https://www.cog-genomics.org/plink/1.9/filter ；系统资源：https://www.cog-genomics.org/plink/1.9/other#memory
- Chang CC, et al. Second-generation PLINK. *GigaScience*. 2015;4:7.

## 常见误区 / 检查清单

- [ ] 是否把 `--allow-extra-chr` 当"万能开闸"，导致 X 语义丢失？
- [ ] VCF 多等位是否被默认静默丢基因型？导入前是否统计并拆分？
- [ ] 表型的 0 值究竟是"缺失"还是"合法数值"？case/control 编码与 `--1` 是否一致？
- [ ] 是否因 FID/IID 顺序或 BOM 导致"没报错但样本全被排除"？
- [ ] `.bed` 魔数错误时是否先确认真实格式，而不是反复重装 PLINK？
- [ ] 自动修复是否只做了 A/B 级？有没有在用户不知情时删样本、放宽阈值或加 `--allow-extra-chr`？
- [ ] 每次修复是否写明"症状/原因/命令/验证结果"，并重跑完整阶段核对 metrics？
- [ ] 家禽 ZW 数据是否误用 XY 性别核查与 X 关联模型？
- [ ] 大文件失败时是否先试 `--memory`/`--threads`/分染色体，而不是直接判定数据损坏？
- [ ] 所有修复后是否核对"输入 − 剔除 = 输出"的计数守恒？
