# 示例报告

`final_report.html` 与 `final_report.md` 是用本仓库的演示数据跑出的完整报告
（数量性状流程：intake → qc → popgen → gwas → calibration → report）。

- 生成方式：`python scripts/selftest.py --keep-run`，取其中一个 run 的 `05_reports/`
- 演示数据：`scripts/simulate_demo.py` 按固定种子 `20261023` 生成——400 样本 × 8000 SNP，
  2 个亚群、10 对重复/亲缘、2% 缺失、30 个已知因果 SNP
- 生成环境：Python 3.12.3、PLINK 1.9.0（2026-09-27 构建）、依赖版本见 `requirements.txt`
- `final_report.html` 是自包含文件（内联 CSS + base64 图片）：下载后用浏览器打开即可；
  在 GitHub 网页上直接点开只能看到源码，请用 Raw / Download 后本地打开

示例对应的自检结果：数量性状 20/20 PASS（见 README"准确性与验证"）。
