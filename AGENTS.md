# Python scRNA-seq 项目代码约定

- 所有可复用函数和分析方法放在 `src/`。
- 项目路径、样本分组、QC、聚类参数和人工注释放在 `config/`。
- `workflow/` 初始为空；每运行一个新分析点，新建 `stepXX_name.py`。
- 已运行的 workflow 不覆盖；参数变化或重跑新建下一步脚本。
- workflow 只导入配置并调用一个清楚的方法，避免复制函数。
- 原始脚本放在 `legacy/original/`，不得直接修改。
- 正式分析只使用 Python 脚本，不使用 notebook 或 Streamlit。
- 不在 `src/` 中写绝对路径、项目编号、样本名、固定 cluster 注释。
- CON/EXP 推断与检验使用样本级汇总，不把细胞当作生物学重复。
- CON 默认颜色 `#72b1b9`，EXP 默认颜色 `#f6ae00`。
- 最终图保存 PNG 和 PDF，绘图输入及结果表保存 CSV。
- 优先保留原始分析意图，代码简洁，避免过度防御性编程。
