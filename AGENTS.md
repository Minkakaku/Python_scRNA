# Python scRNA-seq 项目代码约定

- 跨项目复用函数放在 `src/python/scrna_workflow/`；公共默认放在根目录 `config/`；项目参数和人工决策放在 `projects/<项目>/config/`。
- 项目 `workflow/` 只放编号短入口；已运行脚本不覆盖，改参或重跑时新增下一编号脚本。
- 正式分析只使用 Python 脚本，不使用 notebook 或 Streamlit。
- 原始脚本放在 `legacy/original/`，不得修改；`source_catalog/` 只作为其他人员代码参考。
- 不在复用源码中写绝对路径、项目编号、样本名或固定 cluster 注释。
- CON/EXP 推断与检验使用样本级汇总，不把细胞当作生物学重复。
- CON 默认颜色 `#E64B35`，EXP 默认颜色 `#4DBBD5`。
- 中间 H5AD 写入项目 `data/processed/`；最终图、表和中文说明写入项目编号 `results/` 目录。
- 最终图保存 PNG 和 PDF，绘图输入及结果表保存 CSV。
- 优先保留原始分析意图，代码简洁，避免过度防御性编程。
- 项目目录统一使用小写 `projects/`；默认分类配色使用固定 NPG 色序，连续与双向数值使用 NPG 衍生渐变，项目显式配色可覆盖默认值。
