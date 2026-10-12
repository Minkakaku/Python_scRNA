# 代码分层

`src/python/scrna_workflow/` 是跨项目复用方法，`steps.py` 组合标准分析点，`projects/<项目>/workflow/` 只记录该项目的短执行入口。

- `runtime.py`：目录、PNG/PDF 保存、颜色和 UMAP 坐标。
- `io.py`：样本表、Cell Ranger 矩阵读取和报告复制。
- `qc.py`：线粒体/核糖体 QC、Scrublet 和过滤。
- `analysis.py`：HVG、PCA、Harmony/BBKNN、聚类和 marker。
- `annotation.py`：人工、PanglaoDB 或 CellTypist 注释。
- `visualization.py`：细胞比例、样本级比较和基因表达图。
- `project.py`：公共 YAML 与项目 YAML 的简单覆盖加载。
- `steps.py`：项目级标准分析入口。

换项目仍成立的逻辑放 `src/python/`；项目特异参数、cluster 判断和基因列表放项目 `config/`，一次性改动记录在新的 workflow 脚本中。

- `io.py` 的 `load_expression_data` 扩展四平台读取，旧 10x 接口继续委托同一实现。
- `compatibility.py`：独立 H5AD 兼容副本与可选旧环境验证。
- `vdj.py`：精确细胞匹配、链质控、克隆簇及样本统计，Scirpy/MuData 延迟导入。
