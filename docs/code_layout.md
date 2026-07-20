# 代码分层

`src/scrna_workflow/` 是底层可复用方法，`steps.py` 把方法组合成分析点，`workflow/` 只记录当前项目实际运行的步骤。

- `runtime.py`：目录、PNG/PDF 保存、颜色与 UMAP 坐标。
- `io.py`：样本表、Cell Ranger 矩阵读取和报告复制。
- `qc.py`：线粒体/核糖体 QC、Scrublet 和过滤。
- `analysis.py`：HVG、PCA、Harmony/BBKNN、聚类、marker。
- `annotation.py`：人工、PanglaoDB 或 CellTypist 注释。
- `visualization.py`：细胞比例、样本级比较和基因表达图。
- `steps.py`：项目级标准分析入口。

换一个项目仍成立的逻辑放 `src/`；项目特异 cluster 判断、一次性基因列表和临时比较放新的 workflow；参数放 `config/project.py`。
