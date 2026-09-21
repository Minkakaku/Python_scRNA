# Python/Scanpy 单细胞主流程分析点

| 顺序 | 模块 | 分析点 | 主要输出 |
|---|---|---|---|
| 01 | io.py | Cell Ranger/h5ad 输入和样本 metadata | 输入审计 |
| 02 | qc.py | QC 指标、阈值和过滤 | 过滤前后图表及 h5ad |
| 03 | analysis.py | 标准化、HVG、PCA、kneed 拐点、邻居图、Leiden、UMAP | PCA 表、聚类图和对象 |
| 04 | analysis.py | rank_genes_groups marker | 完整 marker 表 |
| 05 | annotation.py | 人工 cell-type 映射 | 注释 h5ad 和 UMAP |
| 06 | steps.py | 样本级细胞比例 | 比例表与组间图 |
| 07 | visualization.py | UMAP、基因表达和审计图 | 同名 PNG/PDF 与绘图 CSV |

各项目在自己的 `workflow/` 中保留短步骤入口。`cluster_gpu.py` 不属于本仓库，GPU 空间聚类保留在 Xenium/Visium HD 主流程。
