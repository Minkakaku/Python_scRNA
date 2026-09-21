# 方法网页与论文依据

- Scanpy preprocessing/clustering：https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
- scanpy.pp.neighbors：https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html
- scanpy.tl.leiden：https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html
- scanpy.tl.rank_genes_groups：https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html
- kneed KneeLocator：https://kneed.readthedocs.io/
- Scanpy 论文：PMID 29409532

代码决策：

1. h5ad 保留 raw/counts 与 log-normalized 表达职责，不使用中心化矩阵做 marker。
2. KneeLocator 使用 convex/decreasing PCA 方差曲线；返回 None 时记录未找到，不为得到数字反复调参。
3. 自动 knee 仅作 N_PCS 参考，不覆盖项目配置。
4. CON/EXP 比例和评分以样本级汇总为统计单位。
5. UMAP 坐标固定，最终图保存 PNG/PDF，绘图数据保存 CSV。