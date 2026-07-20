# 原始 Python 单细胞代码审查

## 原始文件

`legacy/original/scanpy_pipeline(1).py` 未修改，SHA-256：

`78B65016AC43B55E5A4FBA50874AE4C29FF0F4DEC474D69159EB8C1F128A57A9`

## 保留内容

- Cell Ranger 10x 矩阵读取和样本分组；
- 线粒体、核糖体 QC；
- Scrublet 双细胞识别；
- 标准化、HVG、回归、PCA、Harmony/BBKNN、Leiden/Louvain 和 UMAP；
- cluster marker；
- PanglaoDB、CellTypist 和人工注释路线；
- 细胞比例与基因表达可视化；
- 原函数名作为新的公共 API。

## 已修复的明确问题

1. 移除核心代码对 Streamlit、PIL 和页面 fragment 的依赖，正式流程改为纯 `.py`。
2. 删除固定 `/mnt/NexData...` Scanpy 缓存路径，所有路径来自配置。
3. 用 `shutil` 代替未转义的 Linux `cp`、`nohup cp -r`。
4. 用 `anndata.concat` 代替旧的 `AnnData.concatenate`，细胞 ID 改为 `sample:barcode`。
5. 修正原 `save_fig()` 混用 pyplot 与 Figure 的问题，所有最终图保存 PNG 和 PDF。
6. Scrublet 改为每个样本独立运行，再合并预测结果。
7. 原批次校正参数比较的是带 Markdown 星号的 `**harmony**`、`**bbknn**`；现同时接受普通字符串并明确使用 Harmony PCA。
8. `adata.raw` 改为副本，HVG 子集使用 `.copy()`，避免视图和后续修改混淆。
9. split UMAP 使用一致坐标范围和固定纵横比。
10. marker 不再先递归删除整个输出目录；无显著 marker 时不强行绘制空基因列表。
11. CellTypist 使用完整、log-normalized raw 基因矩阵，模型不再按物种写死。
12. 细胞比例组间统计改为样本级，不以细胞作为生物学重复。
13. 删除重复定义的 `visua_gene_express()`，并移除返回 Streamlit 图片对象的逻辑。
14. 原扩展名为 `.xls` 的制表文本改为真实 CSV。
15. 清除固定服务器目录、UI session 状态和项目专属值。

## 仍需按项目决定

- QC 上下限和预期双细胞率；
- 是否需要 Harmony/BBKNN；
- PC、邻居数和 resolution；
- marker 阈值；
- 人工细胞注释；
- 自动注释数据库及模型是否适合组织。
