# config

新项目开始时：

1. 复制 `sample_info.example.tsv` 为 `sample_info.tsv`；
2. 修改 `project.py` 中的 QC、批次校正、聚类和 marker 参数；
3. marker 完成后复制 `celltype_annotation.example.tsv` 为 `celltype_annotation.tsv`，填写生物学注释。

`sample_info.tsv`：

- `raw`：`data/` 下的 Cell Ranger 样本目录名；
- `sample`：写入 AnnData 的样本名；
- `group`：实验分组，如 CON、EXP。

标准输入路径为 `data/<raw>/outs/filtered_feature_bc_matrix/` 或 `filtered_feature_bc_matrix.h5`。
