# 公共配置

- `default.yaml`：QC、Scrublet、聚类、marker 和注释方法的公共默认参数。
- `reference.yaml`：物种和参考相关默认值。

项目参数放在 `projects/<项目>/config/config.yaml`，同名顶级键覆盖公共值。样本表和人工 cell-type 映射属于项目配置，不放在根目录。

- `input_platform`：`10x`、`c4`、`seek`、`bd`；样本表的可选 `platform` 列可逐行覆盖。
- `matrix_join`：默认 `inner`，项目明确需要基因并集时改为 `outer`。
- `celltype_marker_logfc` / `celltype_marker_min_pct`：细胞类型正向 marker 阈值，FDR 复用 `marker_fdr`。
- `vdj_method`：`TCR` 或 `BCR`；`vdj_metric`：`alignment` 或 `identity`；`vdj_cutoff` 仅作用于 alignment；`vdj_min_network_cells` 控制网络展示阈值。
- 项目可设置 `legacy_python` 为旧环境 Python 的完整路径，供兼容导出回读验证。
