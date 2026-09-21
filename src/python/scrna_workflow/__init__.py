# 集中导出单细胞流程的公开函数；__all__ 指定允许通配导入的接口。

from .analysis import (
    cluster_proportion_table,
    detect_pca_knee,
    save_pca_knee,
    dim_reduction_cluster,
    find_marker_gene,
)
from .annotation import apply_manual_annotation, celltype_annotation
from .io import get_10x_mtx, infer_group, load_10x_data, read_sample_info
from .qc import add_qc_metrics, cell_filtering, data_qc
from .runtime import make_report_dir, save_fig, save_figure
from .visualization import (
    plot_celltype_group_comparison,
    plot_celltype_proportion,
    sample_proportion_table,
    visua_gene_express,
)

__all__ = [
    "make_report_dir",
    "save_fig",
    "save_figure",
    "infer_group",
    "read_sample_info",
    "load_10x_data",
    "get_10x_mtx",
    "add_qc_metrics",
    "data_qc",
    "cell_filtering",
    "dim_reduction_cluster",
    "detect_pca_knee",
    "save_pca_knee",
    "cluster_proportion_table",
    "find_marker_gene",
    "apply_manual_annotation",
    "celltype_annotation",
    "sample_proportion_table",
    "plot_celltype_proportion",
    "plot_celltype_group_comparison",
    "visua_gene_express",
]
