from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPORT_STEPS = (
    "01-CellRanger",
    "02-Data_QC",
    "03-Cell_filter",
    "04-Dimension_reduction_cluster",
    "05-Marker_gene_analysis",
    "06-Cell_type_annotation",
    "07-Cell_proportion",
)


# 功能：建立当前项目的数据、中间对象、结果和日志目录。
def prepare_project_dirs(project):
    for path in (
        project.RAW_DIR,
        project.PROCESSED_DIR,
        project.REFERENCE_DIR,
        project.RESULTS_DIR,
        project.LOGS_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)
    make_report_dir(project.RESULTS_DIR)


# 功能：构建 运行与输出 所需对象或目录（定义：make_report_dir）。
def make_report_dir(out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in REPORT_STEPS:
        (out_dir / name).mkdir(parents=True, exist_ok=True)
    return out_dir


# 功能：保存 运行与输出 产生的对象、表格或图形（定义：save_figure）。
def save_figure(path_no_ext, fig=None, fig_size=None, dpi=300, close=True):
    path_no_ext = Path(path_no_ext)
    path_no_ext.parent.mkdir(parents=True, exist_ok=True)
    fig = fig or plt.gcf()
    if fig_size:
        width = fig_size.get("width", fig.get_figwidth())
        height = fig_size.get("height", fig.get_figheight())
        fig.set_size_inches(width, height)
    fig.tight_layout()
    fig.savefig(f"{path_no_ext}.png", dpi=dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(f"{path_no_ext}.pdf", facecolor="white", bbox_inches="tight")
    if close:
        plt.close(fig)


# 功能：保存 运行与输出 产生的对象、表格或图形（定义：save_fig）。
def save_fig(figure, out_dir, file_name, fig_size=None):
    fig = figure.gcf() if hasattr(figure, "gcf") else figure
    save_figure(Path(out_dir) / file_name, fig=fig, fig_size=fig_size)


# 功能：为 CON、EXP 及其他分组生成稳定配色（定义：group_palette）。
def group_palette(groups, configured=None):
    configured = configured or {}
    fallback = plt.get_cmap("tab10")
    return {
        group: configured.get(group, fallback(i % 10))
        for i, group in enumerate(groups)
    }


# 功能：提取 UMAP 的统一横纵坐标范围（定义：umap_limits）。
def umap_limits(adata):
    coords = np.asarray(adata.obsm["X_umap"])
    return (
        (coords[:, 0].min(), coords[:, 0].max()),
        (coords[:, 1].min(), coords[:, 1].max()),
    )


# 功能：固定 UMAP 比例与坐标范围，保证分面图可比较（定义：fix_umap_axes）。
def fix_umap_axes(axes, xlim, ylim):
    for ax in np.atleast_1d(axes).ravel():
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")
