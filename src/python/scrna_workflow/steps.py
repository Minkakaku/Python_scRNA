import pandas as pd

from .analysis import dim_reduction_cluster, find_marker_gene
from .annotation import celltype_annotation
from .io import load_10x_data
from .qc import cell_filtering, data_qc
from .runtime import make_report_dir, prepare_project_dirs
from .visualization import (
    plot_celltype_group_comparison,
    plot_celltype_proportion,
    visua_gene_express,
)


# 从样本表加载 10x 矩阵并保存原始合并 H5AD。
def step_load(project):
    prepare_project_dirs(project)
    report_dir = project.RESULTS_DIR
    adata = load_10x_data(
        data_dir=project.RAW_DIR,
        sample_info=project.SAMPLE_INFO,
        out_dir=report_dir,
        copy_cellranger_outputs=project.COPY_CELLRANGER_OUTPUTS,
    )
    project.PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    adata.write(project.ADATA_RAW)
    return adata


# 读取原始 H5AD，计算过滤前 QC 并保存带指标的检查点。
def step_qc_before_filter(project):
    import scanpy as sc

    make_report_dir(project.RESULTS_DIR)
    adata = sc.read_h5ad(project.ADATA_RAW)
    adata = data_qc(adata, project.RESULTS_DIR / "02-Data_QC")
    adata.write(project.ADATA_QC)
    return adata


# 使用项目 QC 参数过滤并逐样本识别双细胞，保存过滤后 H5AD。
def step_qc_filter(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_QC)
    adata = cell_filtering(
        adata,
        min_genes=project.MIN_GENES,
        min_cells=project.MIN_CELLS,
        pct_counts_mt=project.MAX_PCT_MT,
        pct_counts_ribo=project.MAX_PCT_RIBO,
        upper_lim=project.MAX_GENES,
        total_counts=project.MAX_COUNTS,
        expected_doublet_rate=project.EXPECTED_DOUBLET_RATE,
        threshold=project.SCRUBLET_THRESHOLD,
        out_dir=project.RESULTS_DIR / "03-Cell_filter",
        seed=project.RANDOM_SEED,
    )
    adata.write(project.ADATA_FILTERED)
    return adata


# 读取过滤后对象，按项目参数运行降维和聚类并保存 H5AD。
def step_global_clustering(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_FILTERED)
    adata = dim_reduction_cluster(
        adata,
        method=project.CLUSTER_METHOD,
        n_pcs=project.N_PCS,
        n_neighbors=project.N_NEIGHBORS,
        resolution=project.RESOLUTION,
        batch_effect=project.BATCH_EFFECT,
        out_dir=project.RESULTS_DIR / "04-Dimension_reduction_cluster",
        n_top_genes=project.N_TOP_GENES,
        regress_covariates=project.REGRESS_COVARIATES,
        seed=project.RANDOM_SEED,
    )
    adata.write(project.ADATA_CLUSTERED)
    return adata


# 读取聚类对象，按配置的检验方法和筛选阈值输出 cluster marker。
def step_find_markers(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_CLUSTERED)
    return find_marker_gene(
        adata,
        method=project.MARKER_METHOD,
        out_dir=project.RESULTS_DIR / "05-Marker_gene_analysis",
        logfc=project.MARKER_LOGFC,
        fdr=project.MARKER_FDR,
        pct=project.MARKER_MIN_PCT,
    )


# 加载人工映射或配置的自动注释方法，保存注释对象及实际使用的映射。
def step_cell_annotation(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_CLUSTERED)
    annotation = None
    if project.ANNOTATION_DATABASE.lower() == "manual":
        table = pd.read_csv(project.ANNOTATION_FILE, sep="\t", dtype=str).dropna()
        annotation = dict(zip(table["cluster"], table["cell_type"]))
    adata, mapping = celltype_annotation(
        adata,
        database=project.ANNOTATION_DATABASE,
        organism=project.ORGANISM,
        out_dir=project.RESULTS_DIR / "06-Cell_type_annotation",
        cell_type=annotation,
        celltypist_model=project.CELLTYPIST_MODEL,
        group_colors=project.GROUP_COLORS,
        group_order=project.GROUP_ORDER,
    )
    adata.write(project.ADATA_ANNOTATED)
    return adata, mapping


# 输出样本和组别组成；多组输入时同时计算样本层面的组间比较。
def step_cell_proportion(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    out_dir = project.RESULTS_DIR / "07-Cell_proportion"
    sample_table = plot_celltype_proportion(adata, out_dir, type="sample")
    group_table = plot_celltype_proportion(adata, out_dir, type="group")
    comparison = None
    if adata.obs["group"].nunique() > 1:
        comparison = plot_celltype_group_comparison(
            adata,
            out_dir,
            group_colors=project.GROUP_COLORS,
            group_order=project.GROUP_ORDER,
        )
    return sample_table, group_table, comparison


# 在已注释对象上绘制指定基因的 UMAP、小提琴图及配套表达表。
def step_visualize_genes(project, genes):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    return visua_gene_express(
        adata,
        genes=genes,
        out_dir=project.RESULTS_DIR / "08-Gene_visualization",
    )
