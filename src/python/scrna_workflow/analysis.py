from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from .runtime import NPG_DIVERGING, fix_umap_axes, save_figure, set_adata_palette, style_dotplot, umap_limits


# 统一批次方法名称：缺省为 none，去除外围星号并转小写。
def _batch_method(value):
    return str(value or "none").strip("*").lower()


# 用 Kneedle 寻找递减凸形 PCA 方差曲线的拐点；无明确拐点时返回 None。
def detect_pca_knee(variance, sensitivity=1.0):
    """Return the Kneedle elbow for a decreasing PCA variance curve."""
    from kneed import KneeLocator

    variance = np.asarray(variance, dtype=float)
    pcs = np.arange(1, variance.size + 1)
    locator = KneeLocator(
        pcs,
        variance,
        curve="convex",
        direction="decreasing",
        S=sensitivity,
    )
    return int(locator.knee) if locator.knee is not None else None


# 保存 PC 方差表、自动拐点和标线图，并把选择信息写入 AnnData。
def save_pca_knee(adata, out_dir, sensitivity=1.0):
    """Save PCA variance, the automatic knee, and a marked elbow plot."""
    out_dir = Path(out_dir)
    variance = np.asarray(adata.uns["pca"]["variance"], dtype=float)
    out_dir.mkdir(parents=True, exist_ok=True)
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"], dtype=float)
    knee_pc = detect_pca_knee(variance, sensitivity=sensitivity)

    table = pd.DataFrame(
        {
            "pc": np.arange(1, variance.size + 1),
            "variance": variance,
            "variance_percent": variance_ratio * 100,
            "cumulative_percent": np.cumsum(variance_ratio) * 100,
        }
    )
    table["is_knee"] = table["pc"].eq(knee_pc)
    table.to_csv(out_dir / "pca_variance.csv", index=False)
    pd.DataFrame(
        [{
            "method": "Kneedle",
            "recommended_pc": knee_pc,
            "curve": "convex",
            "direction": "decreasing",
            "sensitivity": sensitivity,
        }]
    ).to_csv(out_dir / "pca_knee_selection.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(
        table["pc"],
        table["variance_percent"],
        color="#4D4D4D",
        marker="o",
        markersize=3,
        linewidth=1,
    )
    if knee_pc is not None:
        ax.axvline(
            knee_pc,
            color="#E64B35",
            linestyle="--",
            linewidth=1,
            label=f"Kneedle: PC{knee_pc}",
        )
        ax.legend(frameon=False)
    else:
        ax.text(
            0.98,
            0.95,
            "No clear knee identified",
            ha="right",
            transform=ax.transAxes,
        )
    ax.set(
        xlabel="Principal component",
        ylabel="Explained variance (%)",
        title="PCA variance and automatic elbow",
    )
    ax.grid(False)
    save_figure(out_dir / "pca_variance_ratio", fig)

    adata.uns["pca_selection"] = {
        "kneed_pc": knee_pc if knee_pc is not None else -1,
        "curve": "convex",
        "direction": "decreasing",
        "sensitivity": float(sensitivity),
    }
    print(f"Automatic PCA elbow (Kneedle): {knee_pc}")
    return knee_pc, table


# 功能：按样本汇总 cluster 或细胞类型比例（定义：cluster_proportion_table）。
def cluster_proportion_table(adata, cluster_col="cluster"):
    proportions = pd.crosstab(
        adata.obs["sample"],
        adata.obs[cluster_col],
        normalize="index",
    )
    long = proportions.reset_index().melt(
        id_vars="sample",
        var_name=cluster_col,
        value_name="proportion",
    )
    groups = adata.obs[["sample", "group"]].drop_duplicates("sample")
    return long.merge(groups, on="sample", how="left")


# 沿用全对象 UMAP 坐标按组或样本分面，并固定各面板坐标范围。
def _split_umap(adata, split_col, color, out_path):
    levels = adata.obs[split_col].drop_duplicates().tolist()
    if len(levels) < 2:
        return
    import scanpy as sc

    xlim, ylim = umap_limits(adata)
    fig, axes = plt.subplots(1, len(levels), figsize=(len(levels) * 5, 5), squeeze=False)
    for ax, level in zip(axes.ravel(), levels):
        sc.pl.umap(
            adata[adata.obs[split_col] == level],
            color=color,
            title=str(level),
            ax=ax,
            legend_loc="on data",
            legend_fontsize=10,
            legend_fontoutline=2,
            show=False,
        )
    fix_umap_axes(axes, xlim, ylim)
    save_figure(out_path, fig)


# 保留原始计数后依次归一化、筛选 HVG、PCA、批次处理、聚类和 UMAP。
def dim_reduction_cluster(
    _adata,
    method,
    n_pcs,
    n_neighbors,
    resolution,
    batch_effect,
    out_dir,
    n_top_genes=2000,
    regress_covariates=("total_counts", "pct_counts_mt"),
    seed=1234,
):
    import scanpy as sc
    import scanpy.external as sce

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata = _adata.copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # 在筛选 HVG 和缩放前保存全基因 log1p 表达，供后续 marker 与基因绘图使用。
    adata.raw = adata.copy()

    hvg_kwargs = {
        "min_mean": 0.0125,
        "max_mean": 3,
        "min_disp": 0.5,
        "n_top_genes": n_top_genes,
        "flavor": "seurat",
    }
    if adata.obs["sample"].nunique() > 1:
        hvg_kwargs["batch_key"] = "sample"
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)
    sc.pl.highly_variable_genes(adata, show=False)
    save_figure(
        out_dir / "highly_variable_genes",
        fig_size={"width": 16, "height": 8},
    )

    adata = adata[:, adata.var["highly_variable"]].copy()
    covariates = [x for x in regress_covariates if x in adata.obs]
    if covariates:
        sc.pp.regress_out(adata, covariates)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", random_state=seed)

    sc.pl.pca(adata, color="sample", show=False)
    save_figure(out_dir / "pca_sample")
    sc.pl.pca(adata, color="group", show=False)
    save_figure(out_dir / "pca_group")
    save_pca_knee(adata, out_dir)

    # Harmony 校正 PCA，BBKNN 调整邻接图；其他值使用未校正的 PCA 邻接图。
    batch_method = _batch_method(batch_effect)
    if batch_method == "harmony":
        sce.pp.harmony_integrate(
            adata,
            key="sample",
            basis="X_pca",
            adjusted_basis="X_pca_harmony",
        )
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            use_rep="X_pca_harmony",
            random_state=seed,
        )
    elif batch_method == "bbknn":
        sce.pp.bbknn(adata, batch_key="sample", n_pcs=n_pcs)
    else:
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            random_state=seed,
        )

    sc.tl.umap(adata, random_state=seed)
    method = method.lower()
    if method == "leiden":
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added="cluster",
            random_state=seed,
        )
    elif method == "louvain":
        sc.tl.louvain(
            adata,
            resolution=resolution,
            key_added="cluster",
            random_state=seed,
        )
    else:
        raise ValueError("method must be leiden or louvain")

    xlim, ylim = umap_limits(adata)
    for color, name, title in (
        ("cluster", "umap_dimplot", "Cluster"),
        ("sample", "umap_sample_dimplot", ""),
        ("group", "umap_group_dimplot", ""),
    ):
        set_adata_palette(adata, color)
        sc.pl.umap(
            adata,
            color=color,
            title=title,
            legend_loc="on data" if color == "cluster" else "right margin",
            legend_fontsize=10,
            legend_fontoutline=2,
            show=False,
        )
        fix_umap_axes(plt.gca(), xlim, ylim)
        save_figure(out_dir / name)

    metrics = [
        x
        for x in ("cluster", "total_counts", "n_genes_by_counts", "total_counts_mt")
        if x == "cluster" or x in adata.obs
    ]
    sc.pl.umap(
        adata,
        color=metrics,
        ncols=2,
        cmap=LinearSegmentedColormap.from_list("npg_diverging", NPG_DIVERGING),
        use_raw=False,
        show=False,
    )
    for ax in plt.gcf().axes:
        if ax.get_xlabel().startswith("UMAP"):
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_aspect("equal")
    save_figure(
        out_dir / "cluster_total_counts_n_genes_by_counts",
        fig_size={"width": 8, "height": 8},
    )

    _split_umap(adata, "group", "cluster", out_dir / "umap_split_group_cluster")
    _split_umap(adata, "sample", "cluster", out_dir / "umap_split_sample_cluster")

    adata.obs.to_csv(out_dir / "all_cells_features.csv")
    pd.crosstab(adata.obs["cluster"], adata.obs["sample"]).to_csv(
        out_dir / "cluster_cells_statistics.csv"
    )
    cluster_proportion_table(adata).to_csv(
        out_dir / "cluster_proportion_by_sample.csv",
        index=False,
    )
    pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP_1", "UMAP_2"],
    ).to_csv(out_dir / "data_umap.csv")

    sc.pl.correlation_matrix(adata, "cluster", show=False)
    save_figure(out_dir / "cluster_correlation")
    return adata


# 保留分类列的既有 cluster 顺序；普通字符串列则排序后返回。
def _cluster_categories(adata):
    values = adata.obs["cluster"]
    if isinstance(values.dtype, pd.CategoricalDtype):
        return values.cat.categories.tolist()
    return sorted(values.astype(str).unique())


# 按 cluster 检验 marker，依 logFC、FDR 和表达比例筛选并输出表图。
def find_marker_gene(
    _adata,
    method,
    out_dir,
    logfc=1,
    fdr=0.05,
    pct=0.25,
):
    import scanpy as sc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sc.tl.rank_genes_groups(
        _adata,
        groupby="cluster",
        method=method,
        pts=True,
        use_raw=_adata.raw is not None,
    )
    sc.pl.rank_genes_groups(_adata, n_genes=25, sharey=False, show=False)
    save_figure(
        out_dir / "rank_genes_groups",
        fig_size={"width": 20, "height": 12},
    )

    all_diff = []
    for cluster in _cluster_categories(_adata):
        cluster_dir = out_dir / f"cluster_{cluster}"
        cluster_dir.mkdir(parents=True, exist_ok=True)
        markers = sc.get.rank_genes_groups_df(_adata, group=cluster)
        markers = markers[
            (markers["logfoldchanges"].abs() > logfc)
            & (markers["pvals_adj"] < fdr)
            & (markers["pct_nz_group"] >= pct)
        ].copy()
        markers.to_csv(cluster_dir / "differentially_expressed_genes.csv", index=False)
        markers["cluster"] = str(cluster)
        all_diff.append(markers)

        top_genes = markers["names"].head(5).tolist()
        if top_genes:
            sc.pl.rank_genes_groups_violin(
                _adata,
                groups=cluster,
                gene_names=top_genes,
                jitter=False,
                size=0,
                show=False,
            )
            save_figure(cluster_dir / "top_5_genes_violin")

    all_markers = pd.concat(all_diff, ignore_index=True)
    all_markers = all_markers.rename(columns={"names": "Symbol"})
    all_markers.to_csv(out_dir / "all_cluster_marker.csv", index=False)
    all_markers["cluster"].value_counts().sort_index().to_csv(
        out_dir / "cluster_diff_statistics.csv",
        header=["marker_count"],
    )

    sc.pl.rank_genes_groups_dotplot(
        _adata,
        n_genes=3,
        values_to_plot="logfoldchanges",
        vmin=-3,
        vmax=3,
        min_logfoldchange=logfc,
        cmap=LinearSegmentedColormap.from_list("npg_diverging", NPG_DIVERGING),
        smallest_dot=18,
        colorbar_title="log fold change",
        show=False,
    )
    style_dotplot()
    save_figure(
        out_dir / "dotplot_top3_marker_gene",
        fig_size={"width": 20, "height": 8},
    )

    sc.pl.rank_genes_groups_heatmap(
        _adata,
        n_genes=10,
        values_to_plot="logfoldchanges",
        min_logfoldchange=logfc,
        swap_axes=True,
        cmap=LinearSegmentedColormap.from_list("npg_diverging", NPG_DIVERGING),
        standard_scale="var",
        show_gene_labels=True,
        dendrogram=False,
        show=False,
    )
    save_figure(
        out_dir / "heatmap_top10_marker_gene",
        fig_size={"width": 26, "height": 20},
    )
    return _adata
