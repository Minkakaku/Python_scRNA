from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu
from scipy.sparse import issparse

from .runtime import fix_umap_axes, group_palette, save_figure, umap_limits


FEATURE_COLORS = [
    "#FFEFD5",
    "#E6E6FA",
    "#87CEFA",
    "#6495ED",
    "#4169E1",
    "#0000CD",
    "#000080",
]


def sample_proportion_table(adata, celltype_col="cell_type"):
    wide = pd.crosstab(
        adata.obs["sample"],
        adata.obs[celltype_col],
        normalize="index",
    )
    long = wide.reset_index().melt(
        id_vars="sample",
        var_name=celltype_col,
        value_name="proportion",
    )
    groups = adata.obs[["sample", "group"]].drop_duplicates("sample")
    return long.merge(groups, on="sample", how="left")


def plot_celltype_proportion(_adata, out_dir, type="sample", fig_size=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    table = pd.crosstab(
        _adata.obs[type],
        _adata.obs["cell_type"],
        normalize="index",
    )
    table.to_csv(out_dir / f"celltype_proportion_{type}.csv")

    colors = _adata.uns.get("cell_type_colors", plt.get_cmap("tab20").colors)
    palette = dict(zip(table.columns, colors))
    ax = table.plot(
        kind="bar",
        stacked=True,
        color=[palette[x] for x in table.columns],
        edgecolor="black",
        linewidth=0.3,
        width=0.8,
        figsize=(
            fig_size["width"] if fig_size else max(6, len(table) * 1.2),
            fig_size["height"] if fig_size else 6,
        ),
    )
    ax.set_xlabel(type.capitalize())
    ax.set_ylabel("Proportion")
    ax.grid(False)
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        title="Cell types",
        frameon=False,
    )
    save_figure(out_dir / f"celltype_proportion_{type}", ax.figure)
    return table


def plot_celltype_group_comparison(
    adata,
    out_dir,
    group_colors=None,
    group_order=None,
):
    import seaborn as sns

    out_dir = Path(out_dir)
    table = sample_proportion_table(adata)
    table.to_csv(out_dir / "celltype_proportion_by_sample.csv", index=False)
    groups = group_order or table["group"].drop_duplicates().tolist()
    groups = [x for x in groups if x in set(table["group"])]
    palette = group_palette(groups, group_colors)
    celltype_order = table["cell_type"].drop_duplicates().tolist()

    tests = []
    if len(groups) == 2:
        for cell_type, data in table.groupby("cell_type", observed=True):
            x = data.loc[data["group"] == groups[0], "proportion"]
            y = data.loc[data["group"] == groups[1], "proportion"]
            p_value = mannwhitneyu(x, y).pvalue if len(x) and len(y) else np.nan
            tests.append(
                {
                    "cell_type": cell_type,
                    "group_1": groups[0],
                    "group_2": groups[1],
                    "p_value": p_value,
                    "significance": (
                        "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
                    ) if not pd.isna(p_value) else "NA",
                }
            )
    tests = pd.DataFrame(tests)
    tests.to_csv(out_dir / "celltype_proportion_group_tests.csv", index=False)

    fig, ax = plt.subplots(
        figsize=(max(8, table["cell_type"].nunique() * 0.8), 6)
    )
    sns.boxplot(
        data=table,
        x="cell_type",
        y="proportion",
        order=celltype_order,
        hue="group",
        hue_order=groups,
        palette=palette,
        showfliers=False,
        ax=ax,
    )
    sns.stripplot(
        data=table,
        x="cell_type",
        y="proportion",
        order=celltype_order,
        hue="group",
        hue_order=groups,
        palette=palette,
        dodge=True,
        legend=False,
        ax=ax,
    )
    if not tests.empty:
        labels = tests.set_index("cell_type")["significance"]
        global_max = table["proportion"].max()
        for x, cell_type in enumerate(celltype_order):
            local_max = table.loc[
                table["cell_type"] == cell_type,
                "proportion",
            ].max()
            ax.text(x, local_max + global_max * 0.06, labels[cell_type], ha="center")
        ax.set_ylim(top=global_max * 1.18 if global_max else 0.05)

    ax.set_xlabel("")
    ax.set_ylabel("Cell proportion")
    ax.tick_params(axis="x", rotation=75)
    ax.grid(False)
    ax.legend(title="", loc="upper left", bbox_to_anchor=(1.02, 1))
    save_figure(out_dir / "celltype_proportion_group_comparison", fig)
    return table, tests


def plot_annotation_umaps(adata, out_dir):
    import scanpy as sc

    out_dir = Path(out_dir)
    xlim, ylim = umap_limits(adata)
    sc.pl.umap(
        adata,
        color="cell_type",
        title="",
        legend_loc="on data",
        legend_fontsize=9,
        legend_fontoutline=2,
        show=False,
    )
    fix_umap_axes(plt.gca(), xlim, ylim)
    save_figure(out_dir / "umap_dimplot_celltype")

    levels = adata.obs["group"].drop_duplicates().tolist()
    if len(levels) > 1:
        fig, axes = plt.subplots(1, len(levels), figsize=(len(levels) * 5, 5), squeeze=False)
        for ax, group in zip(axes.ravel(), levels):
            sc.pl.umap(
                adata[adata.obs["group"] == group],
                color="cell_type",
                title=str(group),
                ax=ax,
                legend_loc="on data",
                legend_fontsize=10,
                legend_fontoutline=2,
                show=False,
            )
        fix_umap_axes(axes, xlim, ylim)
        save_figure(out_dir / "umap_split_group_celltype", fig)

    sc.pl.umap(
        adata,
        color=["cluster", "cell_type"],
        legend_loc="on data",
        legend_fontsize=10,
        legend_fontoutline=2,
        show=False,
    )
    for ax in plt.gcf().axes:
        if ax.get_xlabel().startswith("UMAP"):
            ax.set_aspect("equal")
    save_figure(out_dir / "umap_cluster_celltype")


def visua_gene_express(
    _adata,
    genes,
    out_dir,
    groupby=None,
    point_size=0.0001,
    vmin="p05",
    vmax="p95",
):
    import scanpy as sc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmap = LinearSegmentedColormap.from_list("hongf_feature", FEATURE_COLORS)
    xlim, ylim = umap_limits(_adata)
    use_raw = _adata.raw is not None
    sc.pl.umap(
        _adata,
        color=genes,
        size=point_size,
        ncols=3,
        cmap=cmap,
        use_raw=use_raw,
        vmin=vmin,
        vmax=vmax,
        show=False,
    )
    for ax in plt.gcf().axes:
        if ax.get_xlabel().startswith("UMAP"):
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_aspect("equal")
    stem = "_".join(genes)
    matrix = _adata.raw[:, genes].X if use_raw else _adata[:, genes].X
    values = matrix.toarray() if issparse(matrix) else np.asarray(matrix)
    metadata = [
        x
        for x in ("sample", "group", "cluster", "cell_type")
        if x in _adata.obs
    ]
    plot_data = _adata.obs[metadata].copy()
    for index, gene in enumerate(genes):
        plot_data[gene] = values[:, index]
    plot_data.to_csv(out_dir / f"{stem}_plot_data.csv")

    save_figure(out_dir / stem)

    groupby = groupby or ("cell_type" if "cell_type" in _adata.obs else "cluster")
    fig, axes = plt.subplots(len(genes), 1, figsize=(10, max(4, len(genes) * 4)), squeeze=False)
    for ax, gene in zip(axes.ravel(), genes):
        sc.pl.violin(
            _adata,
            gene,
            groupby=groupby,
            stripplot=False,
            inner="box",
            ax=ax,
            show=False,
            use_raw=use_raw,
        )
        ax.grid(False)
    save_figure(out_dir / f"{stem}_violin", fig)
    return out_dir / f"{stem}.png", out_dir / f"{stem}.pdf"
