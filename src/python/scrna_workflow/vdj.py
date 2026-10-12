"""提取公共 VDJ 流程的链质控、克隆簇、扩增、共有比例与绘图；复用已有 GEX 坐标。"""
from pathlib import Path

import numpy as np
import pandas as pd

from .runtime import categorical_palette, fix_umap_axes, save_figure, umap_limits


# 用样本和完整 barcode 对齐，不截取前 18 位，也不依赖目录名正则。
def vdj_cell_lookup(gex):
    if not {"sample", "barcode"}.issubset(gex.obs.columns):
        raise ValueError("VDJ matching requires explicit sample and barcode columns in GEX obs")
    keys = gex.obs[["sample", "barcode"]]
    if keys.isna().any().any() or keys.astype(str).eq("").any().any():
        raise ValueError("GEX sample/barcode identifiers must not be empty")
    if not gex.obs_names.is_unique or keys.duplicated().any():
        raise ValueError("GEX sample/barcode pairs and cell identifiers must be unique")
    return pd.Series(gex.obs_names, index=pd.MultiIndex.from_frame(keys.astype(str)))


# 从显式 VDJ 样本 TSV 读取 contig；只保留与 GEX 精确匹配的细胞。
def load_vdj_data(gex, sample_info, base_dir=None):
    import anndata as ad
    import scirpy as ir

    info = pd.read_csv(sample_info, sep="\t", dtype=str) if isinstance(sample_info, (str, Path)) else sample_info.copy()
    if info.empty or not {"sample", "contigs"}.issubset(info):
        raise ValueError("VDJ sample table requires sample and contigs columns")
    if info[["sample", "contigs"]].isna().any().any() or info[["sample", "contigs"]].eq("").any().any():
        raise ValueError("VDJ sample/contigs values must not be empty")
    if info["sample"].duplicated().any():
        raise ValueError("VDJ sample identifiers must be unique")
    base = Path(base_dir or ".")
    lookup = vdj_cell_lookup(gex)
    objects = []
    for record in info.itertuples(index=False):
        airr = ir.io.read_10x_vdj(base / record.contigs)
        keys = pd.MultiIndex.from_arrays([[record.sample] * airr.n_obs, airr.obs_names.astype(str)])
        matched = lookup.reindex(keys)
        keep = matched.notna().values
        if not keep.any():
            raise ValueError(f"No GEX/VDJ cells match for sample {record.sample}")
        airr = airr[keep].copy()
        airr.obs_names = pd.Index(matched[keep].values)
        objects.append(airr)
    airr = ad.concat(objects, join="outer", index_unique=None, merge="same")
    if not airr.obs_names.is_unique:
        raise ValueError("VDJ cell identifiers must be unique")
    # 按原 GEX 顺序重排 AIRR；表达、注释、UMAP 和原对象均不改写。
    cells = gex.obs_names[gex.obs_names.isin(airr.obs_names)]
    airr = airr[cells].copy()
    for key in ("sample", "group", "cell_type", "barcode"):
        if key not in gex.obs:
            raise ValueError(f"GEX is missing required metadata: {key}")
        airr.obs[key] = gex.obs.loc[cells, key].values
    return gex[cells].copy(), airr


# 克隆大小按每个样本独立统计；相似序列克隆簇不等同于严格同一克隆。
def clonotype_tables(obs):
    required = ["sample", "group", "cell_type", "clonotype_cluster"]
    if obs[required].isna().any().any():
        raise ValueError("Clonotype summaries require non-missing sample/group/type/cluster labels")
    if (obs.groupby("sample", observed=True)["group"].nunique() > 1).any():
        raise ValueError("Each biological sample must belong to one group")
    cells = obs.copy()
    cells["clone_size_in_sample"] = cells.groupby(["sample", "clonotype_cluster"], observed=True)["sample"].transform("size")
    counts = cells.groupby(required, observed=True).size().reset_index(name="cell_count")
    diversity = []
    for sample, part in cells.groupby("sample", observed=True):
        sizes = part.groupby("clonotype_cluster", observed=True).size()
        proportions = sizes / sizes.sum()
        entropy = float(-(proportions * np.log(proportions)).sum())
        diversity.append({"sample": sample, "group": part["group"].iloc[0], "vdj_cells": len(part),
                          "clonotype_clusters": len(sizes), "shannon_entropy": entropy,
                          "normalized_shannon_entropy": entropy / np.log(len(sizes)) if len(sizes) > 1 else 0.0})
    # 样本间 Jaccard 基于克隆簇是否存在；不把细胞当作重复进行条件检验。
    presence = pd.crosstab(cells["sample"], cells["clonotype_cluster"]).gt(0).astype(int)
    intersection = presence.values @ presence.values.T
    sizes = presence.sum(axis=1).values
    union = sizes[:, None] + sizes[None, :] - intersection
    overlap = pd.DataFrame(intersection / union, index=presence.index, columns=presence.index)
    cells["expansion_bin"] = np.minimum(cells["clone_size_in_sample"], 6).map(lambda x: "6+" if x == 6 else str(x))
    expansion = pd.crosstab(cells["sample"], cells["expansion_bin"]).reindex(columns=["1", "2", "3", "4", "5", "6+"], fill_value=0)
    return cells, counts, pd.DataFrame(diversity), overlap, expansion


# UMAP 使用原全对象范围与固定面板尺寸；类别图和连续图都使用 NPG。
def _plot_vdj_umaps(airr, coords, limits, out_dir):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    for key in ("sample", "group", "cell_type", "clone_size_in_sample"):
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_axes([0.12, 0.12, 0.66, 0.78])
        if key == "clone_size_in_sample":
            points = ax.scatter(coords[:, 0], coords[:, 1], c=airr.obs[key], s=5, linewidths=0,
                                cmap="npg_sequential", rasterized=False)
            fig.colorbar(points, cax=fig.add_axes([0.84, 0.2, 0.025, 0.55]), label="Clone size in sample")
        else:
            values = airr.obs[key].astype(str)
            palette = categorical_palette(sorted(values.unique()))
            ax.scatter(coords[:, 0], coords[:, 1], c=values.map(palette), s=5, linewidths=0, rasterized=False)
            handles = [Line2D([], [], marker="o", linestyle="", color=color, label=label, markersize=4)
                       for label, color in palette.items()]
            ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
        ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
        fix_umap_axes(ax, *limits)
        # 固定 axes 位置的画布无需 tight_layout，避免图例改变面板大小。
        fig.savefig(out_dir / f"umap_{key}.pdf", bbox_inches="tight", facecolor="white")
        fig.savefig(out_dir / f"umap_{key}.png", dpi=600, bbox_inches="tight", facecolor="white")
        plt.close(fig)


# 标准结果图均从已保存的统计表绘制，避免重复计算或额外生成审计文件。
def _plot_vdj_tables(diversity, overlap, expansion, out_dir):
    import matplotlib.pyplot as plt

    colors = list(categorical_palette(expansion.columns).values())
    ax = expansion.plot.bar(stacked=True, color=colors, figsize=(8, 5))
    ax.set(ylabel="VDJ cells", xlabel="Sample")
    ax.grid(False, which="both")
    ax.legend(title="Clone size in sample", frameon=False)
    save_figure(out_dir / "clonal_expansion", ax.figure)
    ax = diversity.plot.bar(x="sample", y="normalized_shannon_entropy", legend=False,
                            color=list(categorical_palette(diversity["sample"]).values()), figsize=(8, 5))
    ax.set(ylabel="Normalized Shannon entropy", xlabel="Sample")
    ax.grid(False, which="both")
    save_figure(out_dir / "sample_diversity", ax.figure)
    fig, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(overlap, cmap="npg_sequential", vmin=0, vmax=1)
    ax.set_xticks(range(len(overlap)), overlap.columns, rotation=90)
    ax.set_yticks(range(len(overlap)), overlap.index)
    ax.grid(False, which="both")
    fig.colorbar(image, ax=ax, label="Jaccard overlap of clonotype clusters")
    save_figure(out_dir / "repertoire_overlap", fig)


# 复用公共 Scirpy 算法；TCR/BCR 共用入口，不启动 RNA 重聚类或自动抗原推断。
def run_vdj_analysis(gex, sample_info, out_dir, object_path, base_dir=None,
                     method="TCR", metric="alignment", cutoff=15, min_network_cells=3, seed=1234):
    import scirpy as ir
    from mudata import MuData

    method = str(method).upper()
    if method not in {"TCR", "BCR"} or metric not in {"alignment", "identity"}:
        raise ValueError("VDJ method must be TCR/BCR and metric alignment/identity")
    if cutoff < 0 or min_network_cells < 1:
        raise ValueError("VDJ cutoff and minimum network cells are invalid")
    object_path = Path(object_path)
    if object_path.exists():
        raise FileExistsError(f"Preserve the existing VDJ object; choose a new path: {object_path}")
    limits = umap_limits(gex) if "X_umap" in gex.obsm else None
    gex, airr = load_vdj_data(gex, sample_info, base_dir)
    ir.pp.index_chains(airr)
    ir.tl.chain_qc(airr)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    airr.obs.to_csv(out_dir / "chain_qc.csv")
    keep = airr.obs["receptor_type"].eq(method) & ~airr.obs["chain_pairing"].isin(["multichain", "orphan VDJ", "orphan VJ"])
    if not keep.any():
        raise ValueError("No paired receptors remain after chain QC")
    airr = airr[keep].copy()
    gex = gex[airr.obs_names].copy()
    ir.pp.ir_dist(airr, sequence="aa", metric=metric, cutoff=cutoff if metric == "alignment" else None, n_jobs=1)
    ir.tl.define_clonotype_clusters(airr, sequence="aa", metric=metric, receptor_arms="all", dual_ir="any",
                                   key_added="clonotype_cluster", n_jobs=1)
    cells, counts, diversity, overlap, expansion = clonotype_tables(airr.obs)
    airr.obs = cells
    for name, table in [("vdj_cells", cells), ("clonotype_counts", counts), ("sample_diversity", diversity),
                        ("repertoire_overlap", overlap), ("clonal_expansion", expansion)]:
        table.to_csv(out_dir / f"{name}.csv", index=name not in {"clonotype_counts", "sample_diversity"})
    _plot_vdj_tables(diversity, overlap, expansion, out_dir)
    if limits is not None:
        _plot_vdj_umaps(airr, np.asarray(gex.obsm["X_umap"]), limits, out_dir)
    # 没有达到展示阈值的克隆簇时保留完整表，不强行生成空网络。
    if airr.obs["clonotype_cluster"].value_counts().max() >= min_network_cells:
        import matplotlib.pyplot as plt
        from .runtime import set_adata_palette

        ir.tl.clonotype_network(airr, sequence="aa", metric=metric, min_cells=min_network_cells,
                                clonotype_key="clonotype_cluster", random_state=seed)
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        for ax, key in zip(axes, ("sample", "group", "cell_type")):
            palette = set_adata_palette(airr, key)
            ir.pl.clonotype_network(airr, color=key, ax=ax, title="", palette=list(palette.values()),
                                    edges_color="#8491B4", cmap="npg_sequential", base_size=10, label_fontsize=8)
        save_figure(out_dir / "clonotype_network", fig)
    object_path.parent.mkdir(parents=True, exist_ok=True)
    mdata = MuData({"gex": gex, "airr": airr})
    mdata.write_h5mu(object_path)
    return mdata
