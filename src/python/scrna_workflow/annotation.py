from pathlib import Path

import pandas as pd

from .visualization import (
    plot_annotation_umaps,
    plot_celltype_group_comparison,
    plot_celltype_proportion,
)


# 检查人工映射覆盖所有 cluster，再写入 cell_type 分类列。
def apply_manual_annotation(adata, annotation_dict):
    clusters = adata.obs["cluster"].astype(str)
    mapping = {str(key): value for key, value in annotation_dict.items()}
    missing = sorted(set(clusters) - set(mapping))
    if missing:
        raise ValueError(f"Missing annotations for clusters: {', '.join(missing)}")
    adata.obs["cell_type"] = pd.Categorical(clusters.map(mapping))
    return adata, mapping


# 读取物种匹配的 PanglaoDB marker，以 ORA 富集汇总推断各 cluster 标签。
def _panglaodb_annotation(adata, organism):
    import decoupler as dc

    markers = dc.get_resource("PanglaoDB")
    organism_col = "human" if organism.lower() == "human" else "mouse"
    if organism_col in markers:
        markers = markers[markers[organism_col].astype(str).str.lower() == "true"]
    markers = markers.drop_duplicates(["cell_type", "genesymbol"])
    matrix = adata.raw.to_adata() if adata.raw is not None else adata
    dc.run_ora(
        mat=matrix,
        net=markers,
        source="cell_type",
        target="genesymbol",
        min_n=3,
        verbose=True,
        use_raw=False,
    )
    acts = dc.get_acts(matrix, obsm_key="ora_estimate")
    mean_enrichment = dc.summarize_acts(acts, groupby="cluster", min_std=1)
    mapping = {str(k): v for k, v in dc.assign_groups(mean_enrichment).items()}
    return apply_manual_annotation(adata, mapping)


# 使用指定 CellTypist 模型预测并在 cluster 内多数投票，按细胞索引回填标签。
def _celltypist_annotation(adata, model):
    import celltypist

    if not model:
        raise ValueError("CELLTYPIST_MODEL must be configured for CellTypist annotation")
    matrix = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    predictions = celltypist.annotate(
        matrix,
        model=model,
        majority_voting=True,
        over_clustering=adata.obs["cluster"],
    )
    labels = predictions.predicted_labels
    column = "majority_voting" if "majority_voting" in labels else labels.columns[0]
    adata.obs["cell_type"] = pd.Categorical(labels[column].reindex(adata.obs_names))
    mapping = (
        adata.obs.groupby("cluster", observed=True)["cell_type"]
        .agg(lambda x: "|".join(sorted(set(x.astype(str)))))
        .to_dict()
    )
    return adata, {str(k): v for k, v in mapping.items()}


# 选择人工、PanglaoDB 或 CellTypist 注释，输出映射、UMAP 和组成统计。
def celltype_annotation(
    _adata,
    database,
    organism,
    out_dir,
    cell_type=None,
    celltypist_model=None,
    group_colors=None,
    group_order=None,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    database_key = database.lower()

    if cell_type is not None or database_key == "manual":
        if cell_type is None:
            raise ValueError("Manual annotation requires a cluster-to-cell-type mapping")
        adata, annotation_dict = apply_manual_annotation(_adata, cell_type)
    elif database_key == "panglaodb":
        adata, annotation_dict = _panglaodb_annotation(_adata, organism)
    elif database_key == "celltypist":
        adata, annotation_dict = _celltypist_annotation(_adata, celltypist_model)
    else:
        raise ValueError("database must be manual, PanglaoDB, or CellTypist")

    pd.DataFrame(
        sorted(annotation_dict.items()),
        columns=["cluster", "cell_type"],
    ).to_csv(out_dir / "celltype_annotation_used.csv", index=False)
    adata.obs.to_csv(out_dir / "all_cells_features.csv")
    pd.crosstab(adata.obs["cell_type"], adata.obs["sample"]).to_csv(
        out_dir / "celltype_statistics.csv"
    )
    pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP_1", "UMAP_2"],
    ).to_csv(out_dir / "data_umap.csv")

    plot_annotation_umaps(adata, out_dir)
    plot_celltype_proportion(adata, out_dir, type="sample")
    plot_celltype_proportion(adata, out_dir, type="group")
    if adata.obs["group"].nunique() > 1:
        plot_celltype_group_comparison(
            adata,
            out_dir,
            group_colors=group_colors,
            group_order=group_order,
        )
    return adata, annotation_dict
