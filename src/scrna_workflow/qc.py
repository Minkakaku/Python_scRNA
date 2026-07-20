from pathlib import Path

import numpy as np
import pandas as pd

from .runtime import save_figure


QC_COLUMNS = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
]


def add_qc_metrics(adata):
    import scanpy as sc

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    return adata


def _plot_qc_panel(adata, out_dir, prefix):
    import matplotlib.pyplot as plt
    import scanpy as sc

    sc.pl.violin(
        adata,
        QC_COLUMNS,
        jitter=False,
        groupby="sample",
        rotation=45,
        multi_panel=True,
        show=False,
    )
    for ax in plt.gcf().axes:
        ax.grid(False)
    save_figure(
        Path(out_dir) / f"{prefix}_nCount_nFeature_violin",
        fig_size={"width": 25, "height": 10},
    )

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        color="sample",
        legend_loc="none",
        ax=axes[0],
        show=False,
    )
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="sample",
        ax=axes[1],
        show=False,
    )
    for ax in axes:
        ax.grid(False)
    save_figure(Path(out_dir) / f"{prefix}_nCount_nFeature_scatter", fig)


def data_qc(_adata, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _adata = add_qc_metrics(_adata)
    _plot_qc_panel(_adata, out_dir, "raw")
    _adata.obs.to_csv(out_dir / "all_cell_features_raw.csv")

    for sample in _adata.obs["sample"].drop_duplicates():
        sample_dir = out_dir / str(sample)
        sample_dir.mkdir(parents=True, exist_ok=True)
        tmp = _adata[_adata.obs["sample"] == sample].copy()
        _plot_qc_panel(tmp, sample_dir, "raw")
        tmp.obs.to_csv(sample_dir / "cell_features_raw.csv")
    return _adata


def cell_filtering(
    _adata,
    min_genes,
    min_cells,
    pct_counts_mt,
    pct_counts_ribo,
    upper_lim,
    total_counts,
    expected_doublet_rate,
    threshold,
    out_dir,
    seed=1234,
):
    import matplotlib.pyplot as plt
    import scanpy as sc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata = _adata.copy()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    keep = (
        (adata.obs["n_genes_by_counts"] < upper_lim)
        & (adata.obs["pct_counts_mt"] < pct_counts_mt)
        & (adata.obs["pct_counts_ribo"] < pct_counts_ribo)
        & (adata.obs["total_counts"] < total_counts)
    )
    adata = adata[keep].copy()

    import seaborn as sns

    rate = expected_doublet_rate / 100 if expected_doublet_rate > 1 else expected_doublet_rate
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    score_tables = []
    thresholds = []
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample].copy()
        scrublet_kwargs = {
            "expected_doublet_rate": rate,
            "random_state": seed,
        }
        if threshold is not None:
            scrublet_kwargs["threshold"] = threshold
        sc.pp.scrublet(tmp, **scrublet_kwargs)
        adata.obs.loc[tmp.obs_names, "doublet_score"] = tmp.obs["doublet_score"]
        adata.obs.loc[tmp.obs_names, "predicted_doublet"] = tmp.obs["predicted_doublet"]
        score_tables.append(
            tmp.obs[["sample", "doublet_score", "predicted_doublet"]].copy()
        )
        thresholds.append(
            {
                "sample": sample,
                "threshold": tmp.uns.get("scrublet", {}).get("threshold", threshold),
            }
        )

    scores = pd.concat(score_tables)
    scores.to_csv(out_dir / "scrublet_scores.csv")
    pd.DataFrame(thresholds).to_csv(
        out_dir / "scrublet_thresholds.csv",
        index=False,
    )
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(
        data=scores,
        x="doublet_score",
        hue="sample",
        element="step",
        stat="density",
        common_norm=False,
        ax=ax,
    )
    ax.grid(False)
    save_figure(out_dir / "scrublet_score_distribution", fig)
    adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()

    adata.obs.to_csv(out_dir / "all_cell_features_filtered.csv")
    _plot_qc_panel(adata, out_dir, "filtered")
    for sample in adata.obs["sample"].drop_duplicates():
        sample_dir = out_dir / str(sample)
        sample_dir.mkdir(parents=True, exist_ok=True)
        tmp = adata[adata.obs["sample"] == sample].copy()
        _plot_qc_panel(tmp, sample_dir, "filtered")
        tmp.obs.to_csv(sample_dir / "cell_features_filtered.csv")

    retention = pd.DataFrame(
        {
            "before": _adata.obs["sample"].value_counts(),
            "after": adata.obs["sample"].value_counts(),
        }
    ).fillna(0)
    retention["retained_fraction"] = retention["after"] / retention["before"]
    retention.to_csv(out_dir / "cell_retention_by_sample.csv")
    plt.close("all")
    return adata
