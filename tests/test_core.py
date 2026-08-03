import unittest
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scrna_workflow.analysis import cluster_proportion_table, detect_pca_knee
from scrna_workflow.annotation import apply_manual_annotation
from scrna_workflow.io import infer_group
from scrna_workflow.runtime import save_figure
from scrna_workflow.visualization import sample_proportion_table


class CoreWorkflowTests(unittest.TestCase):
    def setUp(self):
        obs = pd.DataFrame(
            {
                "sample": ["Con_1", "Con_1", "Exp_1", "Exp_1"],
                "group": ["CON", "CON", "EXP", "EXP"],
                "cluster": pd.Categorical(["0", "1", "0", "1"]),
            },
            index=["a", "b", "c", "d"],
        )
        self.adata = ad.AnnData(X=np.zeros((4, 2)), obs=obs)

    def test_group_inference(self):
        self.assertEqual(infer_group("Con_1"), "CON")
        self.assertEqual(infer_group("Exp-2"), "EXP")

    def test_pca_knee_detection(self):
        variance = [10, 6, 3.5, 2.2, 1.6, 1.3, 1.1, 1.0]
        self.assertEqual(detect_pca_knee(variance), 4)

    def test_manual_annotation_and_sample_level_proportions(self):
        adata, mapping = apply_manual_annotation(
            self.adata.copy(),
            {"0": "T cell", "1": "B cell"},
        )
        self.assertEqual(mapping["0"], "T cell")
        table = sample_proportion_table(adata)
        totals = table.groupby("sample")["proportion"].sum()
        self.assertTrue(np.allclose(totals, 1))

    def test_cluster_proportions_are_sample_level(self):
        table = cluster_proportion_table(self.adata)
        totals = table.groupby("sample")["proportion"].sum()
        self.assertTrue(np.allclose(totals, 1))

    def test_save_figure_writes_png_and_pdf(self):
        stem = Path(__file__).parent / "_plot_test"
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        try:
            save_figure(stem, fig)
            self.assertTrue(stem.with_suffix(".png").exists())
            self.assertTrue(stem.with_suffix(".pdf").exists())
        finally:
            stem.with_suffix(".png").unlink(missing_ok=True)
            stem.with_suffix(".pdf").unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
