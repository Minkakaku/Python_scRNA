import unittest
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scrna_workflow.analysis import cluster_proportion_table, detect_pca_knee
from scrna_workflow.annotation import apply_manual_annotation
from scrna_workflow.io import infer_group
from scrna_workflow.runtime import NPG_COLORS, categorical_palette, save_figure
from scrna_workflow.visualization import sample_proportion_table


# 验证样本分组、PCA 拐点、人工注释、样本比例和配对图输出。
class CoreWorkflowTests(unittest.TestCase):
    # 功能：构造最小 AnnData，供核心流程单元测试复用。
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

    # 功能：验证样本名中的 CON/EXP 分组可被正确识别。
    def test_group_inference(self):
        self.assertEqual(infer_group("Con_1"), "CON")
        self.assertEqual(infer_group("Exp-2"), "EXP")

    # 功能：验证 PCA 自动拐点算法返回预期 PC 数。
    def test_pca_knee_detection(self):
        variance = [10, 6, 3.5, 2.2, 1.6, 1.3, 1.1, 1.0]
        self.assertEqual(detect_pca_knee(variance), 4)

    # 功能：验证人工注释映射和样本内细胞比例归一化。
    def test_manual_annotation_and_sample_level_proportions(self):
        adata, mapping = apply_manual_annotation(
            self.adata.copy(),
            {"0": "T cell", "1": "B cell"},
        )
        self.assertEqual(mapping["0"], "T cell")
        table = sample_proportion_table(adata)
        totals = table.groupby("sample")["proportion"].sum()
        self.assertTrue(np.allclose(totals, 1))

    # 功能：验证聚类比例按样本计算且每个样本合计为 1。
    def test_cluster_proportions_are_sample_level(self):
        table = cluster_proportion_table(self.adata)
        totals = table.groupby("sample")["proportion"].sum()
        self.assertTrue(np.allclose(totals, 1))

    # 功能：验证最终图形始终成对写出 PNG 和 PDF。
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

    # 功能：验证分类图优先使用固定 NPG 顺序，且超过十类仍可扩展。
    def test_npg_palette_order_and_extension(self):
        labels = [f"C{x}" for x in range(12)]
        palette = categorical_palette(labels)
        self.assertEqual([palette[x] for x in labels[:10]], list(NPG_COLORS))
        self.assertEqual(len(set(palette.values())), 12)
        self.assertEqual(categorical_palette(["C0"], {"C0": "black"})["C0"], "black")
        # 核对实际 Matplotlib 图元的默认颜色，而非仅核对常量。
        fig, ax = plt.subplots()
        try:
            colors = [ax.plot([0, 1], [0, 1])[0].get_color() for _ in range(2)]
            self.assertEqual(colors, list(NPG_COLORS[:2]))
            self.assertEqual(ax.imshow(np.zeros((2, 2))).get_cmap().name, "npg_sequential")
        finally:
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
