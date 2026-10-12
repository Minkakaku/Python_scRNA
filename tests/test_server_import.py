import gzip
import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

from scrna_workflow.compatibility import convert_h5ad_legacy
from scrna_workflow.io import load_expression_data, load_10x_data
from scrna_workflow.vdj import clonotype_tables, vdj_cell_lookup


# 写入实际压缩 MEX，用不同 feature 列布局区分 C4/BD 与 10x。
def write_mex(path, counts, genes, barcodes, platform):
    path.mkdir(parents=True)
    with gzip.open(path / "matrix.mtx.gz", "wb") as handle:
        mmwrite(handle, csr_matrix(counts).T)
    features = [f"ID{i}\t{gene}\tGene Expression" if platform == "10x" else gene for i, gene in enumerate(genes)]
    with gzip.open(path / "features.tsv.gz", "wt") as handle:
        handle.write("\n".join(features) + "\n")
    with gzip.open(path / "barcodes.tsv.gz", "wt") as handle:
        handle.write("\n".join(barcodes) + "\n")


# 真实 C4/BD 混合读入，验证前缀、样本隔离、共有基因及计数数值。
def test_platform_mixed_c4_bd_inner_outer(tmp_path):
    write_mex(tmp_path / "c4/output/filter_matrix", [[1, 2], [3, 4]], ["A", "B"], ["aa", "bb"], "c4")
    write_mex(tmp_path / "bd/filtered_feature_bc_matrix", [[5, 6], [7, 8]], ["B", "C"], ["1", "CellLabel_2"], "bd")
    info = pd.DataFrame({"raw": ["c4", "bd"], "sample": ["C4", "BD"], "group": ["CON", "EXP"], "platform": ["c4", "bd"]})
    inner = load_expression_data(tmp_path, info)
    assert inner.var_names.tolist() == ["B"]
    assert inner.X.toarray().ravel().tolist() == [2, 4, 5, 7]
    assert inner.obs_names.tolist() == ["C4:aa", "C4:bb", "BD:CellLabel_1", "BD:CellLabel_2"]
    outer = load_expression_data(tmp_path, info, join="outer")
    assert outer.var_names.tolist() == ["A", "B", "C"]
    np.testing.assert_array_equal(outer.X.toarray(), [[1, 2, 0], [3, 4, 0], [0, 5, 6], [0, 7, 8]])
    assert list(outer.obs["sample"].cat.categories) == ["C4", "BD"]


# 真正调用 Scanpy 读取 10x/寻因；旧接口仍按 10x 第二列读取 symbol。
def test_platform_10x_and_seek_keep_old_interface(tmp_path):
    pytest.importorskip("scanpy")
    write_mex(tmp_path / "ten/outs/filtered_feature_bc_matrix", [[1, 2]], ["A", "B"], ["same"], "10x")
    write_mex(tmp_path / "seek/Outs/filtered_feature_bc_matrix", [[3, 4]], ["A", "B"], ["same"], "10x")
    info = pd.DataFrame({"raw": ["ten", "seek"], "sample": ["S1", "S2"], "group": ["CON", "EXP"], "platform": ["10x", "seek"]})
    combined = load_expression_data(tmp_path, info, cache=False)
    np.testing.assert_array_equal(combined.X.toarray(), [[1, 2], [3, 4]])
    assert combined.obs_names.tolist() == ["S1:same", "S2:same"]
    old = load_10x_data(tmp_path, info.iloc[:1].drop(columns="platform"), cache=False)
    assert old.var_names.tolist() == ["A", "B"]


# 数据入口不能静默合并重复样本，也不能凭空解释未知平台。
def test_platform_rejects_duplicate_samples(tmp_path):
    info = pd.DataFrame({"raw": ["a", "b"], "sample": ["S", "S"], "group": ["CON", "EXP"]})
    with pytest.raises(ValueError, match="unique"):
        load_expression_data(tmp_path, info)


# H5AD 兼容导出必须保留原文件、稀疏计数、raw、layer、细胞基因顺序与 UMAP。
def test_legacy_conversion_preserves_biological_arrays(tmp_path):
    obs = pd.DataFrame({"cell_type": pd.Categorical(["T", "B"]), "nullable": pd.array([1, None], dtype="Int64")}, index=["cell2", "cell1"])
    obj = ad.AnnData(csr_matrix([[2., 0.], [1., 3.]]), obs=obs, var=pd.DataFrame(index=["G2", "G1"]))
    obj.layers["counts"] = obj.X.copy()
    obj.raw = obj.copy()
    obj.obsm["X_umap"] = np.array([[1., 3.], [2., 4.]])
    obj.uns["nested"] = {"value": 1, "optional": None}
    src, dst = tmp_path / "source.h5ad", tmp_path / "legacy.h5ad"
    obj.write_h5ad(src)
    digest = hashlib.sha256(src.read_bytes()).hexdigest()
    convert_h5ad_legacy(src, dst)
    result = ad.read_h5ad(dst)
    assert hashlib.sha256(src.read_bytes()).hexdigest() == digest
    assert result.obs_names.tolist() == obj.obs_names.tolist()
    assert result.var_names.tolist() == obj.var_names.tolist()
    assert result.obs["cell_type"].tolist() == obj.obs["cell_type"].tolist()
    np.testing.assert_array_equal(result.X.toarray(), obj.X.toarray())
    np.testing.assert_array_equal(result.raw.X.toarray(), obj.raw.X.toarray())
    np.testing.assert_array_equal(result.layers["counts"].toarray(), obj.layers["counts"].toarray())
    np.testing.assert_array_equal(result.obsm["X_umap"], obj.obsm["X_umap"])
    assert "optional" not in result.uns["nested"]
    with pytest.raises(FileExistsError):
        convert_h5ad_legacy(src, dst)
    with pytest.raises(ValueError, match="different"):
        convert_h5ad_legacy(src, src, overwrite=True)


# 旧环境校验失败时已有输出必须完整保留，临时文件须清理。
def test_legacy_failed_validation_preserves_destination(tmp_path, monkeypatch):
    import sys
    import scrna_workflow.compatibility as compatibility

    src, dst = tmp_path / "source.h5ad", tmp_path / "existing.h5ad"
    ad.AnnData(np.ones((2, 2))).write_h5ad(src)
    dst.write_bytes(b"existing output")
    def fail(*args):
        raise RuntimeError("validation failed")
    monkeypatch.setattr(compatibility, "validate_with_legacy", fail)
    with pytest.raises(RuntimeError, match="validation failed"):
        convert_h5ad_legacy(src, dst, overwrite=True, legacy_python=sys.executable)
    assert dst.read_bytes() == b"existing output"
    assert not list(tmp_path.glob(".existing.h5ad.*"))


# 两个样本可以使用同一 barcode，同一样本内重复才是歧义。
def test_vdj_lookup_uses_full_barcode_and_sample():
    obs = pd.DataFrame({"sample": ["Con-1", "Exp-1"], "barcode": ["AAAAAAAAAAAAAAAAAAAA-1"] * 2}, index=["cell_a", "cell_b"])
    gex = ad.AnnData(np.ones((2, 1)), obs=obs)
    lookup = vdj_cell_lookup(gex)
    assert lookup.loc[("Con-1", "AAAAAAAAAAAAAAAAAAAA-1")] == "cell_a"
    assert lookup.loc[("Exp-1", "AAAAAAAAAAAAAAAAAAAA-1")] == "cell_b"
    gex.obs["sample"] = "same"
    with pytest.raises(ValueError, match="unique"):
        vdj_cell_lookup(gex)


# 相同克隆簇在两个样本独立计大小，Jaccard 和多样性以样本为单位。
def test_vdj_clonotype_tables_are_sample_level():
    obs = pd.DataFrame({"sample": ["S1", "S1", "S1", "S2"], "group": ["CON"] * 3 + ["EXP"],
                        "cell_type": ["T"] * 4, "clonotype_cluster": ["A", "A", "B", "A"]})
    cells, counts, diversity, overlap, expansion = clonotype_tables(obs)
    assert cells["clone_size_in_sample"].tolist() == [2, 2, 1, 1]
    assert overlap.loc["S1", "S2"] == 0.5
    assert diversity.set_index("sample").loc["S2", "normalized_shannon_entropy"] == 0
    assert expansion.loc["S1", "2"] == 2
    assert counts["cell_count"].sum() == 4


# 真实运行 Wilcoxon，检查细胞类型正向 marker 能从 raw 全基因中取出。
def test_celltype_markers_use_raw_and_positive_effects(tmp_path, monkeypatch):
    sc = pytest.importorskip("scanpy")
    import scrna_workflow.analysis as analysis

    values = np.array([[10, 0, 2]] * 6 + [[0, 10, 2]] * 6, dtype=float)
    obj = ad.AnnData(values, obs=pd.DataFrame({"cell_type": pd.Categorical(["T/CD4"] * 6 + ["B"] * 6)}, index=[f"c{i}" for i in range(12)]), var=pd.DataFrame(index=["T_marker", "B_marker", "shared"]))
    sc.pp.normalize_total(obj, target_sum=10000)
    sc.pp.log1p(obj)
    obj.raw = obj.copy()
    obj = obj[:, ["shared"]].copy()
    # 本测试检验真实统计与输出表；现有绘图已由核心测试覆盖。
    for name in ["rank_genes_groups", "rank_genes_groups_violin", "rank_genes_groups_dotplot", "rank_genes_groups_heatmap"]:
        monkeypatch.setattr(sc.pl, name, lambda *args, **kwargs: None)
    monkeypatch.setattr(analysis, "save_figure", lambda *args, **kwargs: None)
    monkeypatch.setattr(analysis, "style_dotplot", lambda *args, **kwargs: None)
    analysis.find_marker_gene(obj, "wilcoxon", tmp_path, logfc=0.25, fdr=0.05, pct=0.1, groupby="cell_type", positive_only=True)
    result = pd.read_csv(tmp_path / "all_cell_type_marker.csv")
    assert set(result["Symbol"]) == {"T_marker", "B_marker"}
    assert (result["logfoldchanges"] > 0.25).all()
    assert set(result["cell_type"]) == {"T/CD4", "B"}
    tested = pd.read_csv(tmp_path / "celltype_01/all_tested_genes.csv")
    assert set(tested["Symbol"]) == {"T_marker", "B_marker", "shared"}
    assert tested["cell_type"].nunique() == 1
    assert {"pvals", "pvals_adj", "scores", "logfoldchanges"}.issubset(tested.columns)


# 实际图元核对：统一范围、边界点留白、单位比例及相关面板尺寸一致。
def test_shared_umap_limits_keep_boundary_cells_and_equal_units():
    import matplotlib.pyplot as plt
    from scrna_workflow.runtime import umap_limits, fix_umap_axes

    obj = ad.AnnData(np.ones((2, 1)))
    obj.obsm["X_umap"] = np.array([[-5., -1.], [3., 4.]])
    limits = umap_limits(obj)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    try:
        fix_umap_axes(axes, *limits)
        fig.canvas.draw()
        for ax in axes:
            points = ax.transData.transform(obj.obsm["X_umap"])
            assert ((points[:, 0] > ax.bbox.x0) & (points[:, 0] < ax.bbox.x1)).all()
            assert ((points[:, 1] > ax.bbox.y0) & (points[:, 1] < ax.bbox.y1)).all()
            origin, unit_x, unit_y = ax.transData.transform([[0, 0], [1, 0], [0, 1]])
            assert np.isclose(np.linalg.norm(unit_x - origin), np.linalg.norm(unit_y - origin))
        assert np.isclose(axes[0].bbox.width, axes[1].bbox.width)
        assert np.isclose(axes[0].bbox.height, axes[1].bbox.height)
    finally:
        plt.close(fig)
