import numpy as np
import pandas as pd
import pytest
import anndata as ad


# 构造两个样本的真实 10x contig 表，包含相同 barcode 和可验证的共有克隆。
@pytest.mark.parametrize("metric", ["identity", "alignment"])
def test_vdj_real_scirpy_pipeline(tmp_path, metric):
    ir = pytest.importorskip("scirpy")
    pytest.importorskip("mudata")
    from scrna_workflow.vdj import run_vdj_analysis

    samples = ["Con-1", "Exp-1"]
    barcodes = ["AAAAAAAAAAAAAAAAAAAAA-1", "AAAAAAAAAAAAAAAAAAAAB-1"]
    rows, obs = [], []
    for sample in samples:
        path = tmp_path / sample / "contigs.csv"
        path.parent.mkdir()
        contigs = []
        for barcode in barcodes:
            # 三个细胞为 A 克隆；Exp-1 的第二个细胞为独立 B 克隆。
            alternate = sample == "Exp-1" and barcode == barcodes[1]
            obs.append({"cell_id": sample + ":" + barcode, "sample": sample, "barcode": barcode,
                        "group": "CON" if sample == "Con-1" else "EXP", "cell_type": "T cell"})
            for chain, sequence in [("TRA", "CAVRDTNTNAGKSTF"), ("TRB", "CASSLGQETQYF")]:
                contigs.append({"barcode": barcode, "is_cell": True, "contig_id": barcode + chain,
                                "high_confidence": True, "length": 400, "chain": chain,
                                "v_gene": chain + "V1", "d_gene": "TRBD1" if chain == "TRB" else "",
                                "j_gene": chain + "J1", "c_gene": chain + "C", "full_length": True,
                                "productive": True, "cdr3": "C" + "W" * (len(sequence) - 1) if alternate else sequence,
                                "cdr3_nt": "TGTGCTGTT" if alternate else "TGTGCTGCT", "reads": 20, "umis": 3})
        pd.DataFrame(contigs).to_csv(path, index=False)
        rows.append({"sample": sample, "contigs": str(path)})
    metadata = pd.DataFrame(obs).set_index("cell_id")
    gex = ad.AnnData(np.arange(8).reshape(4, 2).astype(float), obs=metadata,
                    var=pd.DataFrame(index=["G1", "G2"]))
    gex.obsm["X_umap"] = np.array([[-5., -1.], [3., 4.], [2., 3.], [1., 2.]])
    before = gex.X.copy()
    result = run_vdj_analysis(gex, pd.DataFrame(rows), tmp_path / "results", tmp_path / "object.h5mu", metric=metric, min_network_cells=1)
    assert result["airr"].n_obs == 4
    np.testing.assert_array_equal(gex.X, before)
    np.testing.assert_array_equal(result["gex"].obsm["X_umap"], gex.obsm["X_umap"])
    cells = pd.read_csv(tmp_path / "results/vdj_cells.csv", index_col=0)
    assert cells["clone_size_in_sample"].tolist() == [2, 2, 1, 1]
    assert cells["clonotype_cluster"].nunique() == 2
    assert (tmp_path / "results/clonotype_network.pdf").is_file()
    assert (tmp_path / "results/umap_cell_type.png").is_file()
    # 对象回读核对细胞及克隆簇，而非仅检查文件存在。
    from mudata import read_h5mu
    readback = read_h5mu(tmp_path / "object.h5mu")
    assert readback["airr"].obs_names.tolist() == gex.obs_names.tolist()
    assert readback["airr"].obs["clonotype_cluster"].astype(str).tolist() == result["airr"].obs["clonotype_cluster"].astype(str).tolist()
