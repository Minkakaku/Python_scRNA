import re
import shutil
from pathlib import Path

import pandas as pd


DEFAULT_GROUP_PATTERNS = {
    "CON": r"(?i)(^|[_-])con(?:trol)?([_-]|$)",
    "EXP": r"(?i)(^|[_-])exp(?:eriment)?([_-]|$)",
}


# 功能：根据样本名称或样本表推断实验分组（定义：infer_group）。
def infer_group(sample, patterns=None):
    patterns = patterns or DEFAULT_GROUP_PATTERNS
    for group, pattern in patterns.items():
        if re.search(pattern, str(sample)):
            return group
    return str(sample)


# 读取样本 TSV 并检查 raw/sample 列；缺失的分组按样本名规则补全。
def read_sample_info(path, group_patterns=None):
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "raw" not in table or "sample" not in table:
        raise ValueError("sample_info.tsv requires raw and sample columns")
    if "group" not in table:
        table["group"] = ""
    missing = table["group"].eq("")
    table.loc[missing, "group"] = table.loc[missing, "sample"].map(
        lambda x: infer_group(x, group_patterns)
    )
    return table


# 优先寻找 Cell Ranger 过滤 MTX 目录，其次寻找 H5；均缺失时中止。
def _matrix_path(raw_dir):
    outs = Path(raw_dir) / "outs"
    matrix_dir = outs / "filtered_feature_bc_matrix"
    matrix_h5 = outs / "filtered_feature_bc_matrix.h5"
    if matrix_dir.exists():
        return matrix_dir
    if matrix_h5.exists():
        return matrix_h5
    raise FileNotFoundError(f"No filtered Cell Ranger matrix under {outs}")


# 读取单样本矩阵，并用 sample:barcode 构造跨样本唯一细胞索引。
def _read_10x_record(record, data_dir, cache=True):
    import scanpy as sc

    matrix = _matrix_path(Path(data_dir) / record.raw)
    if matrix.is_dir():
        adata = sc.read_10x_mtx(matrix, var_names="gene_symbols", cache=cache)
    else:
        adata = sc.read_10x_h5(matrix)
    adata.var_names_make_unique()
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs["sample"] = record.sample
    adata.obs["group"] = record.group
    adata.obs_names = pd.Index(
        record.sample + ":" + adata.obs["barcode"].astype(str),
        name="cell_id",
    )
    return adata


# 复制存在的 Cell Ranger 报告、指标和过滤矩阵到当前报告目录。
def _copy_cellranger_outputs(record, data_dir, out_dir):
    source = Path(data_dir) / record.raw / "outs"
    target = Path(out_dir) / "01-CellRanger" / record.sample
    target.mkdir(parents=True, exist_ok=True)
    for name in ("metrics_summary.csv", "cloupe.cloupe", "web_summary.html"):
        path = source / name
        if path.exists():
            shutil.copy2(path, target / name)
    matrix = source / "filtered_feature_bc_matrix"
    if matrix.exists():
        shutil.copytree(matrix, target / matrix.name, dirs_exist_ok=True)


# 按样本表读入矩阵，以共同基因合并样本，并可复制 Cell Ranger 结果。
def load_10x_data(
    data_dir,
    sample_info,
    out_dir=None,
    copy_cellranger_outputs=False,
    cache=True,
):
    import anndata

    info = (
        read_sample_info(sample_info)
        if isinstance(sample_info, (str, Path))
        else sample_info.copy()
    )
    records = list(info.itertuples(index=False))
    adatas = [_read_10x_record(record, data_dir, cache=cache) for record in records]
    adata = anndata.concat(adatas, join="inner", merge="same", index_unique=None)
    adata.var_names_make_unique()

    if copy_cellranger_outputs:
        for record in records:
            _copy_cellranger_outputs(record, data_dir, out_dir)
    return adata


# 将 group-to-samples 字典转换为样本表，再复用多样本读取函数。
def get_10x_mtx(data_path, sample_dict, out_dir=None):
    rows = [
        {"raw": sample, "sample": sample, "group": group}
        for group, samples in sample_dict.items()
        for sample in samples
    ]
    return load_10x_data(
        data_dir=data_path,
        sample_info=pd.DataFrame(rows),
        out_dir=out_dir,
        copy_cellranger_outputs=out_dir is not None,
    )
