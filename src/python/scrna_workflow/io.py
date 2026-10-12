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


# 统一平台名称；仅识别已有公共流程支持的四种输入格式。
def _platform_key(platform):
    aliases = {"10x": "10x", "cellranger": "10x", "c4": "c4", "华大": "c4",
               "seek": "seek", "xunyin": "seek", "寻因": "seek", "bd": "bd"}
    key = str(platform).strip().lower()
    if key not in aliases:
        raise ValueError(f"Unsupported input platform: {platform}")
    return aliases[key]


# 按平台定位矩阵；也允许样本表 raw 直接指向一个矩阵目录。
def _matrix_path(raw_dir, platform="10x"):
    raw_dir = Path(raw_dir)
    key = _platform_key(platform)
    candidates = {
        "10x": [raw_dir / "outs" / "filtered_feature_bc_matrix",
                raw_dir / "outs" / "filtered_feature_bc_matrix.h5"],
        "c4": [raw_dir / "output" / "filter_matrix", raw_dir / "filter_matrix"],
        "seek": [raw_dir / "Outs" / "filtered_feature_bc_matrix",
                 raw_dir / raw_dir.name / "Outs" / "filtered_feature_bc_matrix"],
        "bd": [raw_dir / "filtered_feature_bc_matrix"],
    }[key]
    candidates += [raw_dir / "filtered_feature_bc_matrix", raw_dir]
    for path in candidates:
        if path.is_file() and path.suffix == ".h5":
            return path
        if path.is_dir() and any((path / n).is_file() for n in ("matrix.mtx", "matrix.mtx.gz")):
            return path
    raise FileNotFoundError(f"No {key} expression matrix under {raw_dir}")


# C4/BD 的第一列是基因名；不能按 10x 的第二列规则误读。
def _read_platform_matrix(path, platform, cache=True):
    key = _platform_key(platform)
    if key in {"10x", "seek"}:
        import scanpy as sc
        return (sc.read_10x_mtx(path, var_names="gene_symbols", cache=cache)
                if path.is_dir() else sc.read_10x_h5(path))
    import anndata as ad
    from scipy.io import mmread
    from scipy.sparse import csr_matrix

    # 同时兼容压缩和未压缩的 MEX 文件。
    def component(stem):
        for name in (stem, stem + ".gz"):
            if (path / name).is_file():
                return path / name
        raise FileNotFoundError(f"Missing {stem} under {path}")

    genes = pd.read_csv(component("features.tsv"), sep="\t", header=None, dtype=str)
    barcodes = pd.read_csv(component("barcodes.tsv"), sep="\t", header=None, dtype=str)[0]
    if key == "bd":
        # 已有 CellLabel_ 前缀原样保留，不截短或重复加前缀。
        barcodes = barcodes.map(lambda x: x if x.startswith("CellLabel_") else "CellLabel_" + x)
    matrix = csr_matrix(mmread(component("matrix.mtx")).T, dtype="float32")
    if matrix.shape != (len(barcodes), len(genes)):
        raise ValueError("Matrix shape does not match barcode and feature tables")
    return ad.AnnData(matrix, obs=pd.DataFrame(index=pd.Index(barcodes.values)),
                      var=pd.DataFrame({"gene_symbols": genes[0].values}, index=pd.Index(genes[0].values)))


# 保留真实 barcode，再用 sample:barcode 建立跨平台、跨样本唯一索引。
def _read_10x_record(record, data_dir, cache=True, platform="10x"):
    key = _platform_key(getattr(record, "platform", "") or platform)
    matrix = _matrix_path(Path(data_dir) / record.raw, key)
    adata = _read_platform_matrix(matrix, key, cache=cache)
    adata.var_names_make_unique()
    if not adata.obs_names.is_unique:
        raise ValueError(f"Duplicate barcodes in sample {record.sample}")
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs["sample"] = record.sample
    adata.obs["group"] = record.group
    adata.obs["platform"] = key
    adata.obs_names = pd.Index(record.sample + ":" + adata.obs["barcode"].astype(str), name="cell_id")
    return adata


# 按原平台复制矩阵与已有报告；只有项目显式启用时才复制这些上游产物。
def _copy_cellranger_outputs(record, data_dir, out_dir, platform="10x"):
    key = _platform_key(getattr(record, "platform", "") or platform)
    matrix = _matrix_path(Path(data_dir) / record.raw, key)
    source = matrix.parent
    folder = {"10x": "01-CellRanger", "c4": "01-Dnbc4tools",
              "seek": "01-seeksoultools", "bd": "01-Rhapsody"}[key]
    target = Path(out_dir) / folder / record.sample
    target.mkdir(parents=True, exist_ok=True)
    for name in ("metrics_summary.csv", "metrics_summary.xls", "cloupe.cloupe", "web_summary.html",
                 f"{record.raw}_summary.csv", f"{record.raw}_report.html",
                 f"{record.raw}_Metrics_Summary.csv", f"{record.raw}_scRNA_report.html"):
        path = source / name
        if path.is_file():
            shutil.copy2(path, target / path.name)
    if matrix.is_dir():
        shutil.copytree(matrix, target / matrix.name, dirs_exist_ok=True)
    else:
        shutil.copy2(matrix, target / matrix.name)


# 从统一样本表读取四种平台，可逐行指定 platform；outer 显式补零，默认仍为 inner。
def load_expression_data(data_dir, sample_info, out_dir=None,
                         copy_cellranger_outputs=False, cache=True, platform="10x", join="inner"):
    import anndata

    info = read_sample_info(sample_info) if isinstance(sample_info, (str, Path)) else sample_info.copy()
    if info.empty or not {"raw", "sample", "group"}.issubset(info.columns):
        raise ValueError("Sample table requires non-empty raw, sample and group columns")
    if info[["raw", "sample", "group"]].isna().any().any() or info[["raw", "sample", "group"]].eq("").any().any():
        raise ValueError("Sample table contains empty identifiers")
    if info["sample"].duplicated().any():
        raise ValueError("Sample identifiers must be unique")
    if join not in {"inner", "outer"}:
        raise ValueError("matrix join must be inner or outer")
    if copy_cellranger_outputs and out_dir is None:
        raise ValueError("out_dir is required when copying upstream outputs")
    records = list(info.itertuples(index=False))
    adatas = [_read_10x_record(record, data_dir, cache, platform) for record in records]
    adata = anndata.concat(adatas, join=join, merge="same", index_unique=None, fill_value=0)
    if adata.n_vars == 0:
        raise ValueError("Samples have no shared genes; check the input format or explicitly select outer join")
    adata.obs["sample"] = pd.Categorical(adata.obs["sample"], categories=info["sample"], ordered=True)
    adata.obs["group"] = pd.Categorical(adata.obs["group"], categories=info["group"].drop_duplicates(), ordered=True)
    if copy_cellranger_outputs:
        for record in records:
            _copy_cellranger_outputs(record, data_dir, out_dir, platform)
    return adata


# 保留原 10x 公共接口，已有项目和 get_10x_mtx 调用无需改写。
def load_10x_data(data_dir, sample_info, out_dir=None, copy_cellranger_outputs=False, cache=True):
    return load_expression_data(data_dir, sample_info, out_dir, copy_cellranger_outputs, cache)


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
