#!/usr/bin/env python3
"""导出旧版 AnnData 兼容 H5AD；原输入保持不变，技术编码差异仅打印到终端。"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd


DROP = object()


# 解析显式输入输出路径与可选的旧环境验证参数。
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="h5ad-convert-legacy",
        description="Rewrite a modern H5AD while removing encodings unsupported by older AnnData.",
    )
    p.add_argument("input", type=Path, help="Input .h5ad file")
    p.add_argument("output", type=Path, help="Output legacy-compatible .h5ad file")
    p.add_argument(
        "--compression",
        choices=("gzip", "lzf", "none"),
        default="gzip",
        help="Output compression (default: gzip)",
    )
    p.add_argument("--overwrite", action="store_true", help="Replace an existing output file")
    p.add_argument(
        "--legacy-python",
        type=Path,
        help="Optional old-environment Python executable used to validate the output",
    )
    return p.parse_args()


# 递归清理旧版 AnnData 不支持的 null 编码，不改变表达矩阵。
def clean_value(value: Any, path: str, changes: list[str]) -> Any:
    if value is None:
        changes.append(f"DROP null mapping value: {path}")
        return DROP

    if isinstance(value, Mapping):
        cleaned = {}
        for key, child in value.items():
            result = clean_value(child, f"{path}/{key}", changes)
            if result is not DROP:
                cleaned[key] = result
        return cleaned

    if isinstance(value, list):
        out = []
        for i, child in enumerate(value):
            result = clean_value(child, f"{path}[{i}]", changes)
            if result is DROP:
                changes.append(f"REPLACE null list item with empty string: {path}[{i}]")
                result = ""
            out.append(result)
        return out

    if isinstance(value, tuple):
        return tuple(clean_value(list(value), path, changes))

    if isinstance(value, np.ndarray) and value.dtype == object:
        out = value.copy()
        mask = np.frompyfunc(lambda x: x is None, 1, 1)(out).astype(bool)
        if mask.any():
            out[mask] = ""
            changes.append(f"REPLACE {int(mask.sum())} null object-array item(s): {path}")
        return out

    return value


# 把 pandas 可空扩展类型转换为旧格式可读的类型；缺失布尔值转为 NaN。
def normalize_dataframe(df: pd.DataFrame, label: str, changes: list[str]) -> None:
    for column in df.columns:
        series = df[column]
        dtype = series.dtype
        col_path = f"{label}/{column}"

        if isinstance(dtype, pd.StringDtype):
            df[column] = series.astype(object).where(series.notna(), "")
            changes.append(f"CONVERT pandas string dtype to object: {col_path}")
        elif pd.api.types.is_extension_array_dtype(dtype) and not isinstance(
            dtype, pd.CategoricalDtype
        ):
            if pd.api.types.is_integer_dtype(dtype):
                if series.isna().any():
                    df[column] = series.astype("float64")
                else:
                    df[column] = series.astype(dtype.numpy_dtype)
            elif pd.api.types.is_bool_dtype(dtype):
                if series.isna().any():
                    df[column] = series.map({True: 1.0, False: 0.0}).astype("float64")
                else:
                    df[column] = series.astype(bool)
            else:
                df[column] = series.astype(object).where(series.notna(), "")
            changes.append(f"CONVERT extension dtype {dtype}: {col_path}")


# 统一 obs、var、raw 元数据和 uns 的兼容编码。
def normalize_adata(adata: ad.AnnData) -> list[str]:
    changes: list[str] = []
    normalize_dataframe(adata.obs, "obs", changes)
    normalize_dataframe(adata.var, "var", changes)

    if adata.raw is not None:
        normalize_dataframe(adata.raw.var, "raw/var", changes)

    for name, value in list(adata.obsm.items()):
        if isinstance(value, pd.DataFrame):
            normalize_dataframe(value, f"obsm/{name}", changes)

    for name, value in list(adata.varm.items()):
        if isinstance(value, pd.DataFrame):
            normalize_dataframe(value, f"varm/{name}", changes)

    cleaned_uns = clean_value(dict(adata.uns), "uns", changes)
    adata.uns.clear()
    adata.uns.update(cleaned_uns)
    adata.strings_to_categoricals()
    return changes


# 删除写出后残留的 null 节点，保证旧读者不遇到未知编码。
def remove_null_encoded_nodes(path: Path, changes: list[str]) -> None:
    to_delete: list[str] = []
    with h5py.File(path, "r") as handle:
        # 只记录 null 节点路径，关闭读句柄后再删除。
        def visitor(name: str, obj: h5py.Group | h5py.Dataset) -> None:
            encoding = obj.attrs.get("encoding-type", obj.attrs.get("encoding_type", ""))
            if isinstance(encoding, bytes):
                encoding = encoding.decode("utf-8", "replace")
            if encoding == "null":
                to_delete.append("/" + name)
        handle.visititems(visitor)

    if not to_delete:
        return

    with h5py.File(path, "r+") as handle:
        for key in sorted(to_delete, key=lambda x: x.count("/"), reverse=True):
            if key in handle:
                del handle[key]
                changes.append(f"DROP residual null HDF5 node: {key}")


# 以只读 backed 模式回读目标，检查矩阵维度。
def validate_with_current(path: Path, expected_shape: tuple[int, int]) -> None:
    check = ad.read_h5ad(path, backed="r")
    try:
        if check.shape != expected_shape:
            raise RuntimeError(f"shape changed: expected {expected_shape}, got {check.shape}")
    finally:
        check.file.close()


# 可选地用指定旧 Python 回读并检查同样的矩阵维度。
def validate_with_legacy(path: Path, python: Path, expected_shape: tuple[int, int]) -> None:
    code = (
        "import anndata as ad,sys;"
        "a=ad.read_h5ad(sys.argv[1],backed='r');"
        "assert a.shape==tuple(map(int,sys.argv[2:])), 'legacy shape changed';"
        "print('legacy_read_ok shape='+str(a.shape));"
        "a.file.close()"
    )
    result = subprocess.run(
        [str(python), "-c", code, str(path), *map(str, expected_shape)],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode != 0:
        raise RuntimeError("legacy validation failed:\n" + result.stdout)
    print(result.stdout.strip())


# 写入同目录临时文件并回读检查，成功后原子替换目标，禁止覆盖输入本身。
def convert_h5ad_legacy(input_path, output_path, compression="gzip", overwrite=False, legacy_python=None):
    src = Path(input_path).expanduser().resolve()
    dst = Path(output_path).expanduser().resolve()
    if not src.is_file():
        raise FileNotFoundError(src)
    if src == dst or (dst.exists() and os.path.samefile(src, dst)):
        raise ValueError("Input and output must be different files")
    if dst.exists() and not overwrite:
        raise FileExistsError(dst)
    if legacy_python is not None and not Path(legacy_python).is_file():
        raise FileNotFoundError(legacy_python)
    if compression not in {"gzip", "lzf", "none", None}:
        raise ValueError("compression must be gzip, lzf or none")
    dst.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{dst.name}.", suffix=".h5ad", dir=dst.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        adata = ad.read_h5ad(src)
        expected_shape = adata.shape
        changes = normalize_adata(adata)
        adata.write_h5ad(tmp, compression=None if compression == "none" else compression)
        remove_null_encoded_nodes(tmp, changes)
        validate_with_current(tmp, expected_shape)
        if legacy_python is not None:
            validate_with_legacy(tmp, Path(legacy_python).resolve(), expected_shape)
        # 不先删除已有目标，避免验证或替换失败时丢失旧输出。
        os.replace(tmp, dst)
    finally:
        if tmp.exists():
            tmp.unlink()
    for item in changes:
        print(item)
    print(f"conversion_ok: {dst}")
    return dst


# CLI 与项目函数调用共用同一实现；验证失败明确返回非零退出码。
def main():
    args = parse_args()
    try:
        convert_h5ad_legacy(args.input, args.output, args.compression, args.overwrite, args.legacy_python)
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
