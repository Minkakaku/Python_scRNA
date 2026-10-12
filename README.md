# Python-scrna-workflow-template

这是从原始 `scanpy_pipeline(1).py` 整理出的分层 Scanpy 单细胞主流程。公共方法与项目运行记录分开管理；原始脚本逐字节保存在 `legacy/original/`，正式入口只使用 Python 脚本。

## 项目结构

```text
├── src/python/scrna_workflow/   # 跨项目复用函数与标准步骤
├── config/                      # 公共默认与参考配置
├── projects/project1/           # 项目配置、workflow、数据、结果和日志
├── envs/                        # 公共运行环境
├── tests/
├── source_catalog/              # 其他人员代码的语义参考档案
└── legacy/original/             # 未修改的原始代码
```

`config/default.yaml` 和 `config/reference.yaml` 提供公共值，项目 `config/config.yaml` 按同名键覆盖。样本、QC、分组和人工注释均保留在具体项目中。

## 标准分析顺序

1. 读取 Cell Ranger；
2. 过滤前 QC；
3. Scrublet 与过滤；
4. 标准化、批次校正、聚类与 UMAP；
5. marker；
6. 人工或辅助注释；
7. 样本级细胞比例比较；
8. 按需进行基因可视化。

中间 H5AD 写入项目 `data/processed/`，最终图表写入项目下编号 `results/` 目录。

## 安装与运行

```bash
conda env create -f envs/environment.yaml
conda activate python-scrna-workflow
pip install -e .
python projects/project1/workflow/00_setup.py
python projects/project1/workflow/01_load.py
```

复制 `projects/project1/` 即可开始新项目。QC 后人工确认阈值再过滤；marker 后填写 `config/celltype_annotation.tsv` 再注释。已运行的 workflow 不覆盖，改参或重跑时新增下一编号脚本。

完整说明见 [新项目指南](docs/new_project.md) 和 [代码分层](docs/code_layout.md)。

## 统计与绘图

- CON/EXP 比较先汇总到样本，不把细胞当作独立重复。
- CON 使用 `#E64B35`，EXP 使用 `#4DBBD5`。
- split UMAP 使用相同坐标范围和固定纵横比。
- 正式图保存同名 PNG/PDF，结果与绘图输入表保存 CSV。
- PCA 使用 `kneed.KneeLocator()` 输出自动拐点参考，但不覆盖项目配置的 `n_pcs`。
- 人工注释使用生物学名称，主图不把技术 cluster 名作为最终结论。

`source_catalog/` 保留为参考档案，不作为正式运行入口，也不迁入项目 workflow。

## 公共服务器功能同步（2026-10-12）

已将公共 scVision5、cosmicdarahub 的矩阵读取及细胞类型 marker，以及 scRNA-tools 的 H5AD/VDJ 功能接入现有共享源码。已有 QC、Scrublet、Harmony/BBKNN、聚类和注释继续复用，未引入 Streamlit 或带缺失依赖的 Snakefile。

- `load_expression_data()`：10x、华大 C4、寻因、BD 共用样本表。`input_platform` 控制默认平台，样本表可加 `platform` 列；保留原 `load_10x_data()` 接口。默认 `matrix_join: inner`，显式 `outer` 时按基因并集补零。C4/BD 按第一列基因名读取，BD 保留 `CellLabel_`。
- `compatibility.convert_h5ad_legacy()`：导出独立旧格式副本，转换可空元数据和 null 编码；保留 X、raw、layers、坐标及标签。缺失字符串转为空串，可空整数/布尔转为兼容数值类型；可指定旧 Python 回读验证，失败保留已有目标。
- `vdj.run_vdj_analysis()` / `steps.step_vdj()`：TCR/BCR 链质控、氨基酸克隆簇、样本内克隆大小、多样性、Jaccard 共有比例和网络；精确匹配 `sample + barcode`，保留原 UMAP。相似序列克隆簇不等同于严格同一克隆，样本统计只作描述性比较。
- `steps.step_celltype_markers()`：已注释细胞类型的正向 marker，复用现有统计、NPG 与无内框线 DotPlot；输出筛选后的 `all_cell_type_marker.csv` 与各类型的 `celltype_##/all_tested_genes.csv`（全基因统计）。后续富集复用已有 `R-Python-scores-enrichment-workflow-template` 的方法，不复制服务器的硬编码 R 调用。每种细胞类型单独用全基因表，按实际类型-vs-rest解释；不将多类型 marker 混合后作为 GSEA 输入或 ORA 全背景。

新增可选入口为 `08_vdj.py`、`09_export_legacy.py`、`10_celltype_markers.py`；详细输入格式见 [新项目指南](docs/new_project.md)。VDJ 按需执行 `pip install -e ".[vdj]"`，普通流程无需安装这些额外依赖。

功能来源：公共 `/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/common.py`（2026-06-26）；`/mnt/NexData_file/cosmicdarahub/Bioinformatics/Pipeline/scRNA_pipeline/workflow/scripts/1-data_loading.py`、`8-function_enrichment.py`；`/mnt/NexData_file/share/scRNA-tools/scripts/h5ad_convert_legacy.py`（2026-08-06）及 `/mnt/NexData_file/share/scRNA-tools/VDJ/pip_VDJ.py`、`sciana.py`（VDJ 入口 2026-07-30）。此次只提取正式功能，原始档案保持不变。
