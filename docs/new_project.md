# 新项目使用指南

## 1. 建立项目

复制 `projects/project1/`，修改项目名和配置，不复制复用源码：

```bash
cp -r projects/project1 projects/my_project
```

`config/samples.tsv` 的 `raw` 列对应项目 `data/raw/<raw>/outs/` 下的 Cell Ranger 输出。

## 2. 配置覆盖

读取顺序为根目录 `config/default.yaml`、`config/reference.yaml`、项目 `config/config.yaml`。后读取的同名顶级键覆盖前面的值。

## 3. 逐步运行

```bash
python projects/my_project/workflow/00_setup.py
python projects/my_project/workflow/01_load.py
python projects/my_project/workflow/02_qc_before_filter.py
```

检查项目 `results/02-Data_QC/` 后修改项目 `config/config.yaml` 中的 QC 阈值，再运行 `03_qc_filter.py`。Scrublet 按样本分别运行。

聚类后运行 marker，根据结果填写项目 `config/celltype_annotation.tsv`，再执行注释和样本级比例分析。PanglaoDB 和 CellTypist 只作为辅助；CellTypist 模型必须显式配置。

## 4. 按需基因表达

在项目 `config/config.yaml` 填写 `visualization_genes`，再新建下一编号 workflow 调用：

```python
from pathlib import Path
from scrna_workflow.project import load_project
from scrna_workflow.steps import step_visualize_genes

project = load_project(Path(__file__).resolve().parents[1])
step_visualize_genes(project, project.VISUALIZATION_GENES)
```

已运行脚本不要覆盖。改变 PC、邻居数、批次校正或 resolution 时新增下一编号脚本。

## 5. 多平台矩阵

设置项目 `input_platform`，仍使用 `raw/sample/group` 样本表；混合平台可增加 `platform` 列。`raw` 相对项目 `data/raw/`，也可直接指向矩阵目录。

| platform | 样本目录中的默认矩阵位置 |
|---|---|
| 10x | `outs/filtered_feature_bc_matrix/` 或 `.h5` |
| c4 | `output/filter_matrix/` |
| seek | `Outs/filtered_feature_bc_matrix/` 或 `<原目录名>/Outs/filtered_feature_bc_matrix/` |
| bd | `filtered_feature_bc_matrix/` |

不同平台或参考基因集合混用前需核对 gene symbol；默认仅取共同基因，显式 `matrix_join: outer` 才取并集并补零。

## 6. 可选 VDJ

安装 `pip install -e ".[vdj]"`，在项目 `config/vdj_samples.tsv` 写入：

```text
sample	contigs
Con_1	data/raw/vdj_Con_1/filtered_contig_annotations.csv
Exp_1	data/raw/vdj_Exp_1/filtered_contig_annotations.csv
```

路径相对项目目录；sample 与已注释 H5AD 的 `obs.sample` 一致，`obs.barcode` 必须为原完整 barcode。运行 `08_vdj.py`，对象写入 `data/processed/06_vdj.h5mu`，图表写入 `results/08-VDJ/`。所有 UMAP 复用原坐标，统一 1:1，使用全对象范围和相同面板尺寸；不重新拟合嵌入。

默认 alignment/cutoff=15 是相似氨基酸序列克隆簇；需要精确氨基酸序列时显式设 `vdj_metric: identity`。克隆大小和多样性按样本统计，不自动推断抗原特异性，不自动启动 motif 或外部 VDJDB 分析。

## 7. 旧 H5AD 与细胞类型 marker

`09_export_legacy.py` 从已注释对象输出独立旧格式副本；目标已存在时中止。也可直接运行：

```bash
python -m scrna_workflow.compatibility annotated.h5ad annotated_legacy.h5ad --legacy-python /path/to/old-env/bin/python
```

`10_celltype_markers.py` 输出 `results/10-Celltype_markers/all_cell_type_marker.csv`，含 `Symbol` 和真实 `cell_type`。少于两个细胞的类型不参与检验，正向 marker 用于细胞类型描述与后续富集，不替代 CON/EXP 的样本层面差异分析。

每个 `celltype_##/all_tested_genes.csv` 保留该类型对其余类型的全部可检验基因、正负 logFC、P 值与统计量，供后续富集。复用已有 R 富集方法时映射 `--gene-col Symbol --logfc-col logfoldchanges --pvalue-col pvals --padj-col pvals_adj --stat-col scores`，按单个细胞类型处理；原富集模板的 CON/EXP 解释文案需按实际类型-vs-rest含义调整。筛选后的 marker 表不能用作全基因 GSEA 排序或 ORA 的全部背景。
