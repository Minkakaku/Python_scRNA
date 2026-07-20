# 新项目分析顺序

先在项目根目录执行：

```bash
conda activate python-scrna-workflow
pip install -e .
```

每个分析点新建一个 `.py` 文件，已运行的文件不覆盖。文件名使用下划线，并在项目根目录用 `python -m workflow.stepXX_name` 运行。

## Step 01：读取 Cell Ranger

`workflow/step01_load.py`：

```python
from config import project
from scrna_workflow.steps import step_load

step_load(project)
```

## Step 02：过滤前 QC

`workflow/step02_qc_before_filter.py`：

```python
from config import project
from scrna_workflow.steps import step_qc_before_filter

step_qc_before_filter(project)
```

检查 `results/02-Data_QC/` 后修改 `config/project.py` 中的 QC 阈值。

## Step 03：Scrublet 和过滤

`workflow/step03_qc_filter.py`：

```python
from config import project
from scrna_workflow.steps import step_qc_filter

step_qc_filter(project)
```

Scrublet 按样本分别运行，避免不同测序深度的样本共用同一分布。

## Step 04：全局聚类

`workflow/step04_global_clustering.py`：

```python
from config import project
from scrna_workflow.steps import step_global_clustering

step_global_clustering(project)
```

`BATCH_EFFECT` 可设为 `None`、`"harmony"` 或 `"bbknn"`。改变 PC、邻居数或 resolution 时新建下一份 workflow，不覆盖已运行脚本。

## Step 05：marker

`workflow/step05_markers.py`：

```python
from config import project
from scrna_workflow.steps import step_find_markers

step_find_markers(project)
```

根据 marker 结果填写 `config/celltype_annotation.tsv`。

## Step 06：细胞注释

`workflow/step06_annotation.py`：

```python
from config import project
from scrna_workflow.steps import step_cell_annotation

step_cell_annotation(project)
```

默认使用人工 cluster-to-cell-type 表。PanglaoDB 和 CellTypist 仅作为辅助；CellTypist 模型必须在配置中明确指定。

## Step 07：细胞比例

`workflow/step07_cell_proportion.py`：

```python
from config import project
from scrna_workflow.steps import step_cell_proportion

step_cell_proportion(project)
```

组间检验以样本细胞比例为重复，输出样本级表、Mann–Whitney 检验表、PNG 和 PDF。

## 按需：基因表达

```python
from config import project
from scrna_workflow.steps import step_visualize_genes

step_visualize_genes(project, ["CD3D", "CD8A", "MS4A1"])
```

Feature UMAP 使用统一坐标和用户惯用的浅色到深蓝渐变。
