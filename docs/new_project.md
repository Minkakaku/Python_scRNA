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
