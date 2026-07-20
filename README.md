# Python scRNA-seq 分析模板

这是从原始 `scanpy_pipeline(1).py` 整理出的可复用 Scanpy 单细胞分析模板。原始脚本逐字节保存在 `legacy/original/`；正式代码只使用 Python 脚本，不使用 notebook 或 Streamlit。

## 目录约定

- `src/`：可跨项目复用的函数和分析方法。
- `config/`：当前项目的样本表、QC、聚类和人工注释。
- `workflow/`：项目开始时为空；每运行一个分析点，新建 `stepXX_name.py`。
- `legacy/original/`：原始脚本，只追溯、不修改。
- `data/`、`results/`、`logs/`：输入、中间对象和结果，不提交 Git。

已运行的 workflow 不覆盖。参数变化、重聚类或补充分析使用新的步骤文件。

## 环境

```bash
conda env create -f environment.yml
conda activate python-scrna-workflow
pip install -e .
```

## 新项目开始

1. 复制 `config/sample_info.example.tsv` 为 `config/sample_info.tsv`。
2. 修改 `config/project.py`。
3. 把 Cell Ranger 输出放到 `data/<raw>/outs/`。
4. 在 `workflow/` 新建本次步骤脚本。
5. 在项目根目录用模块方式运行，例如 `python -m workflow.step01_load`。

第一个步骤 `workflow/step01_load.py`：

```python
from config import project
from scrna_workflow.steps import step_load

step_load(project)
```

标准顺序：

1. 读取 Cell Ranger；
2. 过滤前 QC；
3. QC、Scrublet 和过滤；
4. 标准化、批次校正、聚类和 UMAP；
5. marker；
6. 人工或辅助注释；
7. 样本级细胞比例比较；
8. 按需进行基因可视化。

完整代码见 [新项目指南](docs/new_project.md)，原代码修正记录见 [审查记录](docs/original_code_review_CN.md)。

## 统计与绘图

- CON/EXP 比较先汇总到样本，不把每个细胞当作独立重复。
- CON 使用 `#72b1b9`，EXP 使用 `#f6ae00`。
- split UMAP 使用相同坐标范围和固定纵横比。
- 最终图保存 PNG 和 PDF，结果与绘图输入表保存 CSV。
- 人工注释使用生物学名称，主图不展示技术性 cluster 名作为最终结论。

真实数据的 QC 阈值、批次校正方法、resolution 和细胞注释仍需按项目决定。
