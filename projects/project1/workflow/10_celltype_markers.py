# 在已注释对象上计算细胞类型 marker，复用共享统计与绘图函数。
from pathlib import Path
from scrna_workflow.project import load_project
from scrna_workflow.steps import step_celltype_markers

step_celltype_markers(load_project(Path(__file__).resolve().parents[1]))
