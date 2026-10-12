# 按项目 VDJ 样本表运行可选的受体链和克隆簇分析。
from pathlib import Path
from scrna_workflow.project import load_project
from scrna_workflow.steps import step_vdj

step_vdj(load_project(Path(__file__).resolve().parents[1]))
