# 根据项目样本表读取并合并表达矩阵，保存原始 H5AD。

from pathlib import Path

from scrna_workflow.project import load_project
from scrna_workflow.steps import step_load

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_load(load_project(Path(__file__).resolve().parents[1]))
