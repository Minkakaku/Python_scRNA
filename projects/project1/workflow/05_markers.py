# 读取聚类对象，计算并输出各 cluster 的 marker。

from pathlib import Path

from scrna_workflow.project import load_project
from scrna_workflow.steps import step_find_markers

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_find_markers(load_project(Path(__file__).resolve().parents[1]))
