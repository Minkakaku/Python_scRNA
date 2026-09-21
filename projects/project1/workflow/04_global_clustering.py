# 按项目参数完成归一化、降维、批次处理与全局聚类。

from pathlib import Path

from scrna_workflow.project import load_project
from scrna_workflow.steps import step_global_clustering

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_global_clustering(load_project(Path(__file__).resolve().parents[1]))
