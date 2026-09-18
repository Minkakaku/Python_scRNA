# 读取注释对象，汇总细胞组成并在样本层面进行组间比较。

from pathlib import Path

from scrna_workflow.project import load_project
from scrna_workflow.steps import step_cell_proportion

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_cell_proportion(load_project(Path(__file__).resolve().parents[1]))
