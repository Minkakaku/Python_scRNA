# 按项目阈值过滤细胞并逐样本识别双细胞，保存过滤后对象。

from pathlib import Path

from scrna_workflow.project import load_project
from scrna_workflow.steps import step_qc_filter

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_qc_filter(load_project(Path(__file__).resolve().parents[1]))
