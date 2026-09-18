# 计算过滤前 QC，输出全样本和逐样本指标图表。

from pathlib import Path

from scrna_workflow.project import load_project
from scrna_workflow.steps import step_qc_before_filter

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_qc_before_filter(load_project(Path(__file__).resolve().parents[1]))
