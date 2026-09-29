from pathlib import Path
from scrna_workflow.project import load_project
from scrna_workflow.steps import step_load

project = load_project(Path(__file__).resolve().parents[1])
step_load(project)
