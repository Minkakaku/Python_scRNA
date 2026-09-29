from pathlib import Path
from scrna_workflow.project import load_project
from scrna_workflow.steps import step_find_markers

project = load_project(Path(__file__).resolve().parents[1])
step_find_markers(project)
