# 导出旧版 AnnData 兼容副本；不覆盖已执行步骤或原注释对象。
from pathlib import Path
from scrna_workflow.compatibility import convert_h5ad_legacy
from scrna_workflow.project import load_project

project = load_project(Path(__file__).resolve().parents[1])
convert_h5ad_legacy(project.ADATA_ANNOTATED, project.PROCESSED_DIR / "05_adata_annotated_legacy.h5ad",
                    legacy_python=getattr(project, "LEGACY_PYTHON", None))
