from scrna_workflow.project import load_project


# 构造临时项目，验证项目参数覆盖公共默认，同时保留其他默认值和路径。
def test_project_config_overrides_public_defaults(tmp_path):
    workspace = tmp_path / "workspace"
    project_dir = workspace / "projects" / "project1"
    (workspace / "config").mkdir(parents=True)
    (project_dir / "config").mkdir(parents=True)
    (workspace / "config" / "default.yaml").write_text("n_pcs: 30\nrandom_seed: 1234\n")
    (workspace / "config" / "reference.yaml").write_text("organism: human\n")
    (project_dir / "config" / "config.yaml").write_text("n_pcs: 25\n")

    project = load_project(project_dir)

    assert project.N_PCS == 25
    assert project.RANDOM_SEED == 1234
    assert project.ORGANISM == "human"
    assert project.SAMPLE_INFO == project_dir / "config" / "samples.tsv"
