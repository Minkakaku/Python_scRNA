# project1

这是可复制的 Scanpy 单细胞项目入口。样本、QC、聚类参数和人工注释放在 `config/`；短入口放在 `workflow/`；中间 H5AD 写入 `data/processed/`；最终图表写入 `results/`。

在仓库根目录安装公共环境和包：

```bash
conda env create -f envs/environment.yaml
conda activate python-scrna-workflow
pip install -e .
```

修改 `config/samples.tsv` 和 `config/config.yaml` 后逐步运行：

```bash
python projects/project1/workflow/00_setup.py
python projects/project1/workflow/01_load.py
python projects/project1/workflow/02_qc_before_filter.py
```

检查 `results/02-Data_QC/` 后再运行过滤；marker 完成后填写 `config/celltype_annotation.tsv`。已运行的 workflow 不覆盖，改参或重跑时新增下一编号脚本。
