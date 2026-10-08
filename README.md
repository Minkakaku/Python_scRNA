# Python-scrna-workflow-template

这是从原始 `scanpy_pipeline(1).py` 整理出的分层 Scanpy 单细胞主流程。公共方法与项目运行记录分开管理；原始脚本逐字节保存在 `legacy/original/`，正式入口只使用 Python 脚本。

## 项目结构

```text
├── src/python/scrna_workflow/   # 跨项目复用函数与标准步骤
├── config/                      # 公共默认与参考配置
├── projects/project1/           # 项目配置、workflow、数据、结果和日志
├── envs/                        # 公共运行环境
├── tests/
├── source_catalog/              # 其他人员代码的语义参考档案
└── legacy/original/             # 未修改的原始代码
```

`config/default.yaml` 和 `config/reference.yaml` 提供公共值，项目 `config/config.yaml` 按同名键覆盖。样本、QC、分组和人工注释均保留在具体项目中。

## 标准分析顺序

1. 读取 Cell Ranger；
2. 过滤前 QC；
3. Scrublet 与过滤；
4. 标准化、批次校正、聚类与 UMAP；
5. marker；
6. 人工或辅助注释；
7. 样本级细胞比例比较；
8. 按需进行基因可视化。

中间 H5AD 写入项目 `data/processed/`，最终图表写入项目下编号 `results/` 目录。

## 安装与运行

```bash
conda env create -f envs/environment.yaml
conda activate python-scrna-workflow
pip install -e .
python projects/project1/workflow/00_setup.py
python projects/project1/workflow/01_load.py
```

复制 `projects/project1/` 即可开始新项目。QC 后人工确认阈值再过滤；marker 后填写 `config/celltype_annotation.tsv` 再注释。已运行的 workflow 不覆盖，改参或重跑时新增下一编号脚本。

完整说明见 [新项目指南](docs/new_project.md) 和 [代码分层](docs/code_layout.md)。

## 统计与绘图

- CON/EXP 比较先汇总到样本，不把细胞当作独立重复。
- CON 使用 `#E64B35`，EXP 使用 `#4DBBD5`。
- split UMAP 使用相同坐标范围和固定纵横比。
- 正式图保存同名 PNG/PDF，结果与绘图输入表保存 CSV。
- PCA 使用 `kneed.KneeLocator()` 输出自动拐点参考，但不覆盖项目配置的 `n_pcs`。
- 人工注释使用生物学名称，主图不把技术 cluster 名作为最终结论。

`source_catalog/` 保留为参考档案，不作为正式运行入口，也不迁入项目 workflow。
