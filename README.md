# 定量蛋白质组下游分析流程

服务器部署代码的分类归档，快照日期 **2026-11-20**。已搜库的 PD-TMT / DIA-NN / 通用定量矩阵 → 质控 → limma 差异 → 富集/PPI/GSEA/GSVA → 中文交互报告。

## 使用入口

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r env/requirements.txt
python proteomics.py --help
python proteomics.py validate -c /绝对路径/config.yaml
python proteomics.py run -c /绝对路径/config.yaml
```

配置文件必须明确指定 Rscript、样本列、物种及比较；完整字段见 [服务器使用手册](README_流程使用说明.md)。
此流程不做原始质谱重搜库、肽段鉴定或上游定量。

### 两套流程分开使用

- **主流程**：`proteomics.py`，支持配置配对比较、FDR结果和交互报告。
- **迁移版**：[alt_chen/](alt_chen/)，bash/R入口 `run_pipeline.sh`，配置模板 `config.user.template.sh`；原始模型不处理同一动物重复采样，默认原始P值筛选，不能等同于配对分析。
- 本仓库上传服务器共享流程的当前快照，不把WLPD猫项目的专用适配脚本混进通用流程。主程序当前标记1.2.0；旧使用手册中的1.1.0为历史文档版本。
- `depends/`仅包含物种注册与配置，基因集/PPI数据库须单独准备；物种目录存在不代表每个物种都完成了真实项目验证。

## 目录

核心目录：`proteomics.py、scripts/、tools/、report_template/、env/、tests/、depends/、alt_chen/`。

- [README_流程使用说明.md](README_流程使用说明.md)：原部署使用说明。
- [部署与资源](docs/DEPLOYMENT.md)：迁移前需配置的路径、软件和外部资源。
- [验证状态](docs/VALIDATION.md)：本次实际执行的检查及未执行事项。
- [来源与署名](NOTICE.md)：来源、修改范围及第三方材料说明。
- [源文件校验清单](docs/source_manifest.json)：已归档文件与原服务器快照的SHA-256。

## 三类组学导航

- [转录组分析流程（Bulk RNA-seq）](https://github.com/black662189/bulk-rnaseq-pipeline)
- [定量蛋白质组下游分析流程](https://github.com/black662189/quantitative-proteomics-pipeline)
- [非靶向代谢组分析流程](https://github.com/black662189/untargeted-metabolomics-pipeline)

## 版本范围

本次发布为部署源码快照，不承诺任意机器零配置运行。大型数据库、参考基因组、索引、Conda/venv环境、真实样本矩阵及项目结果未上传。
服务器原目录未修改；未为上传重新运行真实项目。使用现有说明时，请将服务器绝对路径替换为自己的部署路径。
