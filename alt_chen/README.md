# 定量蛋白质组学标准分析 Pipeline

从蛋白组定量矩阵（DIA-NN TSV 或 xlsx 基因表）出发，一键完成「标准化 → QC 聚类 → 差异分析 → 富集 → GSEA → GSVA → 通路热图 → 网页报告」的全流程标准分析。支持**人类 / 小鼠 / 大鼠**三种物种，组别名与组数可自由填写。

## 运行环境（固定）

- R 4.4.3：`/media/nfs/nfs02/chenyuanbao/miniforge3/envs/R433/bin/Rscript`
- Python 3.11：`/media/nfs/nfs02/chenyuanbao/miniforge3/envs/py311/bin/python`
- 依赖 R 包：limma / fgsea / GSVA / clusterProfiler / ComplexHeatmap / ggplot2 / ggrepel / GSEABase / ggsci / Rtsne / dendextend / pheatmap / ragg / showtext / msigdbr / org.Mm.eg.db / org.Hs.eg.db / KEGGREST 等（均已安装在 R433）
- ImageMagick：`/usr/bin/convert`（QC 图 pdf → png）

## 目录结构

```
Quantitative_Proteomics_Standard_Analysis_test/
├── config.sh               # 人类数据默认配置（全部参数）
├── config.ima94.sh         # 小鼠数据示例配置（可通过 CONFIG= 切换）
├── run_pipeline.sh         # 总控脚本（支持 CONFIG= 与 step0 转换）
├── bin/                    # 参数化分析脚本
│   ├── 00_convert_ima94.py        # 可选: xlsx 基因表 -> 标准 TSV
│   ├── 00_build_gmt.R            # 可选: 生成任意物种基因集 GMT (human/mouse/rat)
│   ├── 01_process_standard.R       # 标准化（过滤/log2/填充 + 逐样本缺失率）
│   ├── 02_cluster_pca.R            # PCA/t-SNE/层次聚类 QC
│   ├── 03_DE_limma.R               # limma 差异分析
│   ├── 03b_volcano.R               # 差异火山图
│   ├── 04_ORA_enrichment.R         # 三库 ORA 富集 + 气泡图
│   ├── 05_full_GSEA.R              # 全量 GSEA
│   ├── 06_GSEA_plots.R             # GSEA 激活/抑制柱状图
│   ├── 07a_prepare_GSVA_input.R    # GSVA 输入准备
│   ├── 07b_GSVA_calc.R             # GSVA 打分
│   ├── 07c_GSVA_top_heatmap.R      # GSVA 两方向热图
│   ├── 08a_pathway_protein_heatmaps.R  # ORA 通路蛋白热图
│   ├── 08b_GSEA_leading_heatmaps.R     # GSEA 前导热图
│   └── make_report.py              # 网页报告生成
├── lib/                    # 依赖文件（打包自包含）
│   ├── databases/GMT/              # 人类 KEGG/GOBP/Reactome GMT
│   ├── databases/GMT_mm/           # 小鼠 KEGG/GOBP/Reactome GMT
│   ├── databases/GMT_rn/           # 大鼠 KEGG/GOBP/Reactome GMT
│   ├── fonts/Arial/                # Arial 字体
│   ├── heatmap/heatmap_light.r     # 热图绘图辅助（--title 标题支持）
│   └── id_map.tsv                  # 通路名→短ID 映射
├── input/                  # 输入数据
├── output/                 # 输出（按分析步骤组织）
└── logs/
```

## 快速使用（单用户 / 直接改 pipeline 内 config）

```bash
cd /media/nfs/nfs02/chenyuanbao/pipeline

# 1. 编辑 config.sh，设置物种(PP_SPECIES)/组别(PP_GROUPS)/输入(PP_MATRIX)/输出(PP_OUTPUT)

# 2. 放置输入数据（默认 ${HOME}/proteomics_analysis/input/ 或自定义 PP_INPUT）

# 3. 运行（CONFIG 指定配置，默认 config.sh；可分步或全跑）
bash run_pipeline.sh all                    # 全流程（human, config.sh）
CONFIG=config.ima94.sh bash run_pipeline.sh all   # mouse 示例配置
bash run_pipeline.sh step1                  # 仅标准化
```

## 多用户共享使用（推荐）

pipeline 为**共享只读**：所有用户共用 `bin/`、`lib/`（三物种基因集），各自的输入/输出/日志落在**自己的目录**，互不干扰。

1. 复制用户模板到自己的目录、重命名并填好 `<...>` 占位项：
   ```bash
   cp /media/nfs/nfs02/chenyuanbao/pipeline/config.user.template.sh ~/my_run.sh
   vi ~/my_run.sh    # 填 PP_SPECIES / PP_GROUPS / PP_MATRIX / PP_OUTPUT 等
   ```
2. 从任意目录以绝对路径指定配置文件运行：
   ```bash
   CONFIG=~/my_run.sh bash /media/nfs/nfs02/chenyuanbao/pipeline/run_pipeline.sh all
   # 或分步 step0..step8; 输出默认在 ~/proteomics_analysis/output/（可用 PP_OUTPUT 改）
   ```
3. 也可用环境变量临时覆盖输入/输出/工作目录，不修改配置文件：
   ```bash
   PP_WORKDIR=/data/你的目录 PP_OUTPUT=/data/结果 CONFIG=~/my_run.sh      bash /media/nfs/nfs02/chenyuanbao/pipeline/run_pipeline.sh all
   ```

## 输入数据格式（两种）

**1) DIA-NN 蛋白矩阵 TSV**（人类数据默认）：前 6 列为注释（`Protein.Group / Protein.Names / Genes / First.Protein.Description / N.Sequences / N.Proteotypic.Sequences`），其后为各样本定量列，**列顺序 = 样本映射顺序**（`部位 × 组别 × 重复`）。

**2) xlsx 基因总定量表**（如小鼠 `IMA94_蛋白总定量表.xlsx`）：列含 `gene_id` + 样本定量列（+ 可选的 Fold.Change/P.Value，会被忽略）。设 `PP_MATRIX_XLSX` 指向它，step0 会用 `00_convert_ima94.py` 自动转成标准 TSV（把 0 值视为缺失、按 `PP_GROUPS` 顺序重排列、补齐 6 注释列）。

## 实验设计（config.sh 关键变量）

| 变量 | 含义 | 示例 |
|---|---|---|
| PP_SPECIES | 物种（一键切换） | `human` / `mouse` / `rat` |
| PP_REGIONS | 部位/批次列表（空格分隔） | `IVS` 或 `IVS LVAW LVPW` |
| PP_REGION_TAG | 单一"部位"标签（xlsx 转换用） | `IMA` |
| PP_GROUPS | 组别列表（**任意命名/任意组数**，对照在前） | `Normal Model` / `CMC HY0094` / `WT KO` |
| PP_REPLICATES | 每组重复数 | `6` 或 `3` |
| PP_CONTRASTS | 对比（`A_vs_B` = A − B） | `Model_vs_Normal` |
| PP_GROUP_CONTROL | GSVA 对照组 | `Normal` 或 `CMC` |
| PP_GROUP_TREAT | GSVA 实验组 | `Model` 或 `HY0094` |
| PP_MATRIX | 输入 TSV 路径 | `input/report.pg_matrix.tsv` |
| PP_MATRIX_XLSX | 输入 xlsx 路径（可选，触发 step0 转换） | `input/IMA94_蛋白总定量表.xlsx` |
| PP_PCUT | pvalue 阈值 | `0.05` |
| PP_FCCUT | \|log2FC\| 阈值（log2 尺度） | `0.2630344058337938`（=log2(1.2)） |
| PP_COLOR_\<GROUP\> | 各组颜色（按组名动态） | `PP_COLOR_CMC` |

## 物种选择（一键切换）

配置里设 `PP_SPECIES=human|mouse|rat`，脚本自动选定对应基因集（org.db + KEGG/GOBP/Reactome 三库 GMT），无需手填路径：

| PP_SPECIES | 基因集目录 | 符号风格 | 说明 |
|---|---|---|---|
| `human` | `lib/databases/GMT/` | HUGO 全大写 | KEGG hsa + C5 GO:BP + C2 REACTOME（原生） |
| `mouse` | `lib/databases/GMT_mm/` | MGI 首字母大写 | KEGG mmu + M5 GO:BP + M2 REACTOME（原生） |
| `rat` | `lib/databases/GMT_rn/` | RGNC 首字母大写 | KEGG rno（原生）+ GO/REACTOME（人类库同源映射） |

三套 GMT 均已生成；如需重建，`PP_SPECIES=<x> PP_GMT_OUTPUT=<dir> Rscript bin/00_build_gmt.R`（需网络，大鼠 GO/Reactome 用 msigdbr 同源映射）。

## 输出目录（按分析步骤）

```
output/
├── 01_raw/                # 原始矩阵副本
├── 02_standardize/        # 标准化（log2 + 填充矩阵 + 样本信息 + QA）
├── 03_qc_pca_cluster/     # PCA / t-SNE / 层次聚类图
├── 04_de/                 # 差异分析（DE 表 + 火山图）
├── 05_ora/                # ORA 富集（enrichment 表 + 气泡图）
├── 06_gsea/               # GSEA（fgsea 结果 + 激活/抑制柱状图）
├── 07_gsva/               # GSVA（input + score + 两方向热图）
├── 08_pathway_heatmaps/   # 通路蛋白热图 + GSEA 前导热图
└── 09_report/             # html_report.html 网页报告
```

通过配置里 `PP_OUTPUT` 可指定独立输出子目录（如 `output/IMA94/`），互不影响。

## 报告

`output/09_report/html_report.html` 双击打开即可，图片直接引用各步骤目录（无冗余副本）。报告 2.2 为「各样本缺失率」。

## 备注

- 所有脚本已参数化（读配置导出的环境变量），新增数据集只需改配置，不改脚本。
- GSVA 方向 = `PP_GROUP_TREAT 均值 − PP_GROUP_CONTROL 均值`（即实验组比对照组）。
- 分组配色按 `PP_COLOR_<GROUP>` 动态读取（未配置则回退调色板），任意组数/组名兼容。
- QC 质控表、缺失率表、富集库路径、对比名等均动态解析，无硬编码部位/组名。
