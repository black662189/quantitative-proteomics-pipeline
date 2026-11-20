# alt_chen —— 移植自 chenyuanbao 的定量蛋白质组流程

版本：v1（移植快照）。移植日期：2026-10-20。

## 定位

这是从 `/media/nfs/nfs02/chenyuanbao/pipeline` 移植而来的**独立** bash+R 风格蛋白质组流程，
作为 `proteomics/` 目录下的一个子模块（`proteomics/alt_chen/`）。

- **与主流程 `proteomics.py` 完全独立**：不改动、不共用代码或状态，单独入口（`run_pipeline.sh`）。
- 主流程（v1.2.0，Python）继续负责 pd_tmt / diann / generic 三种模式 + 8 物种 + 交互式报告；
  本子模块保留 chenyuanbao 原有的 bash+R 用法，供需要该风格输出（t-SNE、层次聚类、
  通路蛋白热图、Reactome 富集等）的场景使用。

## 与主流程的方法学差异（重要）

| 维度 | 主流程 proteomics.py | 本子模块 alt_chen |
|---|---|---|
| 缺失值 | **不填补**，limma 直接处理 | **线性最小值/2 填补**后再 log2 |
| 显著阈值 | BH 校正 P≤0.05 且 FC≥1.2 | **原始 pvalue<0.05** 且 \|log2FC\|≥log2(1.2) |
| 蛋白→基因 | 蛋白组独立，仅唯一 gene 进入富集 | 按 Genes 首个符号，同基因多蛋白**取均值** |
| 物种 | 8 种 | **3 种**（human/mouse/rat）|
| 富集库 | KEGG + GO(UniProt) + PPI | KEGG + GOBP + **Reactome**（msigdbr）|
| 配对/批次协变量 | 支持 | 不支持 |
| 报告 | 交互式(Plotly) + 静态 | 单一静态网页报告 |

同一份数据两者差异数量不同是**结果不同、并非对错**：
- TMT HepG2 实测：主流程 BH 显著 31 个（up 5 / down 26）；本子模块未校正 p 显著 40 个（up 5 / down 35）。

## 环境（固定，继承自 chenyuanbao）

- R 4.4.3：`/media/nfs/nfs02/chenyuanbao/miniforge3/envs/R433/bin/Rscript`
- Python 3.11：`/media/nfs/nfs02/chenyuanbao/miniforge3/envs/py311/bin/python`
- ImageMagick：`/usr/bin/convert`（QC 图 pdf→png）
- R 包：limma / fgsea / GSVA / clusterProfiler / ComplexHeatmap / msigdbr /
  org.Hs/Mm/Rn.eg.db / ape / gmodels / ggsci 等（R433 内已装齐，本子模块不升级不安装）

## 用法

```bash
BASEDIR=/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/proteomics/alt_chen

# 直接跑（默认 config.sh）；分步或全跑
bash "$BASEDIR/run_pipeline.sh" all
bash "$BASEDIR/run_pipeline.sh" step1

# 用自己的配置
CONFIG=/你的/config.sh bash "$BASEDIR/run_pipeline.sh" all
```

配置项见 `config.sh`（含全部注释）。关键变量：
`PP_SPECIES`（human/mouse/rat）、`PP_GROUPS`（对照在前）、`PP_REPLICATES`、
`PP_CONTRASTS`（`A_vs_B`=A−B）、`PP_MATRIX`、`PP_OUTPUT`、阈值 `PP_PCUT/PP_FCCUT`。

## 输入格式（两种）

1. **DIA-NN 蛋白矩阵 TSV**（`PP_MATRIX` 直接指定）：前 6 列注释
   （`Protein.Group/Protein.Names/Genes/First.Protein.Description/N.Sequences/N.Proteotypic.Sequences`）
   + 定量列。样本列名须能被 `{group}_{region}-{rep}` 规则重排（见 `01_process_standard.R`）。

2. **xlsx 蛋白总定量表**（`PP_MATRIX_XLSX` + `PP_INPUT_FORMAT`）：
   - `PP_INPUT_FORMAT=pd_tmt`：PD 蛋白级 TMT 总量表 → 由新增适配器
     `bin/00_convert_pd_tmt.py` 转换成标准 TSV（本子模块新增）。
   - `PP_INPUT_FORMAT=ima94`：基因水平定量表 → 由原 `bin/00_convert_ima94.py` 转换。

## 本子模块相对原 chenyuanbao 流程的改动

- 新增 `bin/00_convert_pd_tmt.py`：PD 蛋白级 TMT xlsx 适配器（原流程无此能力）。
- `run_pipeline.sh` 的 `step0_convert` 增加 `PP_INPUT_FORMAT` 分支（`pd_tmt`/`ima94`）。
- 其余 `bin/`、`lib/`（3 物种 GMT + id_map + Arial 字体 + 热图辅助）未改动。

## 真实数据实测（2026-10-20）

用 `/media/nfs/nfs02/yuzhij/Analysis/Proteomics/20250516_TMT定量蛋白质组学_单兆良_β香茅醇-NASH.zip`
的 `蛋白总定量表.xlsx`（页签 `Proteins(PO+β香茅醇 vs PO)`，人 hsa，2 组×3 重复）跑通 step0→step8：

- 7608 蛋白 → 过滤 7392 → 基因水平 7333。
- 差异（p<0.05 & |log2FC|≥log2(1.2)）：40 个（up 5 / down 35）。
- ORA：KEGG 1 条 / GOBP 34 条 / Reactome 3 条；GSEA：KEGG 显著 17、GOBP 显著 555；
  GSVA、通路蛋白热图、GSEA 领头热图、网页报告均正常生成。

适配要点（可复用）：`PP_INPUT_FORMAT=pd_tmt`，`PP_GROUPS="PO Treatment"`（对照 PO 在前），
`PP_REPLICATES=3`，`PP_CONTRASTS="Treatment_vs_PO"`，`PP_SHEET` 指向全量定量页签。

## 局限

- 仅 human/mouse/rat 三物种；猫（fca）等不在支持范围，富集/GSEA/GSVA 会因缺物种基因集失败。
- 不支持配对设计与批次协变量；样本映射为严格的 `REGIONS × GROUPS × REP`。
- 显著性用未校正 pvalue（非 BH），且分析前填补缺失值；与主流程口径不同，引用结果时需注明。
