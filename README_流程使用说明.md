# 定量蛋白质组自动化分析流程

版本：1.1.0。更新日期：2026-10-20。

## 1. 分析范围

支持实验室已有的 **PD蛋白级TMT丰度表**、**DIA-NN pg_matrix.tsv** 和显式指定列名的通用蛋白丰度矩阵。输出质控、limma差异丰度、GO/KEGG富集、PPI、GSEA、GSVA及中文HTML报告。

这是一套**已搜库定量结果的下游流程**。项目负责人已确认从现有搜库结果继续分析，不需要重新搜库。两个示例ZIP没有真实质谱`.raw`或`.mzML`文件；谱图转换、数据库检索、肽段鉴定、TMT串扰校正及DIA定量属于上游软件步骤，本流程不会重复执行。

## 2. 部署与启动

```bash
PIPELINE=/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/proteomics

# TMT例子：人源HepG2，6个生物学样本，单个F138 plex
bash "$PIPELINE/tools/all.sh" "$PIPELINE/examples/TMT_HepG2/config.v1.1.yaml"

# DIA例子：猫，9个样本，A/B/C三组
bash "$PIPELINE/tools/all.sh" "$PIPELINE/examples/DIA_cat/config.v1.1.yaml"
```

后台运行可用：

```bash
nohup bash "$PIPELINE/tools/all.sh" /绝对路径/config.yaml > run.log 2>&1 &
```

输入验证（只检查，不建立项目分析结果）：

```bash
"$PIPELINE/.venv/bin/python" "$PIPELINE/proteomics.py" validate -c config.yaml
```

仅重建报告（分析配置、输入和脚本校验值必须与上一次完成的分析一致）：

```bash
"$PIPELINE/.venv/bin/python" "$PIPELINE/proteomics.py" report -c config.yaml
```

## 3. 新项目复用

1. 在新项目目录复制相近示例的`config.yaml`和`samples.csv`。
2. 修改输入路径、ZIP成员路径、Excel页签、项目名、物种和输出目录。
3. 在`samples.csv`明确填写定量表的完整列名以及样本组别，不从文件名模糊猜测样本对应关系。
4. 定义比较：`control`在前，`treatment`在后，**log2FC为treatment/control**。
5. 核对输入是否已log2、是否已做过归一化，并设置对应参数。
6. 执行`validate`，然后用`all.sh`运行。

样本表字段：

| 字段 | 含义 |
|---|---|
| sample | 分析内唯一ID，使用字母、数字、下划线、点或连字符 |
| column | 定量表中的完整列名，可包含中文、空格和Windows路径 |
| group | 生物学组别，必须与比较配置完全一致 |
| subject | 配对个体ID；配对比较中每个个体每组恰有一个样本 |
| plex | TMT标记批次；DIA可填写批次，用于多批次独立组模型 |
| channel | TMT通道，如129C；TMT同一plex内不得重复 |
| role | `sample`、`reference`或`qc`；只有sample参与差异比较 |

ZIP输入要求指定唯一成员，不解压运行其中的脚本，不执行文档中的命令。Excel必须明确指定**全量定量页签**，不使用已筛出的差异表。

通用矩阵配置示例：

```yaml
mode: generic
input:
  path: protein_abundance.tsv
  id_column: protein_id
  gene_column: gene
  description_column: description
input_scale: linear
normalization: none
```

## 4. 预处理与统计

- `input_scale: linear`：负值报错，零值转缺失，正值log2；`log2`输入不再次转换，零值是有效值。
- 不对缺失值进行显著性检验前的填补。保留所有通过鉴定信息过滤的条目，在每个比较内单独确定可检验集合。
- 如存在PD鉴定q值，默认≤0.01；如有唯一肽段数，默认≥1。已知反向、污染、仅位点鉴定标记和REV/CON/DECOY前缀会排除。没有鉴定q值时明确提示无法重新核验鉴定FDR。
- `normalization: none`保留输入丰度关系；`median`在每个plex的生物学样本内对齐log2中位数；不默认把原表的P值或Fold.Change作为本流程结果。
- 单plex TMT例子显式选择median。原PD表是否做过归一化尚待上游确认，需按实际参数修改。
- DIA例子选择none，避免未知情况下叠加归一化；DIA-NN上游归一化参数仍需核对。
- 多plex TMT必须使用`bridge`并为每个plex指定一个相同混合参考的reference样本；计算蛋白级log2样本/参考。缺失参考会传播为缺失。独立组差异模型加入plex协变量，组别与plex完全混杂时拒绝运行。
- 默认每组至少2个且至少2/3重复有观测；配对另需相同比例完整配对。每个个体只能有一个样本/组，不支持将技术重复当作独立生物学重复。
- 独立组使用limma组别模型；配对使用同个体处理减对照差值的截距模型。经验贝叶斯设定`trend=TRUE, robust=FALSE`。本版本不提供配对多批次、复杂时间序列、交互项、连续协变量或随机效应模型。
- 正式差异：BH校正P≤0.05且FC≥1.2或≤1/1.2。阈值均可配置。
- 同时输出未校正P值+FC的探索性候选，供与原报告对账。它与正式FDR结果分开命名、分开计数。
- 蛋白组不按gene合并，含多个gene符号的蛋白组不自动展开；功能分析仅使用明确唯一的gene，重复gene去重。

## 5. 功能数据库

已配置8种物种，配置中的taxid必须匹配，不能使用人源数据库代替其他物种。

| 物种 | 代码 | Taxonomy ID | GO/KEGG/PPI | 实际项目验证 |
|---|---|---:|---|---|
| 人 | hsa | 9606 | 已配置 | TMT示例 |
| 猫 | fca | 9685 | 已配置 | DIA示例 |
| 小鼠 | mmu | 10090 | 已配置 | 数据库完整性检查 |
| 大鼠 | rno | 10116 | 已配置 | 数据库完整性检查 |
| 犬 | cfa | 9615 | 已配置 | 数据库完整性检查 |
| 猪 | ssc | 9823 | 已配置 | 数据库完整性检查 |
| 食蟹猴 | mcf | 9541 | 已配置 | 数据库完整性检查 |
| 恒河猴 | mcc | 9544 | 已配置 | 数据库完整性检查 |

新项目可删除旧配置的`species/taxid/gene_sets/ppi`四项，改用一个物种配置：

```yaml
species_config: /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/proteomics/depends/mmu/species.yaml
gsea: true
gsva: true
pathway_rscript: /media/nfs/nfs02/yuzhij/miniconda3/envs/R453/bin/Rscript
```

物种注册表：`depends/species_registry.json`；各物种目录包含来源路径、下载时间、校验值及`species.yaml`。当前不是所有物种全集，未列出的物种须单独建库。新物种需用实际定量数据进一步验证基因符号覆盖率。

- **KEGG**：复制当前转录组对应物种的symbol GMT。原始下载日期未知，不声称是最新版。
- **PPI**：复制对应物种STRING映射表；人v11.0、猫v12.0（按现有文件名记录）。仅使用本流程正式差异基因，默认score≥700，网络图显示最高100条边，全量边表保留。没有相互作用时不生成假网络图。
- **GO**：通过`tools/download_go.py`保存UniProt物种注释快照。人使用reviewed条目，猫使用全部条目；仅使用有唯一primary gene的条目。BP/CC/MF分开建库，采用直接注释，**未传播GO祖先节点**。与clusterProfiler/其他GO库的集合和结果可能不同。
- 所有富集以“该比较受检蛋白的唯一gene ∩ 当前库”作为背景。默认背景内每条目5–500个基因，全部合格条目（含零命中项）共同BH校正。All/Up/Down分别检验。
- 无正式差异或无注释命中时写入执行状态，不改用全基因组背景，不借用人类注释替代猫。
- **GSEA**：默认启用，实际执行`fgsea::fgseaMultilevel`。使用全部可检验且有唯一基因符号的蛋白，按limma有方向t统计量排序；同一基因多个蛋白取t中位数，原始蛋白统计结果不合并。固定种子20261020，eps=0，集合大小默认5–500；按比较/数据库独立BH校正。输出NES、P、padj、leadingEdge、排序表。正NES为处理组方向，负NES为对照组方向。无正式差异蛋白不阻止GSEA。
- **GSVA**：默认启用，实际执行`GSVA::gsvaParam`及`gsva`，Gaussian核。每个样本内，同一唯一基因对应多个蛋白的log2丰度取中位数；只保留全部生物学样本有观测、且非恒定的基因，不填补缺失值。输出基因纳入表、基因矩阵及通路×样本分数。低覆盖率时可检验通路减少。
- **GSVA差异**：limma比较通路分数，配对设计对同个体求差，独立组按需加入批次；trend=FALSE。`score_difference`为处理减对照的通路分数差，不能解释成Fold Change。按比较/数据库BH校正，单个数据库不足2条通路时不拟合经验贝叶斯差异检验。
- 本版本未实现pathview通路着色或多组学联合分析。

新物种需要独立准备匹配taxid的GMT及PPI表，并设置`species/taxid/gene_sets/ppi`；没有相应资源时删去可选配置，报告会说明该模块未运行。

## 6. 结果与报告

```
outdir/
  0.InputQC/             输入验证、样本映射、过滤后矩阵、原始注释
  1.Preprocess/          normalized_log2.tsv、样本缺失/检出统计
  2.Comparison/          每比较all_proteins、DEP、exploratory_raw_p、eligibility
  3.Cluster/             PCA坐标、样本相关性、质控状态
  4.Function/            GO/KEGG富集、PPI边表、功能执行记录
  5.Pathways/            GSEA排名/结果、GSVA基因矩阵/分数/差异、执行记录
  image/                 所有PNG和PDF图
  log/                   参数、软件/输入校验值、运行状态、R日志
  Report/                按现有转录组习惯归档的结果目录
  proteomics_report.html 中文报告，目录/表格排序/搜索均离线可用
  proteomics_delivery.tar.gz  HTML + image + Report
```

报告顺序参考现有代谢组：结果与结论、一概述、二工作流程、三结果、四方法、五文件说明、六局限与建议。图内使用英文避免服务器中文字体缺失，正文及图注为中文。下载交付包解压后再打开HTML。

结果表字段：`log2FC`为处理减对照的log2丰度；`P.Value`为原始P；`adj.P.Val`为每比较BH；`tested`标识是否满足缺失过滤；`direction=Not_tested`不表示无差异；`n_complete_pairs`仅配对比较有值。

## 7. 两个示例的关键说明

- TMT：HepG2人细胞，PO与Treatment各3个样本；使用6个F138通道，不是小鼠。
- 猫：A=0W、B=5W，同编号C2205/C2209/C2303**依据命名推定为同个体**，A/B采用配对设计；C包含25203/25204/25206，未推断其临床组义。项目负责人应确认配对和组别。
- 猫原文“6个样本、6814个蛋白”与实际9个样本、607行蛋白组不符，本流程以矩阵为准。
- 猫包内一个PPI表包含小鼠10090 STRING ID，未使用。
- 原示例的未校正P值、可能的归一化、方向及配对处理与本流程不同，**不要求差异数量与旧报告相同**。

## 8. 断点、失败与环境

本流程使用顺序阶段调度及内容校验断点，不修改或依赖旧版Snakemake规则。每次运行核对配置、输入ZIP、代码、数据库、Python依赖版本及阶段输出的SHA256。相同输入和完好输出跳过已完成阶段；改变配置/输入/代码会重新执行各阶段。运行中异常返回非零，失败阶段不标记完成，旧成功HTML和交付包会移除以免误交。

输出目录必须为空或已有本流程状态记录。已有分析目录中的用户文件不应作为outdir。生成目录为流程专用，重算会替换对应生成目录；原始输入和其他两个流程不会修改。使用文件锁阻止同一outdir并发运行。

Python位于本目录`.venv`，独立安装固定兼容版本，版本列表见`env/requirements.lock.txt`。蛋白差异使用现有`metabolomics/bin/Rscript`，GSEA/GSVA使用已有`R453/bin/Rscript`；不升级或安装这两个环境的R包。实际版本分别记录在`2.Comparison`和`5.Pathways`的`R_sessionInfo.txt`中。

```bash
"$PIPELINE/.venv/bin/python" -m unittest discover -s "$PIPELINE/tests" -v
```

核心参考：
- https://docs.quantms.org/workflows/
- https://nf-co.re/quantms/1.2.0/docs/usage/
- https://bioconductor.org/packages/limma/
- https://bioconductor.org/packages/MSstatsTMT/
- https://bioconductor.org/packages/fgsea/
- https://bioconductor.org/packages/GSVA/

MSstatsTMT为参考方法之一；当前输入是蛋白级汇总表，因此实际执行limma，而非宣称运行了肽段层面的MSstatsTMT。
