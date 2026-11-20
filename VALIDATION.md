# 部署与验收记录

日期：2026-10-20。主机：VScode，账户yuzhij。

部署位置：`/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/proteomics`

## v1.2交互报告验收

项目负责人已查看两页交互原型并确认样式。沿用该样式扩展9个页面：差异蛋白、样本/PCA、质控、蛋白热图、GO/KEGG、GSEA、GSVA、PPI、汇总与文件。主流程报告阶段自动生成`proteomics_report.html`；原章节式HTML另存为`proteomics_static_report.html`，两者随交付包保存。

- 两套原始ZIP重新完成7个分析阶段，8项原有测试全部通过。
- 两例标准化矩阵和全部蛋白差异表与v1.1逐文件SHA256完全一致。正式差异数量仍为TMT 31个、猫三个比较均0个。
- 全部4104条GSEA曲线ES与fgsea原表一致（人3948条，猫156条）；样本丰度和GSVA显示值直接来自已保存结果。
- 两套最终HTML与当前模板一致，9页DOM及事件回调验证通过，包括火山图/PCA点击、蛋白和样本跨页跳转、GSEA命中基因、GSVA、PPI、筛选、比较同步和导出回调。使用jsdom和模拟Plotly渲染器，实际Plotly库作部分schema检查。
- 交互HTML无外部脚本、样式依赖。结果文件和静态报告链接均存在，交付包解压通过。
- 猫没有正式差异蛋白时，默认差异热图与PPI显示空结果；高变集合明确标为探索性展示。
- 两例再次运行均显示7阶段RESUME。

机器记录：`validation/acceptance.v1.2.json`、`validation/interaction.v1.2.json`。可复现的交互测试见`tests/INTERACTIVE_VALIDATION.md`。

输出：`/media/nfs/nfs02/yuzhij/Analysis/Proteomics/pipeline_runs_v1.2_20261020/{TMT_HepG2,DIA_cat}`。旧版本与两页原型保留。

**验收界限**：两页原型样式已由项目负责人查看确认；新增完整报告完成数据、链接和事件逻辑验证。工具受浏览器本地文件访问策略限制，未对完整9页进行真实浏览器渲染、截图或图片下载验收，不将模拟测试描述为浏览器视觉测试。

## v1.1物种与通路模块验收

已按项目负责人确认的范围继续分析已搜库定量结果，不加入原始质谱重搜库。

- 人、猫、小鼠、大鼠、犬、猪、食蟹猴、恒河猴共8种，均具备GO/KEGG/PPI及可复用species.yaml。8种资源的物种编号、GMT非空、PPI格式/分数范围、KEGG与PPI基因符号交集全部通过检查。
- 新增6种完成数据库检查，尚无实际项目数据验收；人和猫完成真实项目全流程验证。不称作“所有物种齐全”。
- GSEA实际使用fgsea 1.36.2；GSVA实际使用GSVA 2.4.9；通路环境为R453（R 4.5.3），BiocParallel 1.44.0，limma 3.66.0。
- 两套真实数据全部7个阶段完成。GSEA排名与蛋白t统计量的基因中位数一致；GSEA及GSVA BH结果独立复核通过；GSVA处理/对照分数差方向复核通过。
- 人TMT：GSEA/GSVA纳入6865个基因；KEGG受检231条，GSEA显著4条、GSVA差异显著11条（BH≤0.05，各模块独立校正）。
- 猫DIA：各比较GSEA纳入244/263/284个基因；GSVA完整观测筛选后172个基因，KEGG仅10条通路。注释与缺失覆盖限制需随结果解释，不能将低覆盖率视作全蛋白组通路结论。
- 模拟数据整合测试通过：预设上调/下调通路的GSEA方向、反向比较的NES符号、GSVA配对差值反向、缺失基因排除、重复基因中位数聚合。
- 原8项单元/统计测试再次全部通过。两例断点恢复均显示7阶段RESUME。
- 新报告链接无缺失，GSEA柱图、GSVA热图静态目视检查通过；浏览器交互验证限制仍同下述说明。

新版输出：`/media/nfs/nfs02/yuzhij/Analysis/Proteomics/pipeline_runs_v1.1_20261020/{TMT_HepG2,DIA_cat}`，旧结果保留。

机器可读记录：`validation/acceptance.v1.1.json`、`validation/species.v1.1.json`、`validation/pathway_synthetic.v1.1.json`。

## 实际验证

两个原始ZIP只读使用，未修改转录组/代谢组流程及其环境。

| 示例 | 输入蛋白组 | 生物学样本 | 比较 | 受检蛋白 | BH显著上调 | BH显著下调 | 未校正P候选 |
|---|---:|---:|---|---:|---:|---:|---:|
| 人HepG2 TMT | 7608 | 6 | PO_vs_Treatment | 6926 | 5 | 26 | 57 |
| 猫DIA | 607 | 9 | A_vs_B（配对） | 409 | 0 | 0 | 25 |
| 猫DIA | 607 | 9 | A_vs_C | 433 | 0 | 0 | 39 |
| 猫DIA | 607 | 9 | B_vs_C | 457 | 0 | 0 | 30 |

统计阈值：BH校正P≤0.05且FC≥1.2或≤1/1.2。TMT鉴定q值过滤后保留7019个蛋白，随后在比较内执行缺失过滤。未校正候选不计入正式显著结果。零显著不证明不存在生物学效应。

已检查：

- 8项单元/统计测试全部通过，覆盖重复ID、样本列匹配、非法丰度、零/缺失处理、歧义gene、配对完整性、物种匹配、BH及配对方向。
- 真实数据的log2FC与独立计算的组间均值差/配对差值均值一致。
- R产生的BH校正P与Python独立实现一致。
- 报告差异数量、全量结果表和差异表一致。
- 两份报告的本地下载/图片链接无缺失，交付包可解压。
- 重复运行全部阶段显示RESUME；代码/模板改动后触发重算。
- 静态火山图、PCA、GO图、PPI图经过目视检查；小网络布局已修正。
- 报告表格JavaScript语法检查通过。

**验证限制**：浏览器本地文件访问策略阻止HTML页面预览，未进行真实浏览器交互和最终页面截图验证；已完成HTML链接、数据和静态图片检查。

## 运行环境

- CentOS Linux 7。
- Python 3.12独立venv，核心固定numpy 1.26.4、pandas 2.2.3、scipy 1.13.1、matplotlib 3.9.2、openpyxl 3.1.5、PyYAML 6.0.2、networkx 3.4.2。
- R 4.5.3，limma 3.66.0，jsonlite 2.0.0；读取现有metabolomics环境，未修改R包。
- 完整Python锁定列表：`env/requirements.lock.txt`。
- 自动验收结果：`validation/acceptance.json`；测试日志：`validation/tests.log`。

## 需要项目负责人核实的分析假设

1. 猫A/B同编号的0W与5W是否确为同一个体；当前按配对分析。
2. A/B/C的正式生物学名称，当前保留原组名。
3. PD及DIA-NN上游归一化参数；示例TMT显式median，DIA不叠加归一化。
4. 输入矩阵无法重建原始谱图鉴定FDR、PSM/肽段完整质控、TMT同位素串扰等信息。

完整使用说明见`README_流程使用说明.md`。
