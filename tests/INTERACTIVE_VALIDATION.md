# 交互报告验收

服务器生产运行只需要原有Python/R环境及随源码交付的Plotly.js，不需要Node.js。

本测试针对TMT_HepG2、DIA_cat两套固定验收项目。测试机需要Node.js 22和jsdom 26.1.0：

```bash
npm ci --prefix tests --ignore-scripts
node tests/verify_interactive.cjs /包含TMT_HepG2和DIA_cat两个目录的路径
```

每个目录中应包含从服务器下载的最终`proteomics_report.html`；测试首先核对HTML中嵌入的脚本与当前模板一致，随后验证9页DOM、点击回调、跨页蛋白/样本联动、比较同步、缺失及空结果、GSEA曲线ES、GSVA原始分数和导出回调。结果写入指定目录的`interaction_validation.json`。

Plotly渲染接口被模拟，实际库仅用于部分图形参数schema验证。Plotly 3.1.0验证器在不显示图例时校验图例子属性会抛异常，因此该检查略去图例样式、隐藏坐标轴样式和unused提示。该测试不等价于真实浏览器布局、图像下载或最终视觉验收。

分析数值和交付包检查在服务器运行：

```bash
.venv/bin/python tools/verify_results.py
.venv/bin/python -m unittest discover -s tests -v
```

`verify_results.py`针对两套v1.2固定配置，将标准化矩阵和蛋白差异表的SHA256与v1.1比较，并独立复核FC、BH校正、通路排名、GSVA差值及文件链接。GSEA全部曲线的ES与fgsea结果核验在每次构建报告时执行。
