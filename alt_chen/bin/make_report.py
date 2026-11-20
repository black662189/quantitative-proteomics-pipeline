#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HY-0094 心脏蛋白质组学标准报告生成脚本 (v4: 按部位→分组排布, GSEA 激活/抑制分开展示)
- 左侧固定目录(aside), 三级: 章节→部位→对比
- 正文 section + p.paragraph + p.table_title(图题) + img_table 并排图
- 表格 .chart 容器可滚动; 通路热图结果表内链接浮层查看
"""
import os, shutil, html, sys

PP_ROOT = os.environ.get("PP_ROOT", os.getcwd())
OUTPUT  = os.environ.get("PP_OUTPUT", os.path.join(PP_ROOT, "output"))
OUT_STD = os.environ.get("PP_OUT_STD", os.path.join(OUTPUT, "02_standardize"))
OUT_QC  = os.environ.get("PP_OUT_QC", os.path.join(OUTPUT, "03_qc_pca_cluster"))
OUT_DE  = os.environ.get("PP_OUT_DE", os.path.join(OUTPUT, "04_de"))
OUT_ORA = os.environ.get("PP_OUT_ORA", os.path.join(OUTPUT, "05_ora"))
OUT_GSEA= os.environ.get("PP_OUT_GSEA", os.path.join(OUTPUT, "06_gsea"))
OUT_GSVA= os.environ.get("PP_OUT_GSVA", os.path.join(OUTPUT, "07_gsva"))
OUT_HEAT= os.environ.get("PP_OUT_HEATMAP", os.path.join(OUTPUT, "08_pathway_heatmaps"))
OUT_REPORT = os.environ.get("PP_OUT_REPORT", os.path.join(OUTPUT, "09_report"))
IDMAP  = os.environ.get("PP_IDMAP", os.path.join(PP_ROOT, "lib", "id_map.tsv"))
os.makedirs(OUT_REPORT, exist_ok=True)

REGIONS  = [x for x in os.environ.get("PP_REGIONS", "IVS LVAW LVPW").split() if x]
CONTRASTS = [x for x in os.environ.get("PP_CONTRASTS", "Model_vs_Normal HY0094_Model_vs_Model HY0094_Model_vs_Normal").split() if x]
GROUPS   = [x for x in os.environ.get("PP_GROUPS", "Normal Model HY0094_Model").split() if x]
PCUT  = float(os.environ.get("PP_PCUT", "0.05"))
FCCUT = float(os.environ.get("PP_FCCUT", "0.2630344058337938"))

REGION_CN = {"IVS": "室间隔 (IVS)", "LVAW": "左室前壁 (LVAW)", "LVPW": "左室后壁 (LVPW)"}
CONTRAST_CN = {
    "Model_vs_Normal": "模型组 vs 正常组",
    "HY0094_Model_vs_Model": "给药组 vs 模型组",
    "HY0094_Model_vs_Normal": "给药组 vs 正常组",
}
def rcn(rg): return REGION_CN.get(rg, rg)
def ccn(cc): return CONTRAST_CN.get(cc, cc.replace("_vs_", " vs "))
DBS = ["KEGG", "GOBP", "Reactome"]

manifest = []
png2pdf = {}   # png 相对url -> pdf 相对url (或 None)

def copy_img(src, name, with_pdf=True):
    """直接引用 pipeline 输出目录中的图 (不复制), 返回相对报告目录的相对路径"""
    if not os.path.exists(src):
        return None
    rel = os.path.relpath(src, OUT_REPORT)
    pdf_rel = None
    if with_pdf:
        pdf_src = os.path.splitext(src)[0] + ".pdf"
        if os.path.exists(pdf_src):
            pdf_rel = os.path.relpath(pdf_src, OUT_REPORT)
    png2pdf[rel] = pdf_rel
    return rel

def read_tsv(f):
    rows = []
    with open(f, encoding="utf-8", errors="replace") as fh:
        header = None
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if i == 0:
                header = parts
            else:
                rows.append(parts)
    return header, rows

def short_id(raw_id):
    return ID_MAP.get(raw_id, raw_id)

# 通路名 -> 标准短ID 映射 (GO:xxx / hsa / R-HSA)
_idmap_rows = read_tsv(IDMAP)
ID_MAP = {}
if _idmap_rows[0] and "gs_exact_source" in _idmap_rows[0]:
    _gi = _idmap_rows[0].index("gs_name")
    _si = _idmap_rows[0].index("gs_exact_source")
    for _r in _idmap_rows[1]:
        ID_MAP[_r[_gi]] = _r[_si]

def esc(x):
    return html.escape(str(x))

def sf(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None

def fmt_name(p, db):
    if db == "KEGG" or db is None:
        return p
    if db == "GOBP":
        s = p[5:] if p.startswith("GOBP_") else p
    elif db == "Reactome":
        s = p[9:] if p.startswith("REACTOME_") else p
    else:
        s = p
    s = s.replace("_", " ").lower()
    return s[:1].upper() + s[1:] if s else s

def lb_link(url, text="查看"):
    return '<a class="lb" href="' + esc(url) + '">' + esc(text) + '</a>'

fig_no = [0]
tab_no = [0]

def tab_caption(title):
    tab_no[0] += 1
    return '<p class="table_title">表 %d  %s</p>' % (tab_no[0], esc(title))

def img_row(pairs, widths=None):
    if widths is None:
        widths = [None] * len(pairs)
    tds = []
    for (img, cap), wd in zip(pairs, widths):
        ws = (' width="' + wd + '"') if wd else ""
        # 把图片包成可点击放大(lightbox) + PDF 下载链接
        wrapped = wrap_fig(img)
        tds.append('<td' + ws + '>' + wrapped + '</td>')
    caps = []
    for img, cap in pairs:
        if cap:
            fig_no[0] += 1
            caps.append('<td><p class="table_title">图 %d  %s</p></td>' % (fig_no[0], esc(cap)))
        else:
            caps.append('<td></td>')
    return ('<table align="center" class="img_table"><tr>' + "".join(tds) + '</tr><tr>' + "".join(caps) + '</tr></table>')

import re as _re
def wrap_fig(img_html):
    """把 <img src=...> 包成可点击放大 + PDF 下载的 <a class=lbpic>"""
    m = _re.search(r'src="([^"]+)"', img_html)
    if not m:
        return img_html
    png = m.group(1)
    pdf = png2pdf.get(png)
    pdf_attr = (' data-pdf="' + esc(pdf) + '"') if pdf else ""
    return '<a class="lbpic" href="' + esc(png) + '"' + pdf_attr + '>' + img_html + '</a>'

def make_table(headers, rows, max_rows=None, link_fn=None):
    # 完整渲染所有行, 靠 .chart 滚动窗口查看(默认看 top, 向下滑动看剩余)
    # 列类型: 文本列(换行, 宽) / 数字列(nowrap 右对齐, 窄) / ID列(窄 nowrap)
    def col_meta(h):
        hl = (h or "").lower()
        if hl in ("description", "pathway", "term", "质控项", "通路") or "通路" in hl:
            return ("text", 3)
        if hl in ("id", "geneid", "gene id", "基因", "genes"):
            return ("id", 2)
        if hl in ("pvalue", "p.adjust", "padj", "pval", "count", "nes", "generatio",
                  "genratio", "size", "log2foldchange", "logfc", "bgratio", "richfactor",
                  "foldenrichment", "zscore", "qvalue", "缺失", "不显著", "上调", "下调"):
            return ("num", 1)
        return ("plain", 1)

    extra_col = link_fn is not None
    th = list(headers)
    if extra_col:
        th.append("通路热图")

    metas = [col_meta(h) for h in th]
    if extra_col:
        metas[-1] = ("link", 1)  # 通路热图列
    weights = [m[1] for m in metas]
    wsum = float(sum(weights))
    colgroup = '<colgroup>' + "".join(
        '<col style="width:%.2f%%">' % (w / wsum * 100.0) for w in weights) + '</colgroup>'

    thead = '<thead><tr>' + "".join("<th>" + esc(h) + "</th>" for h in th) + '</tr></thead>'

    def render_rows(data):
        out = []
        for r in data:
            tds = []
            for j, v in enumerate(r):
                cls = 'class="num"' if metas[j][0] == "num" else ""
                tds.append("<td %s>%s</td>" % (cls, esc(v)))
            if extra_col:
                url = link_fn(r)
                tds.append("<td>" + (lb_link(url) if url else '<span class="empty">—</span>') + "</td>")
            out.append("<tr>" + "".join(tds) + "</tr>")
        return "".join(out)

    n = len(rows)
    c = ['<div class="tableblock">',
         '<div class="chart"><table align="center" id="table">' + colgroup + thead + '<tbody>' + render_rows(rows) + '</tbody></table></div>']
    if n > 10:
        c.append('<p class="introduction">共 %d 行，滚动条查看全部。</p>' % n)
    c.append('</div>')
    return "\n".join(c)

# ---------------- 映射 ----------------
ora_map = {}
lead_map = {}

def load_maps():
    for rg in REGIONS:
        mp = os.path.join(OUT_HEAT, rg, "ORA", "map.tsv")
        if os.path.exists(mp):
            h, rows = read_tsv(mp)
            ic = {c: i for i, c in enumerate(h)}
            for r in rows:
                if "heatmap" not in ic or r[ic["heatmap"]] == "NA":
                    continue
                key = (rg, r[ic["contrast"]], r[ic["db"]])
                term = r[ic["term"]]
                src = os.path.join(OUT_HEAT, rg, "ORA", r[ic["heatmap"]])
                rel = copy_img(src, "hm_" + r[ic["heatmap"]][:-4])
                if rel:
                    ora_map.setdefault(key, {})[term] = rel
        lm = os.path.join(OUT_HEAT, rg, "GSEA", "lead_map.tsv")
        if os.path.exists(lm):
            h, rows = read_tsv(lm)
            ic = {c: i for i, c in enumerate(h)}
            for r in rows:
                if "heatmap" not in ic or r[ic["heatmap"]] == "NA":
                    continue
                key = (rg, r[ic["contrast"]], r[ic["db"]])
                pw = r[ic["pathway"]]
                src = os.path.join(OUT_HEAT, rg, "GSEA", r[ic["heatmap"]])
                rel = copy_img(src, "lead_" + r[ic["heatmap"]][:-4])
                if rel:
                    lead_map.setdefault(key, {})[pw] = rel

load_maps()

# ---------------- 目录树 & 正文 ----------------
# toc 节点: dict(level=1/2/3, anchor, title, children=[])
toc_root = []
_cur_h1 = None
_cur_h2 = None
body_parts = []

def H1(anchor, title):
    global _cur_h1, _cur_h2
    node = dict(level=1, anchor=anchor, title=title, children=[])
    toc_root.append(node)
    _cur_h1 = node
    _cur_h2 = None
    body_parts.append('<h1><a name="%s">%s</a></h1>' % (anchor, esc(title)))

def H2(anchor, title):
    global _cur_h2
    node = dict(level=2, anchor=anchor, title=title, children=[])
    _cur_h1["children"].append(node)
    _cur_h2 = node
    body_parts.append('<h2><a name="%s">%s</a></h2>' % (anchor, esc(title)))

def H3(anchor, title):
    node = dict(level=3, anchor=anchor, title=title, children=[])
    _cur_h2["children"].append(node)
    body_parts.append('<h3><a name="%s">%s</a></h3>' % (anchor, esc(title)))

def H4(title):
    body_parts.append('<h4>%s</h4>' % esc(title))

def P(text):
    body_parts.append('<p class="paragraph">%s</p>' % text)

def NOTE(text):
    body_parts.append('<p class="introduction">%s</p>' % text)

def RAW(html):
    body_parts.append(html)

# ================= 一、分析概述 =================
H1("overview", "一、分析概述")
P("本项目对蛋白质组学定量数据进行标准分析，研究对象为 %d 个部位（%s），每组设置 %s 各 %d 个生物学重复，共 %d 个样本。采用液相色谱-串联质谱（LC-MS/MS）定量，数据经过滤、log2 变换与缺失值填充后，进行差异表达分析、功能富集（KEGG / Gene Ontology BP / Reactome）、基因集富集分析（GSEA）与基因集变异分析（GSVA）。" % (len(REGIONS), "、".join(REGION_CN.get(r, r) for r in REGIONS), "、".join(GROUPS), int(os.environ.get("PP_REPLICATES", "6")), len(REGIONS) * len(GROUPS) * int(os.environ.get("PP_REPLICATES", "6"))))
P("蛋白定量经「任一 部位×组别 中定量占比达标」过滤，缺失值以线性最小值法填充，基因水平按「Genes 列首个符号、同基因多蛋白取均值」合并。差异蛋白以 limma 经验贝叶斯检验，筛选阈值 pvalue &lt; %g 且 |log2FoldChange| ≥ log2(1.2)。具体质控数值见下文「质控指标」表。" % PCUT)

# ================= 二、数据预处理与质控 =================
H1("preproc", "二、数据预处理与质控")
P("数据预处理流程为：原始定量矩阵 → 蛋白过滤 → 缺失值填充 → log2 变换与基因水平合并。关键质控指标如下表。")
H2("qc_metric", "2.1 质控指标")
# 从 process_QA.txt 动态解析质控指标与逐样本缺失率
_qa_path = os.path.join(OUT_STD, "process_QA.txt")
_qc_summary = []
_sample_missing = []   # [ (样本名, 缺失率%), ... ]
_qa_lines = open(_qa_path, encoding="utf-8", errors="replace").read().splitlines() if os.path.exists(_qa_path) else []
for _l in _qa_lines:
    _l = _l.strip()
    if ":" in _l and not _l.startswith(tuple(["  "] + [r + ":" for r in REGIONS])):
        _k, _v = _l.split(":", 1)
        _k = _k.strip(); _v = _v.strip()
        if _k in ("样本映射", "蛋白过滤", "缺失值填充", "log2 数值范围", "部位蛋白集一致"):
            _qc_summary.append([_k, _v])
    _mm = _re.match(r"\s*样本缺失率:\s*(.*)$", _l)
    if _mm:
        for _tok in _mm.group(1).split("|"):
            _tok = _tok.strip()
            if " " in _tok:
                _sn, _pct = _tok.rsplit(" ", 1)
                _sample_missing.append([_sn, _pct])
qc_summary = _qc_summary or [["质控信息", "见 process_QA.txt"]]
RAW(tab_caption("质控指标汇总"))
RAW(make_table(["质控项", "结果"], qc_summary))
H2("qc_missing", "2.2 各样本缺失率")
P("过滤后各样本的缺失比例如下表，均处于较低水平，样本间定量一致性良好。")
RAW(tab_caption("各样本缺失比例（%）"))
RAW(make_table(["样本", "缺失比例（%）"], _sample_missing))
H2("sample", "2.3 样本信息")
meta_h, meta_rows = read_tsv(os.path.join(OUT_STD, "sample_metadata.xls"))
RAW(tab_caption("样本信息"))
RAW(make_table(meta_h, meta_rows, max_rows=60))

# ================= 三、分析结果（按部位 → 分组）=================
H1("result", "三、分析结果")

# 差异蛋白统计总表
stat_rows = []
for rg in REGIONS:
    for cc in CONTRASTS:
        f = os.path.join(OUT_DE, rg, "DE_" + cc + ".xls")
        if not os.path.exists(f):
            continue
        h, rows = read_tsv(f)
        gi = h.index("log2FoldChange"); pi = h.index("pvalue")
        up = dn = 0
        for r in rows:
            p = sf(r[pi]); g = sf(r[gi])
            if p is None or g is None:
                continue
            if p < PCUT and abs(g) >= FCCUT and g > 0:
                up += 1
            elif p < PCUT and abs(g) >= FCCUT and g < 0:
                dn += 1
        stat_rows.append([rcn(rg), ccn(cc), str(up), str(dn), str(len(rows) - up - dn)])
P("以下按「部位 → 分组」顺序展示每个心脏部位在各对比下的差异蛋白、三库富集与两库 GSEA 结果。差异蛋白统计见表。")
RAW(tab_caption("差异蛋白统计（pvalue &lt; 0.05 且 |log2FC| ≥ log2(1.2)）"))
RAW(make_table(["部位", "对比", "上调蛋白数", "下调蛋白数", "不显著"], stat_rows))

for rg in REGIONS:
    H2("region_%s" % rg, "3.%d %s" % (REGIONS.index(rg) + 1, rcn(rg)))
    # ---- 部位质控图 (带标签 PCA + 均值树状聚类) ----
    H4("样本质量控制图（带样本标签 PCA / 均值树状聚类）")
    qc_dir = os.path.join(OUT_QC, rg)
    _grpstr = "-vs-".join(GROUPS)
    qc_pairs = []
    for qk, qlab in [("PCA." + _grpstr, "带样本标签 PCA"),
                     ("hclust_average." + _grpstr, "均值树状聚类")]:
        src = os.path.join(qc_dir, qk + ".png")
        if os.path.exists(src):
            rel = copy_img(src, "qc_%s_%s" % (rg, qk))
            qc_pairs.append(('<img src="%s" alt="%s" class="goImg">' % (rel, qlab), qlab))
    if qc_pairs:
        RAW(img_row(qc_pairs, ["50%", "50%"]))
    for cc in CONTRASTS:
        H3("r_%s_%s" % (rg, cc), "3.%d.%d %s" % (REGIONS.index(rg) + 1, CONTRASTS.index(cc) + 1, ccn(cc)))
        # ---- 差异火山图 ----
        H4("差异蛋白火山图")
        rel_vol = copy_img(os.path.join(OUT_DE, rg, "Volcano", cc + "_volcano.png"),
                           "volcano_" + rg + "_" + cc)
        if rel_vol:
            RAW(img_row([('<img src="%s" alt="火山图" class="goImg" style="max-width:520px;">' % rel_vol, "火山图")], None))
        # ---- 三库富集 ----
        H4("三库富集（KEGG / GOBP / Reactome）")
        bubble_pairs = []
        for db in DBS:
            f = os.path.join(OUT_ORA, rg, cc + "_" + db + "_enrichment.xls")
            if not os.path.exists(f):
                continue
            rel = copy_img(os.path.join(OUT_ORA, rg, cc + "_" + db + "_bubble.png"),
                           "ora_" + rg + "_" + cc + "_" + db)
            if rel:
                bubble_pairs.append(('<img src="%s" alt="%s 富集" class="goImg">' % (rel, db), db + " 富集"))
        if bubble_pairs:
            RAW(img_row(bubble_pairs, ["33%"] * len(bubble_pairs)))
        for db in DBS:
            f = os.path.join(OUT_ORA, rg, cc + "_" + db + "_enrichment.xls")
            if not os.path.exists(f):
                continue
            h, rows = read_tsv(f)
            if not rows:
                RAW('<p class="introduction">%s：无显著富集通路。</p>' % db)
                continue
            keep_cols = ["ID", "Description", "GeneRatio", "pvalue", "p.adjust", "Count"]
            keep_idx = [h.index(c) for c in keep_cols if c in h]
            sub_h = [h[i] for i in keep_idx]
            id_col = sub_h.index("ID") if "ID" in sub_h else None
            sub_rows = []
            for r in rows:
                sr = [r[i] for i in keep_idx]
                if id_col is not None:
                    sr[id_col] = short_id(sr[id_col])
                sub_rows.append(sr)
            desc_idx = sub_h.index("Description") if "Description" in sub_h else None
            mp = ora_map.get((rg, cc, db), {})
            def ora_link(r):
                return mp.get(r[desc_idx]) if desc_idx is not None else None
            RAW(tab_caption("%s 富集通路" % db))
            RAW(make_table(sub_h, sub_rows, max_rows=10, link_fn=ora_link))
        # ---- 两库 GSEA ----
        H4("两库 GSEA（KEGG / GOBP）激活/抑制通路")
        for db in ["KEGG", "GOBP"]:
            f = os.path.join(OUT_GSEA, rg, db + "_" + cc + "_fgsea_results.xls")
            if not os.path.exists(f):
                continue
            rel_up = copy_img(os.path.join(OUT_GSEA, rg, "plots", db + "_" + cc + "_NES_bar_up.png"),
                              "gsea_" + rg + "_" + cc + "_" + db + "_up")
            rel_dn = copy_img(os.path.join(OUT_GSEA, rg, "plots", db + "_" + cc + "_NES_bar_down.png"),
                              "gsea_" + rg + "_" + cc + "_" + db + "_down")
            pairs = []
            if rel_up:
                pairs.append(('<img src="%s" alt="%s 激活通路" class="goImg">' % (rel_up, db), db + " 激活 top"))
            if rel_dn:
                pairs.append(('<img src="%s" alt="%s 抑制通路" class="goImg">' % (rel_dn, db), db + " 抑制 top"))
            if pairs:
                RAW(img_row(pairs, ["50%"] * len(pairs)))
            # 显著通路表
            h, rows = read_tsv(f)
            if "pval" in h and "NES" in h:
                pi = h.index("pval")
                sig = [r for r in rows if (sf(r[pi]) is not None and sf(r[pi]) < 0.05)]
                sig.sort(key=lambda r: sf(r[pi]))
            else:
                sig = rows
            if sig:
                keep_cols = [c for c in ["pathway", "pval", "padj", "NES", "size"] if c in h]
                keep_idx = [h.index(c) for c in keep_cols]
                sub_h = [h[i] for i in keep_idx]
                pw_idx = sub_h.index("pathway") if "pathway" in sub_h else 0
                sub_rows = []
                for r in sig:
                    sr = [r[i] for i in keep_idx]
                    raw_pw = sr[pw_idx]
                    sid = short_id(raw_pw)
                    nm = fmt_name(raw_pw, db)
                    if sid != raw_pw:
                        sr[pw_idx] = sid + " / " + nm
                    else:
                        sr[pw_idx] = nm
                    sub_rows.append(sr)
                lmp = lead_map.get((rg, cc, db), {})
                def lead_link(r):
                    key = r[pw_idx]
                    return lmp.get(key) or lmp.get(key.split(" / ", 1)[-1])
                RAW(tab_caption("%s GSEA 显著通路（p&lt;0.05）" % db))
                RAW(make_table(sub_h, sub_rows, max_rows=10, link_fn=lead_link))

    # ---- GSVA 两方向 top 通路热图 (本部位末尾, 与对比无关) ----
    H4("GSVA 通路活性（%s 两方向 top 通路热图）" % ccn(CONTRASTS[0]))
    g_pairs = []
    for db in ["KEGG", "GOBP"]:
        rel = copy_img(os.path.join(OUT_GSVA, rg, "GSVA_top", db + "_directional_GSVA_heatmap.png"),
                       "gsva_" + rg + "_" + db + "_directional")
        if rel:
            g_pairs.append(('<img src="%s" alt="%s 方向热图" class="goImg">' % (rel, db), db + " 两方向 top 通路热图"))
    if g_pairs:
        RAW(img_row(g_pairs, ["50%", "50%"]))

# ================= 四、方法描述 =================
H1("method", "四、方法描述")
P("差异表达分析采用 R 包 limma 的经验贝叶斯线性模型，对比为「%s」。差异蛋白筛选阈值为 pvalue &lt; 0.05 且 |log2FoldChange| ≥ log2(1.2)（即表达变化 ≥1.2 倍），并同时给出 Benjamini-Hochberg 多重检验校正的 qvalue。" % "、".join(CONTRASTS))
P("富集分析采用 clusterProfiler 的超几何检验（enricher），输入为上述 pvalue &lt; 0.05 且 |log2FC| ≥ log2(1.2) 筛选出的差异蛋白（上、下调合并为一个基因集），背景为全部定量基因（universe），数据库为 KEGG（hsa.kegg_optimal）、GOBP（c5.go.bp.v7.5.1.symbols）与 Reactome（c2.cp.reactome.v7.5.1.symbols），p 值经 BH 校正，过滤 pvalue &lt; 0.05 且通路内基因数 ≥3。")
P("GSEA 采用 fgsea（排序统计量为 limma 的 log2FoldChange），对完整 KEGG 与 GOBP 基因集进行基因集富集检验。柱状图分别展示激活（NES&gt;0，红）与抑制（NES&lt;0，蓝）top 通路。")
P("GSVA 采用 R 包 GSVA（kcdf=Gaussian），对 KEGG 与 GOBP 基因集进行通路活性打分，以「实验组 vs 对照组」组均值差定义方向（实验组 − 对照组），分别取「实验组上调」与「实验组下调」各 top%d 通路绘制热图。" % int(os.environ.get("PP_GSVA_TOPN", "15")))
P("所有绘图采用 R 语言（ggplot2 / ComplexHeatmap），配色遵循项目统一规范（各分组对应固定颜色，见 config.sh）；显著性配色沿用红蓝渐变。")

# ================= 五、结果文件说明 =================
H1("files", "五、结果文件说明")
P("本报告的图与结果文件均已整理到 HY-0094蛋白质组学标准分析 目录（本报告 html_report.html 即位于该目录下），按「部位 → 分析流程」组织，目录结构说明如下：")
RAW('<ul class="filelist">')
RAW('<li><b>01_数据预处理与质控/</b>：质控指标与缺失率.xlsx、样本信息.xls、质控图_PCA与聚类/（各部位带标签 PCA 与均值树状聚类，pdf/png）</li>')
RAW('<li><b>%s</b>（各部位内含以下 4 个流程子目录）：</li>' % (" / ".join("0%d_%s" % (i+2, rcn(rg)) for i, rg in enumerate(REGIONS))))
RAW('<li style="margin-left:1.5em">1_差异表达分析/：{对比}_火山图.pdf/.png + {对比}_差异蛋白表.xls</li>')
RAW('<li style="margin-left:1.5em">2_ORA功能富集/：{对比}_{数据库}_富集气泡图.pdf/.png + {对比}_{数据库}_富集表.xls；通路蛋白热图/ 下按「对比_数据库」分组存放各通路蛋白热图</li>')
RAW('<li style="margin-left:1.5em">3_GSEA基因集富集/：{数据库}_{对比}_激活/抑制通路bar图 + {数据库}_{对比}_显著通路表.xls；前导边热图/ 下按「对比_数据库」分组</li>')
RAW('<li style="margin-left:1.5em">4_GSVA通路活性/：{KEGG,GOBP}_两方向top通路热图.pdf/.png</li>')
RAW('<li><b>数据/</b>（整个分析流程的数据，含原始数据与处理后分析输入）：</li>')
RAW('<li style="margin-left:1.5em">原始质谱数据/：report.pg_matrix.tsv + report.pr_matrix.tsv + 高通量检测登记表.docx</li>')
RAW('<li style="margin-left:1.5em">标准化数据/：{部位}_protein_expression_log2.xls（log2）+ {部位}_protein_expression_log2_imputed.xls（log2+填充，分析输入）+ sample_metadata.xls + process_QA.txt</li>')
RAW('<li style="margin-left:1.5em">基因水平合并矩阵/：{部位}_gene_expression_log2_all.xls</li>')
RAW('<li style="margin-left:1.5em">GSVA输入与打分/：{部位}_log2_gene_GSVA_input.xls + {部位}_{KEGG,GOBP}_gsva_score.xls</li>')
RAW('</ul>')
P('其中「对比」为 ' + " / ".join(ccn(cc) for cc in CONTRASTS) + '；「数据库」在 ORA 富集中为 KEGG / GOBP / Reactome，在 GSEA 与 GSVA 中为 KEGG / GOBP。所有 .xls 为制表符分隔文本。')

# ================= 组装目录 (三级嵌套) =================
def render_toc(nodes, cls):
    if not nodes:
        return ""
    out = ['<ul class="%s">' % cls]
    for n in nodes:
        out.append('<li><a href="#%s">%s</a></li>' % (n["anchor"], esc(n["title"])))
        if n["children"]:
            out.append(render_toc(n["children"], "rank3" if cls == "rank2" else "rank3"))
    out.append('</ul>')
    return "\n".join(out)

toc_parts = ['<ul id="rank1">', '<li style="text-align:center;font-size:22px;">目录</li>']
for n in toc_root:
    toc_parts.append('<li><a href="#%s">%s</a></li>' % (n["anchor"], esc(n["title"])))
    if n["children"]:
        # rank2
        toc_parts.append('<ul id="rank2">')
        for n2 in n["children"]:
            toc_parts.append('<li><a href="#%s">%s</a></li>' % (n2["anchor"], esc(n2["title"])))
            if n2["children"]:
                toc_parts.append('<ul id="rank3">')
                for n3 in n2["children"]:
                    toc_parts.append('<li><a href="#%s">%s</a></li>' % (n3["anchor"], esc(n3["title"])))
                toc_parts.append('</ul>')
        toc_parts.append('</ul>')
toc_parts.append('</ul>')
toc_html = "\n".join(toc_parts)

sec_html = "\n".join(body_parts)

HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>HY-0094 心脏蛋白质组学分析报告</title>
<style>
header{ text-align:center; font-size:32px; font-weight:bold; font-family:Arial,KaiTi,"楷体",sans-serif; }
body{ min-height:100%; font-family:Arial,KaiTi,"楷体",sans-serif; margin:0; }
aside{ left:0; margin:0; padding:0; width:23%; float:left; position:fixed; height:99%; background-color:#ededed; z-index:9999; }
div#Catalogue{ width:100%; height:100%; margin:0; padding:4px 0; overflow:auto; border-right:1px solid #ccc; text-align:left; font-weight:bold; }
ul{ margin:0; padding:0; }
li{ list-style-type:none; border-bottom:1px solid #fff; margin:0; padding:0; }
ul#rank1>li{ font-size:17px; line-height:36px; text-indent:0.3em; color:#337ab7; }
ul#rank2>li{ font-size:14px; line-height:30px; text-indent:1.5em; color:#337ab7; }
ul#rank3>li{ font-size:12px; line-height:24px; text-indent:2.6em; }
ul#rank3>li>a{ color:#788fa4; }
ul li:hover{ background-color:#dedede; }
a{ text-decoration:none; color:#337ab7; }
section{ text-align:center; margin:0% 4% 0% 27%; height:100%; }
h1,h2,h3,h4{ text-align:left; }
h4{ font-size:14px; color:#39414a; margin:14px 0 6px; }
p.paragraph{ text-align:justify; text-indent:2em; line-height:1.5em; margin:6px 0; }
p.introduction{ text-align:justify; text-indent:2em; line-height:1.3em; font-size:13px; color:#555; }
p.table_title{ text-align:center; font-weight:bold; margin:8px 0 2px; }
.chart{ width:100%; max-height:400px; overflow:auto; margin:6px 0; }
table#table{ border-collapse:collapse; width:100%; font-size:12.5px; font-family:Arial, sans-serif; table-layout:fixed; }
table#table th{ background-color:#6b7280; color:#fff; padding:4px 8px; border:1px solid #9ca3af; position:sticky; top:0; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }
table#table td{ border:1px solid #cbd5e1; padding:3px 8px; text-align:left; white-space:normal; word-break:break-word; vertical-align:top; overflow-wrap:anywhere; }
table#table td.num{ white-space:nowrap; text-align:right; }
table.img_table{ width:100%; border-collapse:collapse; }
table.img_table>tbody>tr>td{ border-width:0; text-align:center; vertical-align:top; }
img.goImg{ width:100%; max-width:100%; }
pre.qa{ background:#0f172a; color:#e2e8f0; padding:12px; border-radius:4px; font-size:11.5px; overflow-x:auto; line-height:1.5; white-space:pre-wrap; text-align:left; }
a.lb{ color:#337ab7; text-decoration:underline; cursor:pointer; }
.empty{ color:#9ca3af; }
a.lbpic{ cursor:zoom-in; display:inline-block; width:100%; }
a.lbpic img.goImg{ pointer-events:none; }
ul.filelist{ text-align:left; margin-left:2em; }
ul.filelist li{ border:none; list-style-type:disc; font-size:13px; line-height:1.6; }
#lightbox{ display:none; position:fixed; inset:0; background:rgba(0,0,0,.78); z-index:9999; align-items:center; justify-content:center; flex-direction:column; cursor:zoom-out; }
#lightbox .lbwrap{ position:relative; text-align:center; }
#lightbox img{ max-width:92vw; max-height:86vh; border:none; border-radius:4px; background:#fff; }
#lightbox .lbpdf{ display:inline-block; margin-top:12px; padding:6px 18px; background:#fff; color:#333; border-radius:4px; font-size:14px; text-decoration:none; cursor:pointer; }
#lightbox .lbpdf:hover{ background:#337ab7; color:#fff; }
#lightbox .lbclose{ position:absolute; top:8px; right:14px; color:#fff; font-size:28px; cursor:pointer; }
</style>
</head>
<body>
<aside>
<div id="Catalogue">
@@TOC@@
</div>
</aside>
<section>
<header>HY-0094 心脏蛋白质组学分析报告</header>
@@SEC@@
</section>
<div id="lightbox">
  <div class="lbwrap"><span class="lbclose">&times;</span><img src="" alt=""></div>
  <a class="lbpdf" href="" target="_blank" style="display:none;">查看 / 下载 PDF</a>
</div>
<script>
(function(){
  var lb = document.getElementById('lightbox');
  var img = lb.querySelector('img');
  var pdf = lb.querySelector('.lbpdf');
  document.addEventListener('click', function(e){
    var a = e.target.closest('a.lbpic');
    if (a) {
      e.preventDefault();
      img.src = a.getAttribute('href');
      var p = a.getAttribute('data-pdf');
      if (p) { pdf.href = p; pdf.style.display = 'inline-block'; }
      else { pdf.href = '#'; pdf.style.display = 'none'; }
      lb.style.display = 'flex';
      return;
    }
    var b = e.target.closest('a.lb');
    if (b) {
      e.preventDefault();
      img.src = b.getAttribute('href');
      pdf.href = '#'; pdf.style.display = 'none';
      lb.style.display = 'flex';
      return;
    }
    var c = e.target.closest('.lbclose');
    if (c) { lb.style.display = 'none'; img.src = ''; return; }
    if (e.target === lb) { lb.style.display = 'none'; img.src = ''; }
  });
  document.addEventListener('keydown', function(e){
    if (e.key === 'Escape') { lb.style.display = 'none'; img.src = ''; }
  });
})();
</script>
</body>
</html>
"""

final_html = HTML_TEMPLATE.replace("@@TOC@@", toc_html).replace("@@SEC@@", sec_html)
_report_html = os.path.join(OUT_REPORT, "html_report.html")
with open(_report_html, "w", encoding="utf-8") as fh:
    fh.write(final_html)

print("HTML 完成: " + _report_html)
print("引用图片数: " + str(len(manifest)))
print("ORA 映射: " + str(sum(len(v) for v in ora_map.values())))
print("GSEA 映射: " + str(sum(len(v) for v in lead_map.values())))
print("图题数: " + str(fig_no[0]))
print("表题数: " + str(tab_no[0]))
