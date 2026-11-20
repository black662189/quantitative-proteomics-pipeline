#!/usr/bin/env python3
# ============================================================================
# PD 蛋白级 TMT 定量表 (xlsx) -> chenyuanbao pipeline 期望的 DIA-NN 风格 TSV
#
# 输入:  PD 导出的"蛋白总量表" xlsx (sheet 含 Accession/Description/.../Gene Symbol/
#        以及 N 个 "*_Abundances (Grouped): ..." 定量列)
# 输出:  前 6 列注释 (Protein.Group/Protein.Names/Genes/First.Protein.Description/
#        N.Sequences/N.Proteotypic.Sequences) + 定量列, 样本名 = {group}_{region}-{rep}
#        列顺序 = PP_GROUPS 顺序 (对照在前, 实验在后)
#
# 环境变量:
#   PP_MATRIX_XLSX   输入 xlsx 路径
#   PP_MATRIX        输出 TSV 路径 (= step1 的 PP_MATRIX)
#   PP_GROUPS        组别 (空格分隔, 对照在前)
#   PP_REPLICATES    每组重复数
#   PP_REGION_TAG    部位标签 (默认 IVS)
#   PP_SHEET         页签名 (默认第一个 sheet)
#   PP_ID_COL        蛋白 ID 列名 (默认 Accession)
#   PP_GENE_COL      基因符号列名 (默认 Gene Symbol)
#   PP_DESC_COL      描述列名 (默认 Description)
#   PP_ABUND_PATTERN 定量列名含有的子串 (默认 "Abundances")
# ============================================================================
import openpyxl, math, os, sys

XLSX   = os.environ["PP_MATRIX_XLSX"]
os.makedirs(os.path.dirname(os.path.abspath(os.environ["PP_MATRIX"])), exist_ok=True)
OUTTSV = os.environ["PP_MATRIX"]
GROUPS = os.environ["PP_GROUPS"].split()
REP    = int(os.environ["PP_REPLICATES"])
REGION = os.environ.get("PP_REGION_TAG", "IVS")
SHEET  = os.environ.get("PP_SHEET", "")
ID_COL = os.environ.get("PP_ID_COL", "Accession")
GENE_COL= os.environ.get("PP_GENE_COL", "Gene Symbol")
DESC_COL= os.environ.get("PP_DESC_COL", "Description")
ABUND  = os.environ.get("PP_ABUND_PATTERN", "Abundances")

wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb[SHEET] if SHEET and SHEET in wb.sheetnames else wb.worksheets[0]
rows = list(ws.iter_rows(values_only=True))
hdr = list(rows[0])
data = rows[1:]

hi = {name: i for i, name in enumerate(hdr)}
for need in (ID_COL, GENE_COL, DESC_COL):
    if need not in hi:
        raise SystemExit(f"缺少关键列: {need} (表头={hdr[:8]}...)")

# 定量列: 含指定子串且非 Fold.Change/P.Value 类统计列
abund_idx = [i for i, h in enumerate(hdr)
             if isinstance(h, str) and ABUND in h
             and not h.startswith(("Fold.Change", "P.Value", "q-value"))]

def isnum(v):
    return v is not None and not (isinstance(v, float) and math.isnan(v))

# 定量列 -> 组别/重复序号: 依赖列名顺序 = GROUPS 顺序 x 组内重复顺序 (PD 常规导出顺序)
# 若无法从列名推断, 按 PP_GROUPS 顺序均分 (对照前 REP 列, 实验后 REP 列)
if len(abund_idx) != REP * len(GROUPS):
    raise SystemExit(f"定量列数 {len(abund_idx)} != {REP}x{len(GROUPS)} (检查 PP_ABUND_PATTERN)")

order = []  # (group, rep, orig_col_idx)
for gi, g in enumerate(GROUPS):
    seg = abund_idx[gi * REP:(gi + 1) * REP]
    for rep, ci in enumerate(seg, start=1):
        order.append((g, rep, ci))

# 输出头: 前 6 注释 + 定量列
newhdr = ["Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
          "N.Sequences", "N.Proteotypic.Sequences"]
newhdr += [f"{g}_{REGION}-{rep}" for g, rep, _ in order]

n_out = 0
with open(OUTTSV, "w", encoding="utf-8") as f:
    f.write("\t".join(newhdr) + "\n")
    for r in data:
        pid = r[hi[ID_COL]]
        if pid is None or str(pid).strip() == "":
            continue
        pid = str(pid).strip()
        gene = r[hi[GENE_COL]]
        gene = "" if gene is None else str(gene).strip()
        desc = r[hi[DESC_COL]]
        desc = "" if desc is None else str(desc).replace("\t", " ").strip()
        # 数量列 (不可用时写空, 与 DIA-NN 一致为可选内容)
        ns = 1
        np = 1
        if "# Unique Peptides" in hi:
            v = r[hi["# Unique Peptides"]]
            if isnum(v):
                ns = str(v)
        def fmt(v):
            if not isnum(v) or v == 0:
                return ""
            return str(v)
        vals = [fmt(r[ci]) for _g, _rep, ci in order]
        f.write("\t".join([pid, pid, gene, desc, str(ns), str(ns)] + vals) + "\n")
        n_out += 1

print(f"PD TMT 转换: {len(data)} 行 -> 写出 {n_out} 蛋白 x {len(order)} 样本 -> {OUTTSV}", file=sys.stderr)
