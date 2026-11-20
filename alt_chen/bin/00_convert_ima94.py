#!/usr/bin/env python3
# ============================================================================
# IMA94 适配: 蛋白总定量表.xlsx -> pipeline 期望的"蛋白水平"TSV
# 输入: 基因水平定量表 (列: gene_id + 6样本 + Fold.Change + P.Value)
# 输出: 4 注释列 (Protein.Group/Names/Genes/First.Description) + 定量列(按对照在前)
# 规则: 丢弃全空行; Fold.Change/P.Value 列忽略(由 pipeline 重新计算);
#       region 名 = PP_REGION_TAG (默认 IMA), 样本名 = {group}_{region}-{rep}
# ============================================================================
import openpyxl, math, os, sys

XLSX   = os.environ["PP_MATRIX_XLSX"]
OUTTSV = os.environ["PP_MATRIX"]
REGION = os.environ.get("PP_REGION_TAG", "IMA")
GROUPS = os.environ["PP_GROUPS"].split()          # 对照在前, 实验在后
REP    = int(os.environ["PP_REPLICATES"])

# 原始 xlsx 样本列名 -> (组, 重复序号)
# 判定组: 按 config GROUPS 顺序; 样本名前缀匹配组名
def norm(x):
    return x.replace("-", "").upper()   # 去连字符统一比较 (HY-0094 -> HY0094)

G_NORM = [norm(g) for g in GROUPS]
def map_sample(col):
    c0 = col.split("-")[0]              # "HY-0094-3" -> "HY"
    full = norm(col.split("-")[0] + "-".join(col.split("-")[:2]) if col.count("-") >= 2 else col)
    # 直接: 取 col 去掉末尾 -rep 后的组名部分
    parts = col.split("-")
    base = "-".join(parts[:-1]) if len(parts) > 1 else col   # "HY-0094-3" -> "HY-0094"
    nb = norm(base)
    for i, gn in enumerate(G_NORM):
        if nb == gn:
            return GROUPS[i], int(parts[-1])
    raise SystemExit(f"无法识别样本列 {col} 属于哪个组 (GROUPS={GROUPS}, 规范化={nb})")

wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb["Sheet 1"] if "Sheet 1" in wb.sheetnames else wb.worksheets[0]
rows = list(ws.iter_rows(values_only=True))
hdr = list(rows[0])
data = rows[1:]

# 定位样本列 (排除 gene_id / Fold.Change / P.Value)
sample_cols = [(i, c) for i, c in enumerate(hdr)
               if c not in ("gene_id", "Fold.Change", "P.Value")]

def isnum(v):
    return v is not None and not (isinstance(v, float) and math.isnan(v))

# 每行保留: 6 样本全有值
def all_present(r):
    return all(isnum(r[i]) for i, _c in sample_cols)

kept = [r for r in data if r[0] and all_present(r)]
print(f"IMA94: {len(data)} 总行 -> 保留有值行 {len(kept)}", file=sys.stderr)

# 目标样本顺序: GROUPS 顺序(对照在前), 组内按原始列名末尾序号升序
# (xlsx 中实验组重复编号可能是 3/4/5 而非 1/2/3, 故按实际列名末尾数字排序后重新编号 1..REP)
order = []  # (组, 新rep, 原始列idx)
for g in GROUPS:
    # sample_cols 元素为 (idx, 列名)
    grp_cols = [(i, c) for i, c in sample_cols if map_sample(c)[0] == g]
    grp_cols_sorted = sorted(grp_cols, key=lambda ic: int(ic[1].split("-")[-1]))
    if len(grp_cols_sorted) != REP:
        raise SystemExit(f"组 {g} 样本列数 {len(grp_cols_sorted)} != REP {REP}")
    for newrep, (_i, c) in enumerate(grp_cols_sorted, start=1):
        order.append((g, newrep, _i))

# 6 注释列, 与 step1 的 -(1:6) 假设对齐 (N.Sequences/N.Proteotypic.Sequences 置 1)
newhdr = ["Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
          "N.Sequences", "N.Proteotypic.Sequences"]
newhdr += [f"{g}_{REGION}-{rep}" for g, rep, _ in order]

out = [newhdr]
for r in kept:
    gid = str(r[0]).strip()
    # 0 值视为"未检出" -> 输出为空(NaN), 交由 step1 的缺失值填充逻辑处理
    def fmt(v):
        if v is None or (isinstance(v, float) and math.isnan(v)) or v == 0:
            return ""
        return str(v)
    vals = [fmt(r[i]) for _g, _rep, i in order]
    row = [gid, gid, gid, "", "1", "1"] + vals
    out.append(row)

with open(OUTTSV, "w", encoding="utf-8") as f:
    for row in out:
        f.write("\t".join(row) + "\n")

print(f"写出 {len(out)-1} 基因 x {len(order)} 样本 -> {OUTTSV}", file=sys.stderr)
