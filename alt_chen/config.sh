#!/bin/bash
# ============================================================================
# 定量蛋白质组学标准分析 pipeline 配置
# 修改这里的变量即可适配新的数据集，无需改动脚本
# ============================================================================

# ---- 运行环境 (固定 R433) ----
export PP_RSCRIPT="/media/nfs/nfs02/chenyuanbao/miniforge3/envs/R433/bin/Rscript"
export PP_PYTHON="/media/nfs/nfs02/chenyuanbao/miniforge3/envs/py311/bin/python"

# ---- pipeline 根目录 (共享只读, 自动定位) ----
export PP_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # 自动指向 pipeline 目录所在位置
export PP_BIN="$PP_ROOT/bin"
export PP_LIB="$PP_ROOT/lib"

# ---- 用户工作目录 (每位用户自己的输入/输出/日志; 默认在家目录) ----
export PP_WORKDIR="${PP_WORKDIR:-$HOME/proteomics_analysis}"
export PP_LOGS="$PP_WORKDIR/logs"

# ---- 依赖路径 ----
export PP_GMT_BASE="$PP_LIB/databases"          # 三个物种 GMT 的总目录
export PP_FONT_DIR="$PP_LIB/fonts/Arial"
export PP_IDMAP="$PP_LIB/id_map.tsv"
export PP_HEATMAP_SCRIPT="$PP_LIB/heatmap/heatmap_light.r"

# ---- 物种 (human / mouse / rat) ----
export PP_SPECIES="human"
case "$PP_SPECIES" in
  human)
    export PP_ORGDB="org.Hs.eg.db"
    export PP_GMT_KEGG="$PP_GMT_BASE/GMT/hsa.kegg_optimal.symbol.gmt"
    export PP_GMT_GOBP="$PP_GMT_BASE/GMT/c5.go.bp.v7.5.1.symbols.gmt"
    export PP_GMT_REACTOME="$PP_GMT_BASE/GMT/c2.cp.reactome.v7.5.1.symbols.gmt"
    ;;
  mouse)
    export PP_ORGDB="org.Mm.eg.db"
    export PP_GMT_KEGG="$PP_GMT_BASE/GMT_mm/mmu_kegg.symbol.gmt"
    export PP_GMT_GOBP="$PP_GMT_BASE/GMT_mm/mmu_go_bp.symbol.gmt"
    export PP_GMT_REACTOME="$PP_GMT_BASE/GMT_mm/mmu_reactome.symbol.gmt"
    ;;
  rat)
    export PP_ORGDB="org.Rn.eg.db"
    export PP_GMT_KEGG="$PP_GMT_BASE/GMT_rn/rno_kegg.symbol.gmt"
    export PP_GMT_GOBP="$PP_GMT_BASE/GMT_rn/rno_go_bp.symbol.gmt"
    export PP_GMT_REACTOME="$PP_GMT_BASE/GMT_rn/rno_reactome.symbol.gmt"
    ;;
  *) echo "PP_SPECIES 必须是 human/mouse/rat, 当前=$PP_SPECIES" >&2; exit 1 ;;
esac

# ---- 输入数据 (PP_INPUT 可用环境变量覆盖; 默认在 PP_WORKDIR 下) ----
export PP_INPUT="${PP_INPUT:-$PP_WORKDIR/input}"
export PP_MATRIX="$PP_INPUT/report.pg_matrix.tsv"   # DIA-NN 蛋白组矩阵
# 可选: 样本登记表 CSV (含 region_cn/group_cn/animal_id 三列, 行序=样本编号顺序; 留空则跳过核验)
export PP_REG_CSV=""

# ============================================================================
# 实验设计 (组别名/组数均可自由填写, 任意命名)
# ============================================================================
export PP_REGIONS="IVS"                       # 部位/批次列表(空格分隔)
export PP_GROUPS="Normal Model"               # 组别列表(空格分隔, 顺序=样本排列顺序; 对照在前, 实验在后)
export PP_REPLICATES=6                        # 每组生物学重复数
export PP_CONTRASTS="Model_vs_Normal"         # 对比列表(空格分隔; A_vs_B = A - B)

# GSVA 方向: 实验组比对照组 (effect = 实验组均值 - 对照组均值)
export PP_GROUP_CONTROL="Normal"              # 对照组名(须是 PP_GROUPS 之一)
export PP_GROUP_TREAT="Model"                 # 实验组名(须是 PP_GROUPS 之一)

# ---- 组别颜色 (按组名 PP_COLOR_<组名>; 缺省则回退调色板) ----
export PP_COLOR_NORMAL="#008B45"
export PP_COLOR_MODEL="#EE0000"

# ---- 分析阈值 ----
export PP_PCUT=0.05                # pvalue 显著性阈值
export PP_FCCUT="0.2630344058337938"  # |log2FC| 阈值 = log2(1.2)
export PP_ORA_TOPN=15              # ORA 气泡图 top N
export PP_ORA_HEATMAP_TOPN=30      # 通路蛋白热图 top N
export PP_GSEA_BAR_TOPN=10         # GSEA 柱状图 top N
export PP_GSEA_LEAD_TOPN=30        # GSEA 领头热图 top N
export PP_GSVA_TOPN=15             # GSVA 两方向 top N

# ---- 输出目录 (按分析步骤组织; PP_OUTPUT 可用环境变量覆盖) ----
export PP_OUTPUT="${PP_OUTPUT:-$PP_WORKDIR/output}"
export PP_OUT_RAW="$PP_OUTPUT/01_raw"
export PP_OUT_STD="$PP_OUTPUT/02_standardize"
export PP_OUT_QC="$PP_OUTPUT/03_qc_pca_cluster"
export PP_OUT_DE="$PP_OUTPUT/04_de"
export PP_OUT_ORA="$PP_OUTPUT/05_ora"
export PP_OUT_GSEA="$PP_OUTPUT/06_gsea"
export PP_OUT_GSVA="$PP_OUTPUT/07_gsva"
export PP_OUT_HEATMAP="$PP_OUTPUT/08_pathway_heatmaps"
export PP_OUT_REPORT="$PP_OUTPUT/09_report"

export PP_THREADS=8
