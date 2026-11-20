#!/bin/bash
# ============================================================================
# 定量蛋白质组学标准分析 pipeline —— 用户配置模板
# 用法: 复制本文件到自己的目录, 改名为任意名(如 my_run.sh), 填好下面的 <...> 项,
#       然后:  CONFIG=/你的路径/my_run.sh  bash /media/nfs/nfs02/chenyuanbao/pipeline/run_pipeline.sh all
# 注意: 本模板的 PP_ROOT 指向【共享只读的 pipeline 安装目录】, 不要改它。
#       你只需要填: PP_SPECIES(物种) / PP_GROUPS(组别) / PP_MATRIX(输入) / PP_OUTPUT(输出)。
# ============================================================================

# ---- 运行环境 (固定 R433 / Python311, 服务器上已装好, 不要改) ----
export PP_RSCRIPT="/media/nfs/nfs02/chenyuanbao/miniforge3/envs/R433/bin/Rscript"
export PP_PYTHON="/media/nfs/nfs02/chenyuanbao/miniforge3/envs/py311/bin/python"

# ---- pipeline 安装目录 (共享只读, 不要改) ----
export PP_ROOT="/media/nfs/nfs02/chenyuanbao/pipeline"
export PP_BIN="$PP_ROOT/bin"
export PP_LIB="$PP_ROOT/lib"

# ---- 用户工作目录 (你自己的输入/输出/日志) ----
export PP_WORKDIR="${PP_WORKDIR:-$HOME/proteomics_analysis}"
export PP_LOGS="$PP_WORKDIR/logs"

# ---- 依赖路径 (pipeline 自带, 不要改) ----
export PP_GMT_BASE="$PP_LIB/databases"
export PP_FONT_DIR="$PP_LIB/fonts/Arial"
export PP_IDMAP="$PP_LIB/id_map.tsv"
export PP_HEATMAP_SCRIPT="$PP_LIB/heatmap/heatmap_light.r"

# ============================================================================
# 【1】物种: human / mouse / rat  (三选一, 一键切换基因集)
# ============================================================================
export PP_SPECIES="<填 human 或 mouse 或 rat>"
case "$PP_SPECIES" in
  human)
    export PP_GMT_KEGG="$PP_GMT_BASE/GMT/hsa.kegg_optimal.symbol.gmt"
    export PP_GMT_GOBP="$PP_GMT_BASE/GMT/c5.go.bp.v7.5.1.symbols.gmt"
    export PP_GMT_REACTOME="$PP_GMT_BASE/GMT/c2.cp.reactome.v7.5.1.symbols.gmt"
    ;;
  mouse)
    export PP_GMT_KEGG="$PP_GMT_BASE/GMT_mm/mmu_kegg.symbol.gmt"
    export PP_GMT_GOBP="$PP_GMT_BASE/GMT_mm/mmu_go_bp.symbol.gmt"
    export PP_GMT_REACTOME="$PP_GMT_BASE/GMT_mm/mmu_reactome.symbol.gmt"
    ;;
  rat)
    export PP_GMT_KEGG="$PP_GMT_BASE/GMT_rn/rno_kegg.symbol.gmt"
    export PP_GMT_GOBP="$PP_GMT_BASE/GMT_rn/rno_go_bp.symbol.gmt"
    export PP_GMT_REACTOME="$PP_GMT_BASE/GMT_rn/rno_reactome.symbol.gmt"
    ;;
  *) echo "PP_SPECIES 必须是 human/mouse/rat" >&2; exit 1 ;;
esac

# ============================================================================
# 【2】输入数据  (填你的文件路径, 绝对路径最稳)
# ============================================================================
export PP_INPUT="<填输入目录, 如 /home/你的用户名/我的数据>"
export PP_MATRIX="$PP_INPUT/report.pg_matrix.tsv"   # DIA-NN 蛋白组矩阵
# export PP_MATRIX_XLSX="$PP_INPUT/我的数据.xlsx"    # 若输入是 xlsx 基因表, 取消注释并填路径

# ============================================================================
# 【3】实验设计 —— 组别信息在这里自由填写 (任意组名/任意组数)
# ============================================================================
export PP_REGIONS="<填部位/批次, 如 IVS 或 IVS LVAW LVPW>"
export PP_GROUPS="<填组别, 空格分隔, 对照在前 实验在后, 如 Normal Model>"
export PP_REPLICATES="6"
export PP_CONTRASTS="<填对比, 如 Model_vs_Normal>"
export PP_GROUP_CONTROL="<填对照组名, 需在 PP_GROUPS 内>"
export PP_GROUP_TREAT="<填实验组名, 需在 PP_GROUPS 内>"
# export PP_COLOR_对照组名="#008B45"   # 可选: 每组颜色 PP_COLOR_<组名>, 缺省自动配色

# ============================================================================
# 【4】输出目录  (填你想存结果的地方, 绝对路径)
# ============================================================================
export PP_OUTPUT="<填输出目录, 如 /home/你的用户名/我的结果>"
export PP_OUT_RAW="$PP_OUTPUT/01_raw"
export PP_OUT_STD="$PP_OUTPUT/02_standardize"
export PP_OUT_QC="$PP_OUTPUT/03_qc_pca_cluster"
export PP_OUT_DE="$PP_OUTPUT/04_de"
export PP_OUT_ORA="$PP_OUTPUT/05_ora"
export PP_OUT_GSEA="$PP_OUTPUT/06_gsea"
export PP_OUT_GSVA="$PP_OUTPUT/07_gsva"
export PP_OUT_HEATMAP="$PP_OUTPUT/08_pathway_heatmaps"
export PP_OUT_REPORT="$PP_OUTPUT/09_report"

# ---- 分析阈值 (一般不用改) ----
export PP_PCUT=0.05
export PP_FCCUT="0.2630344058337938"
export PP_ORA_TOPN=15
export PP_ORA_HEATMAP_TOPN=30
export PP_GSEA_BAR_TOPN=10
export PP_GSEA_LEAD_TOPN=30
export PP_GSVA_TOPN=15
export PP_THREADS=8
