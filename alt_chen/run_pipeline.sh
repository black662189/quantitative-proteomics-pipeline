#!/bin/bash
# ============================================================================
# 定量蛋白质组学标准分析 pipeline 总控脚本
# 用法: bash run_pipeline.sh [step1|step2|...|all]
# 按分析步骤组织输出 (见 config.sh 的 PP_OUT_*)
# ============================================================================
set -euo pipefail
# 定位 run_pipeline.sh 所在目录(共享 pipeline 目录), 用于默认 config
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# CONFIG 可为绝对路径或相对路径(相对当前工作目录); 默认用 pipeline 目录下的 config.sh
CFG="${CONFIG:-$SELF_DIR/config.sh}"
[ -f "$CFG" ] || { echo "配置文件不存在: $CFG"; exit 1; }
source "$CFG"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
log "使用配置: $CFG"

mkdir -p "$PP_LOGS" "$PP_OUTPUT" \
  "$PP_OUT_RAW" "$PP_OUT_STD" "$PP_OUT_QC" "$PP_OUT_DE" \
  "$PP_OUT_ORA" "$PP_OUT_GSEA" "$PP_OUT_GSVA" "$PP_OUT_HEATMAP" "$PP_OUT_REPORT"
R="$PP_RSCRIPT"

step0_convert() {
  # 可选: 将非标准输入 (xlsx) 转成统一 TSV
  # PP_INPUT_FORMAT: pd_tmt(PD蛋白级TMT总量表) | ima94(基因水平定量表) | 空/其他(跳过)
  if [ -n "${PP_MATRIX_XLSX:-}" ] && [ -f "$PP_MATRIX_XLSX" ]; then
    local fmt="${PP_INPUT_FORMAT:-}"
    case "$fmt" in
      pd_tmt)
        log "步骤0: PD TMT 总量表转换 ($PP_MATRIX_XLSX -> $PP_MATRIX)"
        "$PP_PYTHON" "$PP_BIN/00_convert_pd_tmt.py"
        ;;
      ima94)
        log "步骤0: ima94 基因表转换 ($PP_MATRIX_XLSX -> $PP_MATRIX)"
        "$PP_PYTHON" "$PP_BIN/00_convert_ima94.py"
        ;;
      *)
        log "步骤0: 有 xlsx 但未设置 PP_INPUT_FORMAT, 跳过 (请设 pd_tmt 或 ima94)"
        ;;
    esac
  else
    log "步骤0: 无 xlsx 输入, 跳过"
  fi
}

step1_standardize() {
  log "步骤1: 数据标准化 (原始矩阵 -> log2/填充矩阵)"
  mkdir -p "$PP_OUT_RAW"
  cp -f "$PP_MATRIX" "$PP_OUT_RAW/" 2>/dev/null || true
  "$R" --vanilla "$PP_BIN/01_process_standard.R"
}

step2_qc() {
  log "步骤2: PCA / t-SNE / 层次聚类 QC 出图"
  for rg in $PP_REGIONS; do
    local d="$PP_OUT_QC/$rg"
    mkdir -p "$d"
    ( cd "$d" && "$R" --vanilla "$PP_BIN/02_cluster_pca.R" \
        "$PP_OUT_STD/${rg}_protein_expression_log2_imputed.xls" $PP_GROUPS )
    # 报告展示需要的 png (ImageMagick 转换)
    local grpstr=$(echo $PP_GROUPS | sed 's/ /-vs-/g')
    convert -density 200 "$d/PCA.$grpstr.pdf" "$d/PCA.$grpstr.png" 2>/dev/null || true
    convert -density 200 "$d/hclust_average.$grpstr.pdf" "$d/hclust_average.$grpstr.png" 2>/dev/null || true
    log "  $rg QC 完成"
  done
}

step3_de() {
  log "步骤3: limma 差异分析 + 火山图"
  for rg in $PP_REGIONS; do
    "$R" --vanilla "$PP_BIN/03_DE_limma.R" "$rg"
    "$R" --vanilla "$PP_BIN/03b_volcano.R" "$rg"
    log "  $rg 差异完成"
  done
}

step4_ora() {
  log "步骤4: ORA 富集分析 (KEGG/GOBP/Reactome)"
  for rg in $PP_REGIONS; do
    "$R" --vanilla "$PP_BIN/04_ORA_enrichment.R" "$rg"
    log "  $rg ORA 完成"
  done
}

step5_gsea() {
  log "步骤5: GSEA (KEGG/GOBP)"
  for rg in $PP_REGIONS; do
    "$R" --vanilla "$PP_BIN/05_full_GSEA.R" "$rg"
    "$R" --vanilla "$PP_BIN/06_GSEA_plots.R" "$rg"
    log "  $rg GSEA 完成"
  done
}

step6_gsva() {
  log "步骤6: GSVA (KEGG/GOBP)"
  for rg in $PP_REGIONS; do
    "$R" --vanilla "$PP_BIN/07a_prepare_GSVA_input.R" "$rg"
    for db in GOBP KEGG; do
      local gmt="$PP_GMT_GOBP"; [ "$db" = "KEGG" ] && gmt="$PP_GMT_KEGG"
      "$R" --vanilla "$PP_BIN/07b_GSVA_calc.R" "$gmt" \
        "$PP_OUT_GSVA/${rg}_log2_gene_GSVA_input.xls" \
        "$PP_OUT_GSVA/${rg}_${db}_gsva_score.xls" "$PP_THREADS"
    done
    "$R" --vanilla "$PP_BIN/07c_GSVA_top_heatmap.R" "$rg"
    log "  $rg GSVA 完成"
  done
}

step7_heatmap() {
  log "步骤7: 通路蛋白热图 + GSEA 前导热图"
  for rg in $PP_REGIONS; do
    "$R" --vanilla "$PP_BIN/08a_pathway_protein_heatmaps.R" "$rg"
    "$R" --vanilla "$PP_BIN/08b_GSEA_leading_heatmaps.R" "$rg"
    log "  $rg 热图完成"
  done
}

step8_report() {
  log "步骤8: 生成网页报告"
  "$PP_PYTHON" "$PP_BIN/make_report.py"
  log "  报告完成: $PP_OUT_REPORT/html_report.html"
}

run_all() {
  step0_convert
  step1_standardize
  step2_qc
  step3_de
  step4_ora
  step5_gsea
  step6_gsva
  step7_heatmap
  step8_report
  log "全流程完成"
}

case "${1:-all}" in
  step0) step0_convert ;;
  step1) step1_standardize ;;
  step2) step2_qc ;;
  step3) step3_de ;;
  step4) step4_ora ;;
  step5) step5_gsea ;;
  step6) step6_gsva ;;
  step7) step7_heatmap ;;
  step8) step8_report ;;
  all)   run_all ;;
  *) echo "用法: bash run_pipeline.sh {step0|step1|step2|step3|step4|step5|step6|step7|step8|all}"; exit 1 ;;
esac
