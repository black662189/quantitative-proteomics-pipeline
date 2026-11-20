#!/usr/bin/env bash
set -euo pipefail
PIPELINE="/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/proteomics"
CONFIG="${1:-config.yaml}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec "$PIPELINE/.venv/bin/python" "$PIPELINE/proteomics.py" run -c "$CONFIG"
