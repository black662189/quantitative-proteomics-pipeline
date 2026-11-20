#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$ROOT/env" "$ROOT/validation"
"$ROOT/.venv/bin/pip" freeze > "$ROOT/env/requirements.lock.txt"
"$ROOT/.venv/bin/python" -m unittest discover -s "$ROOT/tests" -v > "$ROOT/validation/tests.log" 2>&1
"$ROOT/.venv/bin/python" "$ROOT/proteomics.py" run -c "$ROOT/examples/TMT_HepG2/config.yaml" > "$ROOT/validation/TMT_HepG2.run.log" 2>&1
"$ROOT/.venv/bin/python" "$ROOT/proteomics.py" run -c "$ROOT/examples/DIA_cat/config.yaml" > "$ROOT/validation/DIA_cat.run.log" 2>&1
echo "BOTH_EXAMPLES_COMPLETE"
