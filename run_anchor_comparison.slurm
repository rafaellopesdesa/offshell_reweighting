#!/bin/bash
set -euo pipefail
setup=$1
repository=$2
results=$3
reference=$4
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec python3 -u "$repository/compare_asimov_hepdata.py" \
    --fits "$results" --sample-subdirectory fit --normalization-mode anchors \
    --reference "$reference/table_2.json" --crosscheck-table "$reference/table_1.json" \
    --reference-csv "$reference/table_2.csv" --output "$results/comparison"
