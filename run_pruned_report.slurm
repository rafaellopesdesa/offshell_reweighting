#!/bin/bash
# PYTHON SOURCE_REPO RUN RESULTS_ROOT
set -euo pipefail
interpreter=$1
repository=$2
run=$3
results=$4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="$run/matplotlib_cache"
exec "$interpreter" -u "$repository/report_pruned_asimov.py" \
  --run "$run" --results "$results" --output "$run/report"
