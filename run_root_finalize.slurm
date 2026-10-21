#!/bin/bash
# Usage: sbatch ... run_root_finalize.slurm SETUP_ROOT REPOSITORY INPUT OUTPUT mu...
set -euo pipefail
setup=$1
repository=$2
input=$3
output=$4
shift 4
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
parts=("$output"/scan_*/scan_points.json)
python3 -u "$repository/root_workspace_fit.py" finalize \
    --input "$input" --output "$output" --weighted-asimov \
    --range 0 1.5 --grid "$@" --scan-parts "${parts[@]}"
exec python3 -u "$repository/root_workspace_validation.py" \
    --input "$input" --weighted-asimov --range 0.01 1.5 \
    --bestfit "$output/bestfit.json" --scan-points "$output/scan_points.json" \
    --output "$output/likelihood_validation.json"
