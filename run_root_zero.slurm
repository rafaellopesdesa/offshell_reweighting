#!/bin/bash
# Usage: sbatch ... run_root_zero.slurm SETUP_ROOT REPOSITORY INPUT OUTPUT BESTFIT
set -euo pipefail
setup=$1
repository=$2
input=$3
output=$4
bestfit=$5
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$output"
python3 -u "$repository/root_workspace_fit.py" scan \
    --input "$input" --output "$output" --bestfit "$bestfit" \
    --weighted-asimov --range 0 1.5 --allow-zero --grid 0
exec python3 -u "$repository/root_workspace_validation.py" \
    --input "$input" --weighted-asimov --range 0 1.5 \
    --bestfit "$bestfit" --scan-points "$output/scan_points.json" \
    --output "$output/likelihood_validation.json"
