#!/bin/bash
# Usage: sbatch ... run_root_scan.slurm SETUP_ROOT REPOSITORY INPUT OUTPUT BESTFIT mu...
set -euo pipefail
setup=$1
repository=$2
input=$3
output=$4
bestfit=$5
shift 5
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$output"
exec python3 -u "$repository/root_workspace_fit.py" scan \
    --input "$input" --output "$output" --bestfit "$bestfit" \
    --weighted-asimov --range 0 1.5 --resume --grid "$@"
