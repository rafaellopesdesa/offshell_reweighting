#!/bin/bash
# Submit with explicit, currently valid Unity account/partition/resources.
# Usage: sbatch ... run_root_fit.slurm SETUP_ROOT REPOSITORY INPUT OUTPUT [fit options]
set -euo pipefail
setup=$1
repository=$2
input=$3
output=$4
shift 4
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$output"
exec python3 -u "$repository/root_workspace_fit.py" fit \
    --input "$input" --output "$output" --weighted-asimov --hesse "$@"
