#!/bin/bash
# SETUP_ROOT SOURCE_REPO OUTPUT
# One independent process and Slurm array task per coordinate gives attributable RSS.
set -euo pipefail
setup=$1
repository=$2
output=$3
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$output/build_cache"
coordinate=$(python3 -c 'import json,os,sys; grid=json.load(open(sys.argv[1])); print(grid[int(os.environ["SLURM_ARRAY_TASK_ID"])])' "$output/grid.json")
exec python3 -u "$repository/root_workspace_fit.py" scan \
  --input "$output/workspace/HZZ_4l_data_like_asimov_4096_pruned54.root" \
  --output "$output/scan/point_${SLURM_ARRAY_TASK_ID}" --bestfit "$output/fit/bestfit.json" \
  --weighted-asimov --range 0 2 --allow-zero --grid "$coordinate"
