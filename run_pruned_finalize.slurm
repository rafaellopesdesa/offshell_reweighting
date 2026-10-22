#!/bin/bash
# SETUP_ROOT SOURCE_REPO OUTPUT
set -euo pipefail
setup=$1
repository=$2
output=$3
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$output/build_cache"
parts=()
for index in $(seq 0 41); do
  parts+=("$output/scan/point_$index/scan_points.json")
done
bash "$repository/run_shifted_finalize.slurm" "$setup" "$repository" \
  "$output/workspace/HZZ_4l_data_like_asimov_4096_pruned54.root" \
  "$output/final" "$output/fit/bestfit.json" "$output/build_cache" "${parts[@]}"
