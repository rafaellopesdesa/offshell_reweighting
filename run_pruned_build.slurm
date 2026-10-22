#!/bin/bash
# SETUP_ROOT SOURCE_REPO INPUT RANKING OUTPUT
set -euo pipefail
setup=$1
repository=$2
input=$3
ranking=$4
output=$5
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$output/build_cache"
python3 -u "$repository/prune_data_like_asimov.py" --input "$input" \
  --ranking "$ranking" --keep 54 --output "$output/workspace"
python3 -u "$repository/validate_pruned_asimov.py" \
  --input "$output/workspace/HZZ_4l_data_like_asimov_4096_pruned54.root" \
  --output "$output/workspace/validation.json"
