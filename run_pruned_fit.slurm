#!/bin/bash
# SETUP_ROOT SOURCE_REPO OUTPUT
set -euo pipefail
setup=$1
repository=$2
output=$3
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$output/build_cache"
input="$output/workspace/HZZ_4l_data_like_asimov_4096_pruned54.root"
python3 -u "$repository/root_workspace_fit.py" fit --input "$input" \
  --output "$output/fit" --weighted-asimov --starts 1.4 --fit-strategy 0 --fit-eps 0.0001
python3 -u "$repository/validate_pruned_asimov.py" --input "$input" \
  --bestfit "$output/fit/bestfit.json" --output "$output/fit/global_validation.json"
