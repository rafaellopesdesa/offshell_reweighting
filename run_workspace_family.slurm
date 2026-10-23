#!/bin/bash
# SETUP_ROOT FROZEN_SOURCE RUN_DIRECTORY STAGE MODEL_INDEX
set -euo pipefail
setup=$1
repository=$2
run_directory=$3
stage=$4
model_index=$5
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$run_directory/build_cache"
exec python3 -u "$repository/workspace_family_scan.py" "$stage" \
  --run "$run_directory" --model-index "$model_index"
