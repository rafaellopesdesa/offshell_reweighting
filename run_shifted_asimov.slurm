#!/bin/bash
# SETUP_ROOT SOURCE_REPO TRUE_WORKSPACE OBSERVED_BESTFIT OUTPUT
set -euo pipefail
setup=$1
repository=$2
original=$3
observed=$4
output=$5
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$output/build_cache"
mkdir -p "$output/fit"
python3 -u "$repository/build_shifted_asimov.py" --source "$original" \
  --observed-bestfit "$observed" --output "$output/workspace"
input="$output/workspace/HZZ_4l_data_like_asimov_4096.root"
python3 -u "$repository/validate_shifted_asimov.py" --input "$input" \
  --output "$output/workspace/validation.json"
python3 -u "$repository/root_workspace_fit.py" fit --input "$input" --output "$output/fit" \
  --weighted-asimov --starts 1 1.4 --fit-strategy 0 --fit-eps 0.0001
python3 -u "$repository/validate_shifted_asimov.py" --input "$input" \
  --bestfit "$output/fit/bestfit.json" --output "$output/fit/global_validation.json"
