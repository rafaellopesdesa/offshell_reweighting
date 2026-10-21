#!/bin/bash
# Submit one job per sample with explicit Unity account, partition, and resources.
# Arguments: ROOT_SETUP REPO STAGED TEMPLATE PREVIOUS_WORKSPACE OUTPUT SIZE
set -euo pipefail
setup=$1
repository=$2
staged=$3
template=$4
previous=$5
output=$6
size=$7
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
workspace="$output/workspace/HZZ_4l_asimov_$size.root"
python3 -u "$repository/build_root_asimov.py" --template "$template" \
  --staged "$staged" --size "$size" --output "$output/workspace" --normalize-systematic-anchors
python3 -u "$repository/validate_root_asimov.py" --input "$workspace" \
  --staged "$staged" --output "$output/workspace/validation.json"
python3 -u "$repository/root_workspace_fit.py" fit --input "$workspace" \
  --output "$output/fit" --weighted-asimov --fit-strategy 0 --fit-eps 0.0001 --starts 1 1.4
python3 -u "$repository/validate_anchor_global_fit.py" --input "$workspace" \
  --bestfit "$output/fit/bestfit.json" --previous-workspace "$previous" \
  --output "$output/fit/global_validation.json"
