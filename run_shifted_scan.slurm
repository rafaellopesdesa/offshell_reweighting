#!/bin/bash
# SETUP_ROOT SOURCE_REPO INPUT OUTPUT BESTFIT BUILD_CACHE mu...
set -euo pipefail
setup=$1
repository=$2
input=$3
output=$4
bestfit=$5
build_cache=$6
shift 6
source "$setup"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$build_cache"
exec python3 -u "$repository/root_workspace_fit.py" scan --input "$input" \
  --output "$output" --bestfit "$bestfit" --weighted-asimov \
  --range 0 2 --allow-zero --resume --grid "$@"
