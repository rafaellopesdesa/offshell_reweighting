#!/bin/bash
# SETUP_ROOT SOURCE_REPO INPUT BESTFIT OUTPUT BUILD_CACHE
set -euo pipefail
source "$1"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$6"
exec python3 -u "$2/compute_fit_covariance.py" --input "$3" --bestfit "$4" \
    --output "$5" --weighted-asimov
