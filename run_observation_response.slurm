#!/bin/bash
# SETUP_ROOT SOURCE_REPO INPUT BESTFIT COVARIANCE OUTPUT BUILD_CACHE
set -euo pipefail
source "$1"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$7"
exec python3 -u "$2/compute_observation_response.py" --input "$3" --bestfit "$4" \
    --covariance "$5" --output "$6" --weighted-asimov
