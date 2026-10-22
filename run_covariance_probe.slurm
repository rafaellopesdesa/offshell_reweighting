#!/bin/bash
# SETUP_ROOT SOURCE_REPO INPUT BESTFIT OUTPUT BUILD_CACHE
set -euo pipefail
source "$1"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OFFSHELL_ROOT_BUILD="$6"
python3 -u "$2/compute_observation_response.py" --input "$3" --bestfit "$4" \
    --covariance "$5/covariance.json" --output "$5" --weighted-asimov
exec python3 -u "$2/validate_covariance_curvature.py" --input "$3" --bestfit "$4" \
    --covariance "$5/covariance.json" --output "$5/curvature_validation.json" --weighted-asimov
