#!/bin/bash
# SETUP_ROOT REPOSITORY INPUT OUTPUT BESTFIT BUILD_CACHE scan_parts...
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
mapfile -t grid < <(python3 -c 'print("\n".join(map(str,[0.,.01]+[round(i*.05,2) for i in range(1,41)])))')
python3 -u "$repository/root_workspace_fit.py" finalize --input "$input" \
  --output "$output" --bestfit "$bestfit" --weighted-asimov --range 0 2 \
  --allow-zero --grid "${grid[@]}" --scan-parts "$@"
covariance_options=()
if python3 -c 'import json,sys; sys.exit(0 if json.load(open(sys.argv[1]))["covariance_quality"] != 3 else 1)' "$bestfit"; then
  covariance_options+=(--allow-approximate-covariance)
fi
python3 -u "$repository/root_workspace_validation.py" --input "$input" \
  --output "$output/likelihood_validation.json" --bestfit "$bestfit" --weighted-asimov \
  --range 0 2 --scan-points "$output/scan_points.json" "${covariance_options[@]}"
