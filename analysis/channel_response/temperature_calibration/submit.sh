#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
for arg in "$@"; do
    case "$arg" in -h|--help|--dry-run) exec bash "$HERE/run.sh" "$@";; esac
done
bash "$HERE/run.sh" "$@" --dry-run
logs="$HERE/../../log/channel_response/temperature_calibration"
mkdir -p -- "$logs"
exec sbatch --job-name=ps_temperature --chdir="$(pwd -P)" \
    --output="$logs/%j.out" --error="$logs/%j.err" \
    "$HERE/../execute.sl" bash "$HERE/run.sh" "$@"
