#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
logs="$HERE/../../log/channel_response/mip_closure"
mkdir -p -- "$logs"
exec sbatch --job-name=mip_closure --output="$logs/%j.out" --error="$logs/%j.err" \
    "$HERE/../execute.sl" bash "$HERE/run.sh" "$@"
