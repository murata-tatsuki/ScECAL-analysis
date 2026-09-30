#!/bin/bash
set -euo pipefail
HERE=$(dirname -- "$(realpath -e -- "${BASH_SOURCE[0]}")")
tasks=${1:?Missing task manifest};shift
index=${SLURM_ARRAY_TASK_ID:?Run this script through the plot job array}
[[ $index =~ ^[0-9]+$ ]]||exit 1
line=$(sed -n "$((10#$index+1))p" "$tasks")
IFS=$'\t' read -r energy layer <<< "$line"
[[ -n $energy && $layer =~ ^[0-9]+$ ]]||{ echo 'Invalid/missing array task' >&2;exit 1; }
printf 'Plot task %s: %s GeV, layer %s\n' "$index" "$energy" "$layer"
exec bash "$HERE/run.sh" --energy "$energy" "$@" --plots-only --layers "$layer"
