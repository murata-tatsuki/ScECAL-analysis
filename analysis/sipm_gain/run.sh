#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0
(flock 9; make --no-print-directory -s -C "$HERE" all) 9>"$HERE/.build.lock"
exec "$HERE/Analyze" "$@"
