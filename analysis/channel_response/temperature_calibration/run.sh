#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec "${PYTHON:-python3.12}" "$HERE/run.py" "$@"
