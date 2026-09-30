#!/bin/bash
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --job-name=comparison

set -euo pipefail
ulimit -c 0
exec "$@"
