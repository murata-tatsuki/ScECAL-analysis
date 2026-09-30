#!/bin/bash
# Shared Slurm launcher; submit.sh supplies the analysis-specific job name/logs.
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=channel_response
set -euo pipefail
ulimit -c 0
exec "$@"
