#!/bin/bash
# Shared Slurm launcher; submit.sh supplies the analysis-specific job name/logs.
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=channel_response
#SBATCH --output=/home/murata_t/ScECAL_BeamTest/analysis/log/channel_response/slurm-%j.out
#SBATCH --error=/home/murata_t/ScECAL_BeamTest/analysis/log/channel_response/slurm-%j.err
set -euo pipefail
ulimit -c 0
exec "$@"
