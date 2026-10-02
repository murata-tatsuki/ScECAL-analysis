#!/bin/bash
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --job-name=cog_channel_scan
#SBATCH --output=/home/murata_t/ScECAL_BeamTest/analysis/log/cog_channel_scan/slurm-%j.out
#SBATCH --error=/home/murata_t/ScECAL_BeamTest/analysis/log/cog_channel_scan/slurm-%j.err

set -euo pipefail
ulimit -c 0
exec "$@"
