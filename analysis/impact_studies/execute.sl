#!/bin/bash
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --job-name=impact_studies
#SBATCH --output=/home/murata_t/ScECAL_BeamTest/analysis/log/impact_studies/slurm-%j.out
#SBATCH --error=/home/murata_t/ScECAL_BeamTest/analysis/log/impact_studies/slurm-%j.err
set -euo pipefail
ulimit -c 0
exec "$@"
