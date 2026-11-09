#! /bin/bash
slurm_log_dir=$(realpath -m -- "$(dirname -- "${BASH_SOURCE[0]}")/../log/calibration_parameters")
mkdir -p -- "$slurm_log_dir"

source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config


mkdir -p ../result/calibrations/figures/mip/raw
mkdir -p ../result/calibrations/figures/pedestal/raw
mkdir -p ../result/calibrations/figures/threshold/raw
mkdir -p ../result/calibrations/figures/gain/raw


sbatch -o "${slurm_log_dir}/test-%A.out" --error="${slurm_log_dir}/test-%A.err" execute.sl ./calibration_drawing ../result/calibrations/out.root
