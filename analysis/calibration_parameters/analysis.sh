#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config


mkdir -p ../result/calibrations/figures/mip/raw
mkdir -p ../result/calibrations/figures/pedestal/raw
mkdir -p ../result/calibrations/figures/threshold/raw
mkdir -p ../result/calibrations/figures/gain/raw


sbatch -o jobs/test-%A.out --error="jobs/test-%A.err" execute.sl ./calibration_drawing ../result/calibrations/out.root
