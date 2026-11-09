#!/bin/bash
###############################################################################
# A sample SLURM script to run analyzer on 1 node
#
#  usage:
#        (0) Edit this script
#        (1) Submit the job.
#            % cd $MEG2SYS/analyzer
#            % sbatch example/analyzer.sl
#            % sbatch -o test-%A_%a.out --error="test-%A_%a.err" jobThrowing.sl
#
###############################################################################

#SBATCH --partition=all
#SBATCH --ntasks=1       
#SBATCH --job-name=analyzer   
#SBATCH --output="output/analyzer%A.out"
#SBATCH --error="ioutput/analyzer%A.err"

MY_HOST=`hostname`
MY_DATE=`date`
#
# Disable core dump. (Remove this line if you want to create core files)
#
ulimit -c 0

echo "Running on $MY_HOST at $MY_DATE"
echo "================================================================"

#
# Print envirionment
#
#echo "Running environment":
#env | sort
#echo "================================================================"

#
# Replace configuration file by yours
#
#./meganalyzer -i example/bartenderMode.xml -r 1 -b
# bash track_fit_.sh
# bash track_fit.sh

"$@"

# /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/bin/track_fit /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/simpleCalib/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run171_20250605_162329.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/trackFit/mipScan/layer_except_10_11_14_15_20_21_26_27/ 8



RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE"
