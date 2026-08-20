#!/bin/bash
###############################################################################
# A sample SLURM script to run analyzer on 1 node
#
#  usage:
#        (0) Edit this script
#        (1) Submit the job.
#            % cd $MEG2SYS/analyzer
#            % sbatch example/analyzer.sl
#            % sbatch -o test-%A.out --error="test-%A.err" execute.sl
#            % sbatch -o jobs/test-%A.out --error="jobs/test-%A.err" execute.sl
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


track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/trackFit/mu-/100GeV
for dat_name in $(ls ${track_path})
do
  prefix=`basename $dat_name`
  mkdir -p ../../result/trackFit/efficiency/sps/${prefix}
done

track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/trackFit/mu-/10GeV/
for dat_name in $(ls ${track_path})
do
  prefix=`basename $dat_name`
  mkdir -p ../../result/trackFit/efficiency/ps/${prefix}
done

mkdir -p ../../result/trackFit/efficiency/sps/all
mkdir -p ../../result/trackFit/efficiency/ps/all


./channel_efficiecy


RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE"
