#!/bin/bash
###############################################################################
# A sample SLURM script to run analyzer on 1 node
#
#  usage:
#        (0) Edit this script
#        (1) Submit the job.
#            % cd $MEG2SYS/analyzer
#            % sbatch -o test-%A_%a.out --error="test-%A_%a.err" execute_analyzer_ssa.sl
#
###############################################################################

#SBATCH --partition=all
#SBATCH --ntasks=1       
#SBATCH --job-name=analyzer   
#SBATCH --output="analyzer%A.out"
#SBATCH --error="analyzer%A.err"

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
./ssa.sh

RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE"
