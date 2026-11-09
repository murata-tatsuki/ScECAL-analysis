#!/bin/bash
###############################################################################
# A sample SLURM script to run analyzer on 1 node
#
#  usage:
#        (0) Edit this script
#        (1) Submit the job.
#            % cd $MEG2SYS/analyzer
#            % sbatch -a 0-20%10 example/analyzer_array.sl
#
###############################################################################

#SBATCH --partition=all
#SBATCH --ntasks=1       
#SBATCH --job-name=test
#SBATCH --array=0
#SBATCH --output="test-%A_%a.out"
#SBATCH --error="test-%A_%a.err"

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
# Replace the start run number and configuration file by yours
#
START_RUNNUM=1
RUNNUM=$(( START_RUNNUM + SLURM_ARRAY_TASK_ID ))
#./meganalyzer -i example/bartenderMode.xml -r $RUNNUM -b
echo "Process ID: $SLURM_ARRAY_TASK_ID"
sleep 90
echo "Process ID end: $SLURM_ARRAY_TASK_ID"

RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE ($SECONDS sec)"
