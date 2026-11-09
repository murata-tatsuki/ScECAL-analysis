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
# bash calib.sh

source ./global_config

file_name=$1
pede_factor=$2
inter_factor=$3
mip_factor=$4
calib_o=$5

echo "Calibration ${5}"
echo " ... ... "
#${bin_dir}/calib ${decode_o}/${file_name} ${pedestal_o_2}/${file_name} ${inter_factor} ${mip_factor} ${calib_o}/${prefix}.root #> ${log_dir}/${prefix}.log
${bin_dir}/calib ${1} ${2} ${3} ${4} ${5} #> ${log_dir}/${prefix}.log
echo "    Output dir ===> ${5} "

echo "=============================="


RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE"
