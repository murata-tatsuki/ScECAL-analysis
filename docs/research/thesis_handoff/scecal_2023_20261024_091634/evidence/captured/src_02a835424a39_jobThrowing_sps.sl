#!/bin/bash
###############################################################################
# A sample SLURM script to run analyzer on 1 node
#
#  usage:
#        (0) Edit this script
#        (1) Submit the job.
#            % cd $MEG2SYS/analyzer
#            % sbatch -a 0-20%10 example/analyzer_array.sl
#             sbatch -o out/test-%A_%a.out --error="err/test-%A_%a.err" jobThrowing_sps.sl 
#
###############################################################################

#SBATCH --partition=all
#SBATCH --ntasks=1       
#SBATCH --job-name=analyzer
#SBATCH --array=0
#SBATCH --output="analyzer-%A_%a.out"
#SBATCH --error="analyzer-%A_%a.err"

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

particle=e-
energy=$1
ijob=$2
angle=$3

path="$PWD"
condor_mac="$path/mac"
condor_root="$path/../../Result_MC/generate/${particle}/mu_track/${energy}GeV"
mkdir -p $condor_root

data_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/trackFit/beamsize_track/ps/mu-


# posX=$((-0.1846+50))
# posY=$((0.7962+50))
i=0
# posX=$((-50*i))
# posY=$((-50+j))

for txt in $( ls ${data_path} ); do
	# prefix=`basename $txt`
  prefix="${txt:0:-21}"
	echo ${prefix}

  /home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/CEPCScECAL/build/CEPCScECAL ${condor_mac}/${particle}/${energy}GeV/${prefix}.mac $condor_root/${prefix}.root

done

# 		let j++
# 	done
# 	let i++
# done

RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE ($SECONDS sec)"
