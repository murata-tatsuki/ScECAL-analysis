#!/bin/bash
set -euo pipefail
path="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$path"

# beam optimzation search simulation

particle=e-
energy=100
tagName=point_noTriggerSci  # point point_noTriggerSci
condor_out="$path/out"
condor_err="$path/err"
condor_log="$path/log"
condor_cmd="$path/cmd"
condor_sh="$path/run"
condor_mac="$path/mac"
# condor_root="$path/../../Result_MC/generate/${particle}/${energy}GeV"

# # rm -rf $condor_out $condor_err $condor_log $condor_cmd $condor_sh $condor_mac $condor_root
# mkdir -p $condor_out $condor_err $condor_log $condor_cmd $condor_sh $condor_mac $condor_root
# mkdir -p ${condor_mac}/${particle}/${energy}GeV

# condor_root="$path/../../Result_MC/generate/${particle}/mu_track/${energy}GeV"
# mkdir -p $condor_root

# data_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/trackFit/beamsize_ssa/sps/e-
# data_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/trackFit/beamsize_ssa/ps/e-


# posX=−6.7
# posY=4.9


# for txt in $( ls ${data_path}/*GeV* ); do
	# echo ${txt}
	# prefix=`basename $txt`
  # prefix="${prefix:0:-5}"
	# # echo ${prefix}

	# energy="${prefix:0:-3}"
	# echo $energy

	# if [ "$energy" -ne 150 ] && [ "$energy" -ne 200 ] && [ "$energy" -ne 250 ]; then
  #   continue
	# fi

	echo $energy
	prefix="${energy}GeV"

	condor_root="$path/../../Result_MC/generate/${particle}/beam_scan/$tagName/${energy}GeV"
	# rm -rf $condor_out $condor_err $condor_log $condor_cmd $condor_sh $condor_mac $condor_root
	mkdir -p $condor_out $condor_err $condor_log $condor_cmd $condor_sh $condor_mac $condor_root
	mkdir -p ${condor_mac}/${particle}/${energy}GeV
	condor_root="$path/../../Result_MC/generate/${particle}/beam_scan/$tagName/${energy}GeV"
	mkdir -p $condor_root


	nfiles=100
	events_per_file=1000
		for ((i=0; i<nfiles; i++)); do
  	  outfile="$condor_root/${prefix}_${i}.root"
  	  if [[ -e "$outfile" ]]; then
  	      echo "既存ファイルがあります: $outfile" >&2
  	      exit 1
  	  fi
		done


	i=0
	while [ "$i" -lt "$nfiles" ]
  do

	seed1=$((310000 + 2 * i))
	seed2=$((310001 + 2 * i))

#--------Create mac file------------------------------
	# echo "Create ${n} mac file"
	echo \
	"
	/control/verbose 0
	/run/verbose 0
	/event/verbose 0

	/random/setSeeds ${seed1} ${seed2}

	/mydet/generator fGParticleSource
	/gps/particle ${particle}
	/gps/energy ${energy} GeV
	/gps/direction 0 0 1
	/gps/pos/type Point

	/gps/pos/centre -6.7 4.9 -5 mm
	
	/run/beamOn ${events_per_file}

	" > $condor_mac/${particle}/${energy}GeV/${prefix}_${i}.mac

	sbatch -o out/test-%A.out --error="err/test-%A.err" execute.sl /home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/CEPCScECAL/build/CEPCScECAL ${condor_mac}/${particle}/${energy}GeV/${prefix}_${i}.mac $condor_root/${prefix}_${i}.root
	# sbatch -o out/test-%A_%a.out --error="err/test-%A_%a.err" execute.sl /home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/CEPCScECAL/build/CEPCScECAL ${condor_mac}/${particle}/${energy}GeV/${prefix}.mac $condor_root/${prefix}.root

		i=$((i + 1))
	done
# done
