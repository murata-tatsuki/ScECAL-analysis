#!/bin/sh
path="$PWD"
cd $path

# beam optimzation search simulation

particle=mu-
energy=100



condor_out="$path/out"
condor_err="$path/err"
condor_log="$path/log"
condor_cmd="$path/cmd"
condor_sh="$path/run"
condor_mac="$path/mac"
# condor_root="$path/../../Result_MC/generate/${particle}/${energy}GeV"
condor_root="$path/../../Result_MC/generate/${particle}/mu_track_slope_2_14/${energy}GeV"

# rm -rf $condor_out $condor_err $condor_log $condor_cmd $condor_sh $condor_mac $condor_root
mkdir -p $condor_out $condor_err $condor_log $condor_cmd $condor_sh $condor_mac $condor_root
mkdir -p ${condor_mac}/${particle}/${energy}GeV

# mkdir -p $condor_root

data_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result.bak/trackFit/beamsize_track/sps/mu-
# data_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result.bak/trackFit/beamsize_track/sp/mu-


# posX=$((-0.1846+50))
# posY=$((0.7962+50))
i=0
# posX=$((-50*i))
# posY=$((-50+j))

for txt in $( ls ${data_path} ); do
	# prefix=`basename $txt`
  prefix="${txt:0:-21}"
	# echo ${prefix}

#--------Create mac file------------------------------
	# echo "Create ${n} mac file"
	echo \
	"
	/control/verbose 0
	/run/verbose 0
	/event/verbose 0
	
	/mydet/generator fGParticleSource
	/gps/particle ${particle}
	/gps/energy ${energy} GeV
	/gps/direction 0 0 1
	/gps/direction 0.002 0.014 1

	/mydet/hitmapFile ${data_path}/${txt}
	
	/run/beamOn 100000

	" > $condor_mac/${particle}/${energy}GeV/${prefix}.mac

	sbatch -o out/test-%A.out --error="err/test-%A.err" execute.sl /home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/CEPCScECAL/build/CEPCScECAL ${condor_mac}/${particle}/${energy}GeV/${prefix}.mac $condor_root/${prefix}.root
	# sbatch -o out/test-%A_%a.out --error="err/test-%A_%a.err" execute.sl /home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/CEPCScECAL/build/CEPCScECAL ${condor_mac}/${particle}/${energy}GeV/${prefix}.mac $condor_root/${prefix}.root

	
#--------Create .sh file end------------------------------
	# echo \
	# "#!/bin/sh
	
	# source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc11-opt/setup.sh
	# /home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/CEPCScECAL/build/CEPCScECAL  $condor_mac/kek_${particle}_${energy}.mac $condor_root/kek_${particle}_${energy}.root
	
	# " > $condor_sh/${result}.sh
	# chmod +x $condor_sh/${result}.sh	
	# hep_sub -os CentOS7  $condor_sh/${result}.sh -g higgs -o $condor_out -e $condor_err
	# let "n+=1"

	let i++
done
# hadd merge_root electron_10000_40GeV.root  (from) generate_root 

#--------Creat cmd file------------------------------
#	if test -s $condor_cmd/${result}.cmd; 
#	then
#		echo "CMD File is not Empty"
#		rm  $condor_cmd/${result}.cmd
#	fi
#	cat >> $condor_cmd/${result}.cmd << EOF
#	  Universe             = vanilla
#	  Notification         = Never
#	  GetEnv               = True
#	  Executable           = $condor_sh/${result}.sh
#	  Output               = $condor_out/${result}.out
#	  Error                = $condor_err/${result}.err
#	  Log                  = $condor_log/${result}.log
#	  Queue								 
#EOF
#	
#	chmod +x $condor_sh/${result}.sh
#	chmod +x $condor_cmd/${result}.cmd 
#	
#	echo "condor_submit $condor_cmd/${result}.cmd"
#	condor_submit $condor_cmd/${result}.cmd 
#	
#	let "n+=1"
#done
