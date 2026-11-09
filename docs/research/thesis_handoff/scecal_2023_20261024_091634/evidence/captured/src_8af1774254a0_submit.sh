#!/bin/bash
path="$PWD"
cd $path

particle=e-
condor_out="$path/out"
condor_err="$path/err"
condor_log="$path/log"
condor_cmd="$path/cmd"
condor_sh="$path/run"
condor_mac="$path/mac"

# data_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/trackFit/beamsize_ssa/sps/e-
data_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/trackFit/beamsize_ssa/ps/e-

for txt in $( ls ${data_path}/*GeV* ); do
    prefix=`basename $txt`
    prefix="${prefix:0:-5}"
    energy="${prefix:0:-3}"

		# if [ "$energy" -ne 150 ] && [ "$energy" -ne 200 ] && [ "$energy" -ne 250 ]; then
	  #   continue
		# fi

		echo $energy

    condor_root="$path/../../Result_MC/generate/${particle}/e_ssa/${energy}GeV"
    mkdir -p $condor_out $condor_err $condor_log $condor_cmd $condor_sh $condor_mac $condor_root
    mkdir -p ${condor_mac}/${particle}/${energy}GeV

    # 1つのtxtファイルにつき100個のアレイジョブ(0〜99)を一括投入
    sbatch --array=0-99 \
           -o "out/e/test-%A_%a.out" \
           --error="err/e/test-%A_%a.err" \
           execute_e.sl "$txt" "$particle" "$energy" "$prefix" "$condor_root" "${condor_mac}/${particle}/${energy}GeV"
done