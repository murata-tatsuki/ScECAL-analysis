#/bin/bash

workDir=${PWD}
echo "========== Work Dir : ${workDir}"
particleTagl=e-/e_ssa
# particleTag=mu-/mu_track
# energyTag=GeV
particleTag=e-/beam_scan/point
particleTag=e-/beam_scan/point_noTriggerSci

#resultDir="/home/wangjx/HEP/CEPC2023/Result_MC"
resultDir="/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC"
# truthDir="${resultDir}/truth/${particleTag}/${energyTag}"
# # digiDir="${resultDir}/decode/${particleTag}_new_edit/${energyTag}"
# digiDir="${resultDir}/decode/${particleTag}/threScan/${energyTag}"

shareDir="${workDir}/share"
pedFile="${shareDir}/pedestal2023_SPS.root"
HLRatioFile="${shareDir}/all_hl_electron2023_hlratio_v3.root"
MIPFile="${shareDir}/all_auto_muon_v4_trackfit.root"
#MIPFile="${shareDir}/all_auto_muon_v4.root"
SiPMGainFile="${shareDir}/20230504_1425_Calib_LED_all.root"
# ThresholdFile="${shareDir}/threshold_rev.root"
# ThresholdFile="notreverse"
ThresholdFile="/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/threshold/threshold.root"

## cosmic ray parameters
# shareDir="${workDir}/share/cosmic_parameters"
# pedFile="${shareDir}/2020ped.root"
# HLRatioFile="${shareDir}/20201216_20210308_InterCalib.root"
# MIPFile="${shareDir}/20201216_20210308.root"
# SiPMGainFile="${shareDir}/20210702_160312_LED_Calib.root"

echo "... ..."

# ${workDir}/build/MCAnalysis ${mergeDir}/${fileName} ${truthDir}/${fileName} ${energyAnaTruthDir}/${fileName} ${energyAnaTruthHistDir}/${fileName} ${digiDir}/${fileName} ${energyAnaDigiDir}/${fileName} ${energyAnaDigiHistDir}/${fileName} ${pedFile} ${HLRatioFile} ${MIPFile} ${SiPMGainFile}
# ${workDir}/build/MCAnalysis  ${truthDir}/${fileName} ${pedFile} ${HLRatioFile} ${MIPFile} ${SiPMGainFile} ${digiDir}/${fileName}

# for i in `seq 0 10`; do
for energy in $( ls ${resultDir}/generate/${particleTag} ); do   
  echo ${energy}
  energyTag=${energy}


  truthDir="${resultDir}/truth/${particleTag}/${energyTag}"
  digiDir="${resultDir}/decode/${particleTag}/${energyTag}"
	mkdir -p "$truthDir"
	mkdir -p "$digiDir"

  jobOut=job/scan/${particleTag}/${energyTag}/out
  jobErr=job/scan/${particleTag}/${energyTag}/err
  rm -r ${jobOut}
  rm -r ${jobErr}
  mkdir -p ${jobOut}
  mkdir -p ${jobErr}

  for fileName in `ls ${resultDir}/generate/${particleTag}/${energyTag}`; do

    echo ${fileName}
    # cp ${resultDir}/generate/${particleTag}/${energyTag}/${fileName} ${truthDir}

    # for i in `seq 0 10`; do
      # a=$((i*10))
      mkdir -p job/sl/scan/${particleTag}/${energyTag}
      # mkdir -p ${digiDir}/threshold
      # outsl=job/sl/scan/${particleTag}/sps/${energyTag}/${fileName}_thre.sl
      # cp job/sl/example.sl ${outsl}

    	outfileName=${fileName}
    	# outfileName=${particleTag}_${energyTag}_100000_thre_${a}.root        # pedestal + ${a}

      Arg7=${digiDir}/${fileName}

      if [[ -e ${Arg7} ]]
      then
        echo "${Arg7} has been decoded, skip ..."
        # continue
      fi
      
      sbatch -o ${jobOut}/${fileName}-%A.out --error="${jobErr}/${fileName}-%A.err" execute.sl ${workDir}/build/MCAnalysis ${resultDir}/generate/${particleTag}/${energyTag}/${fileName} ${pedFile} ${HLRatioFile} ${MIPFile} ${SiPMGainFile} ${ThresholdFile} ${Arg7} 0
      # sbatch -o ${jobOut}/${fileName}-%A.out --error="${jobErr}/${fileName}-%A.err" execute.sl ${workDir}/build/MCAnalysis ${truthDir}/${fileName} ${pedFile} ${HLRatioFile} ${MIPFile} ${SiPMGainFile} ${ThresholdFile} ${Arg7} 0
      
    # done
  done
  sleep 3
done
