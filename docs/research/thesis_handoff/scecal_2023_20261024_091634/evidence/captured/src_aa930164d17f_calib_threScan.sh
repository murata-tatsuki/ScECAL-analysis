#! /bin/bash
source ./global_config

#------ I/O directories for each processor -------
# Decode: Raw to ROOT
decode_i=${raw_path}
decode_o=${decode_path}

# Pedestal
pedestal_i=${decode_path}
pedestal_o_1=${ped_hist}
pedestal_o_2=${ped_fit}

# MIP Fit
mip_i=${decode_path}
mip_o_1=${mip_hist}
mip_o_2=${mip_fit}

# Calibration
calib_o=${calib_path}

#------ Create Directores -----
#mkdir -p ${calib_o}

#------ Decode -----
#for i in $(ls ${decode_o})
#do
#  file_name=`basename $i`
#  prefix=${file_name%.root}
#  echo "Calibration ${prefix}"
#  echo " ... ... "
#  ${bin_dir}/calib ${decode_o}/${file_name} ${pede_factor} ${inter_factor} ${mip_factor} ${calib_o}/${prefix}.root #> ${log_dir}/${prefix}.log
#  echo "    Output dir ===> ${calib_o}/${prefix}.root "
#
#  echo "=============================="
#done


if [[ "${energy}" != "all" ]]
then
  mkdir -p ${calib_o}
  for i in $(ls ${decode_o})
  do
    file_name=`basename $i`
    prefix=${file_name%.root}

    if [[ -e ${calib_o}/${prefix}.root ]]
      then
        echo "${calib_o}/${prefix} has been decoded, skip ..."
        continue
      fi
    
    sbatch -o output/test-%A_%a.out --error="output/test-%A_%a.err" execute_analyzer_calib.sl ${decode_o}/${file_name} ${pede_factor} ${inter_factor} ${mip_factor} ${calib_o}/${prefix}.root

    # echo "Calibration ${prefix}"
    # echo " ... ... "
    # #${bin_dir}/calib ${decode_o}/${file_name} ${pedestal_o_2}/${file_name} ${inter_factor} ${mip_factor} ${calib_o}/${prefix}.root #> ${log_dir}/${prefix}.log
    # ${bin_dir}/calib ${decode_o}/${file_name} ${pede_factor} ${inter_factor} ${mip_factor} ${calib_o}/${prefix}.root #> ${log_dir}/${prefix}.log
    # echo "    Output dir ===> ${calib_o}/${prefix}.root "

    # echo "=============================="
  done
fi

ene=0
if [[ "${energy}" == "all" ]]
then
  decode_o=`dirname ${decode_o}`
  for i in $(ls ${decode_o}/*)
  do
    if [[ ${i} == *"GeV:" ]]
    then
      ene=`basename $i`
      continue
    fi

    # file_name=`basename $i`
    # prefix=${file_name%.root}
#     calib_ou=`dirname ${calib_o}`/${ene/%?/}/threshold
    # decode_ou=${decode_o}/${ene/%?/}
    # pedestal_ou_2=`dirname ${pedestal_o_2}`/${ene/%?/}
    # echo $pedestal_ou_2

    threName=`basename $i`
    echo ${threName}

    calib_ou=`dirname ${calib_o}`/${threName}/${ene/%?/}
    decode_ou=${decode_o}/${ene/%?/}/${threName}

    mkdir -p ${calib_ou}

    jobOut=output/scan/${date_tag}/${ene/%?/}/out
    jobErr=output/scan/${date_tag}/${ene/%?/}/err
    # rm -r ${jobOut}
    # rm -r ${jobErr}
    mkdir -p ${jobOut}
    mkdir -p ${jobErr}

    for rootfile in $(ls ${decode_ou}/*)
    do
      file_name=`basename $rootfile`
      prefix=${file_name%.root}

      if [[ -e ${calib_ou}/${prefix}.root ]]
      then
        echo "${calib_ou}/${prefix} has been decoded, skip ..."
        continue
      fi

      sbatch -o ${jobOut}/${prefix}-%A_%a.out --error="${jobErr}/${prefix}-%A_%a.err" execute_analyzer_calib.sl ${decode_ou}/${file_name} ${pede_factor} ${inter_factor} ${mip_factor} ${calib_ou}/${prefix}.root

    done
  #   echo "Calibration ${prefix}"
  #   echo " ... ... "
  #   ${bin_dir}/calib ${decode_ou}/${file_name} ${pede_factor} ${inter_factor} ${mip_factor} ${calib_ou}/${prefix}.root #> ${log_dir}/${prefix}.log
  #   ## ${bin_dir}/calib ${decode_ou}/${file_name} ${pedestal_ou_2}/${file_name} ${inter_factor} ${mip_factor} ${calib_ou}/${prefix}.root #> ${log_dir}/${prefix}.lo
  # #       ${bin_dir}/calib ${decode_o}/${file_name} ${pede_factor} ${inter_factor} ${mip_factor} ${calib_o}/${prefix}.root #> ${log_dir}/${prefix}.log
  #   echo "    Output dir ===> ${calib_ou}/${prefix}.root "

  #   echo "=============================="
  done
fi





# if [[ -e ${decode_o}/${date_tag}.root ]]
# then
#   rm ${decode_o}/${date_tag}.root
#   hadd ${decode_o}/${date_tag}.root ${decode_o}/*
#   rm -i ${decode_o}/*_*.root
# else
#   hadd ${decode_o}/${date_tag}.root ${decode_o}/*
#   rm -i ${decode_o}/*_*.root
# fi
