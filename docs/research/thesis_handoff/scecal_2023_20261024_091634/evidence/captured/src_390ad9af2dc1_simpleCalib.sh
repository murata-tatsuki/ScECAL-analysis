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

simpleCalib_o=${simpleCalib_path}

# MIP Fit
mip_i=${simpleCalib_o}

#------ Create Directores -----

#mkdir -p ${simpleCalib_o}

#------ simple calibration -----
#for i in $(ls ${decode_o})
#do
#        file_name=`basename $i`
#        prefix=${file_name%.root}
#        if [[ -e ${mip_i}/${prefix}.root ]]
#        then
#                echo "${mip_i}/${prefix} has been decoded, skip ..."
#                continue
#        fi
#
#        echo "simple calibration: ${file_name}"
#        echo " ... ... "
#        ${bin_dir}/simpleCalib ${decode_o}/${prefix}.root "${pedestal_o_2}/" "${simpleCalib_o}/" # > ${log_dir}/${prefix}.mip.log
#
#        echo "=============================="
#done


if [[ "${energy}" != "all" ]]
then
        mkdir -p ${simpleCalib_o}
        for i in $(ls ${decode_o})
        do
                file_name=`basename $i`
                prefix=${file_name%.root}

                if [[ -e ${simpleCalib_o}/${prefix}.root ]]
                  then
                          echo "${simpleCalib_o}/${prefix} has been decoded, skip ..."
                          continue
                  fi

                echo "simple Calibration ${prefix}"
                echo " ... ... "
                sbatch -o output/test-%A.out --error="output/test-%A.err" execute.sl ${bin_dir}/simpleCalib ${decode_o}/${prefix}.root "${pede_factor}" "${simpleCalib_o}/" # > ${log_dir}/${prefix}.mip.log
                # echo ${bin_dir}/simpleCalib ${decode_o}/${prefix}.root "${pedestal_o_2}/" "${simpleCalib_o}/" # > ${log_dir}/${prefix}.mip.log
                echo "    Output dir ===> ${simpleCalib_o}/${prefix}.root "

                echo "=============================="
        done
fi

ene=0
if [[ "${energy}" == "all" ]]
then
        decode_o=`dirname ${decode_o}`
        for i in $(ls ${decode_o}/*)
        do
                if [[ ${i} == *"GeV"* ]]
                then
                        ene=`basename $i`
                        continue
                fi

                file_name=`basename $i`
                prefix=${file_name%.root}
                simpleCalib_ou=`dirname ${simpleCalib_o}`/${ene/%?/}
                decode_ou=${decode_o}/${ene/%?/}
                pedestal_ou_2=`dirname ${pedestal_o_2}`/${ene/%?/}

                mkdir -p ${simpleCalib_ou}

                if [[ -e ${simpleCalib_ou}/${prefix}.root ]]
                  then
                          echo "${simpleCalib_ou}/${prefix} has been decoded, skip ..."
                          continue
                  fi

                echo "simple Calibration ${prefix}"
                echo " ... ... "
                ${bin_dir}/simpleCalib ${decode_ou}/${prefix}.root "${pedestal_ou_2}/" "${simpleCalib_ou}/" # > ${log_dir}/${prefix}.mip.log
                echo "    Output dir ===> ${simpleCalib_ou}/${prefix}.root "

                echo "=============================="
        done
fi
