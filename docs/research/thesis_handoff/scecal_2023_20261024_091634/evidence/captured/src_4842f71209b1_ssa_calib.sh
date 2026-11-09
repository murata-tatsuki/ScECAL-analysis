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

# Strip Split
stripsplit_o=${stripSprit_path}
stripsplit_o=${calibStripSprit_path}

#------ Create Directores -----
#mkdir -p ${stripsplit_o}

##------ Decode -----
#for i in $(ls ${calib_o})
#do
#	file_name=`basename $i`
#	prefix=${file_name%.root}
#
#	if [[ -e ${stripsplit_o}/${prefix}.root ]]
#          then
#                  echo "${stripsplit_o}/${prefix} has been decoded, skip ..."
#                  continue
#        fi
#
#	echo "Strip Splitting Algorithm ${prefix}"
#	echo " ... ... "
#	${bin_dir}/ssa ${file_name} "${calib_o}/" "${stripsplit_o}/"
#	echo "    Output dir ===> ${stripsplit_o}/${prefix}.root "
#
#	echo "=============================="
#done
#
simpleCalib_o=${simpleCalib_path}
#------ SSA -----
#for i in $(ls ${simpleCalib_o})
#do
#        file_name=`basename $i`
#        prefix=${file_name%.root}
#
#        if [[ -e ${stripsplit_o}/${prefix}.root ]]
#          then
#                  echo "${stripsplit_o}/${prefix} has been decoded, skip ..."
#                  continue
#        fi
#
#        echo "Strip Splitting Algorithm ${prefix}"
#        echo " ... ... "
#        ${bin_dir}/ssa ${file_name} "${simpleCalib_o}/" "${stripsplit_o}/"
#        echo "    Output dir ===> ${stripsplit_o}/${prefix}.root "
#
#        echo "=============================="
#done



# if [[ -e ${decode_o}/${date_tag}.root ]]
# then
# 	rm ${decode_o}/${date_tag}.root
# 	hadd ${decode_o}/${date_tag}.root ${decode_o}/*
# 	rm -i ${decode_o}/*_*.root
# else
# 	hadd ${decode_o}/${date_tag}.root ${decode_o}/*
# 	rm -i ${decode_o}/*_*.root
# fi








if [[ "${energy}" != "all" ]]
then
        mkdir -p ${stripsplit_o}
        #for i in $(ls ${simpleCalib_o})
        for i in $(ls ${calib_o})
        do
                file_name=`basename $i`
                prefix=${file_name%.root}

                if [[ -e ${stripsplit_o}/${prefix}.root ]]
                  then
                          echo "${stripsplit_o}/${prefix} has been decoded, skip ..."
                          continue
                  fi

                echo "strip splitting ${prefix}"
                echo " ... ... "
                sbatch -o output/test-%A.out --error="output/test-%A.err" execute.sl ${bin_dir}/ssa ${file_name} "${calib_o}/" "${stripsplit_o}/"
                echo "    Output dir ===> ${stripsplit_o}/${prefix}.root "

                echo "=============================="
        done
fi

echo $stripsplit_o
echo $calib_o

ene=0
if [[ "${energy}" == "all" ]]
then
        calib_o=`dirname ${calib_o}`
        #simpleCalib_o=`dirname ${simpleCalib_o}`
        #for i in $(ls ${simpleCalib_o}/*)
        for i in $(ls ${calib_o}/*)
        do
                if [[ ${i} == *"GeV:"* ]]
                then
                        ene=`basename $i`
                        continue
                fi

                file_name=`basename $i`
                prefix=${file_name%.root}
                stripsplit_ou=`dirname ${stripsplit_o}`/${ene/%?/}
                #simpleCalib_ou=${simpleCalib_o}/${ene/%?/}
                calib_ou=${calib_o}/${ene/%?/}
                # pedestal_ou_2=`dirname ${pedestal_o_2}`/${ene/%?/}

                mkdir -p ${stripsplit_ou}

                if [[ -e ${stripsplit_ou}/${prefix}.root ]]
                  then
                          echo "${stripsplit_ou}/${prefix} has been decoded, skip ..."
                          continue
                  fi

                echo "strip splitting ${prefix}"
                echo " ... ... "
                sbatch -o output/test-%A.out --error="output/test-%A.err" execute.sl ${bin_dir}/ssa ${file_name} "${calib_ou}/" "${stripsplit_ou}/"
                echo "    Output dir ===> ${stripsplit_ou}/${prefix}.root "

                echo "=============================="
        done
fi

