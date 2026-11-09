#! /bin/bash
# source ./global_config
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/seino_config

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

#------ Create Directores -----


if [ "${date_tag}" == "pedestal" ] || [ "${date_tag}" == "Pedestal" ]
then
	# pedestal_o_1=`dirname ${pedestal_o_1}`
        # pedestal_o_2=`dirname ${pedestal_o_2}`
	# pedestal_i=`dirname ${pedestal_i}`
	mkdir -p ${pedestal_o_1}
        mkdir -p ${pedestal_o_2}
	for i in $(ls ${pedestal_i})
	do
		file_name=`basename $i`
		prefix=${file_name%.root}

		if [[ -e ${pedestal_o_1}/${prefix}.root && -e ${pedestal_o_2}/${prefix}.root ]]
                then
                        echo "\[${prefix}\] has been decoded, skip ..."
                        continue
                fi

                if [[ ! -e ${pedestal_i}/${prefix}.root ]]
                then
                        echo "${pedestal_i}/${prefix}.root does NOT exist, skip ..."
                        continue
                fi

                #------ Diagnose -----
                echo "Fitting Pedestal: ${file_name}"
                echo " ... ... "
                ${bin_dir}/pedestal_pedestal ${pedestal_i}/${prefix}.root "${pedestal_o_1}/" "${pedestal_o_2}/" #> ${log_dir}/${prefix}.ped.log
                echo "    Output histogram ===> ${pedestal_o_1}/${prefix}.root"
                echo "    Output fit result ===> ${pedestal_o_2}/${prefix}.root"

                echo "=============================="
        done
fi

if [ "${date_tag}" != "pedestal" ] && [ "${date_tag}" != "Pedestal" ]
then
	echo "please change date_tag to pedestal or Pedestal"
fi


# if [[ -e ${decode_o}/${date_tag}.root ]]
# then
# 	rm ${decode_o}/${date_tag}.root
# 	hadd ${decode_o}/${date_tag}.root ${decode_o}/*
# 	rm -i ${decode_o}/*_*.root
# else
# 	hadd ${decode_o}/${date_tag}.root ${decode_o}/*
# 	rm -i ${decode_o}/*_*.root
# fi
#！/bin/sh

