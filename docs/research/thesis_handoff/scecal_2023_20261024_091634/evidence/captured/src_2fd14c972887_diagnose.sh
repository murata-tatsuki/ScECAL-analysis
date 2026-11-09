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


#------ Create Directores -----
mkdir -p ${decode_o}

mkdir -p ${pedestal_o_1}
mkdir -p ${pedestal_o_2}

mkdir -p ${mip_o_1}
mkdir -p ${mip_o_2}

mkdir -p ${log_dir}

#------ Decode -----
for i in $(ls ${raw_path})
do
	file_name=`basename $i`
	prefix=${file_name%.dat}
	if [[ -e ${decode_o}/${prefix}.root ]]
	then
		echo "\[${prefix}\] has been decoded, skip ..."
		continue
	fi
	echo "Decoding ${prefix}"
	echo " ... ... "
	${bin_dir}/decode ${decode_i}/${file_name} ${decode_o}/${prefix}.root > ${log_dir}/${prefix}.log
	echo "    Output dir ===> ${decode_o}/${prefix}.root "

	#------ Diagnose -----
	echo "Fitting Pedestal: ${file_name}"
	echo " ... ... "
	${bin_dir}/pedestal ${pedestal_i}/${prefix}.root "${pedestal_o_1}/" "${pedestal_o_2}/" > ${log_dir}/${prefix}.ped.log
	echo "    Output histogram ===> ${pedestal_o_1}/${prefix}.root"
	echo "    Output fit result ===> ${pedestal_o_2}/${prefix}.root"

	echo "Fitting MIP: ${file_name}"
	echo " ... ... "
	${bin_dir}/mip ${mip_i}/${prefix}.root "${mip_o_1}/" "${mip_o_2}/" > ${log_dir}/${prefix}.mip.log
	echo "    Output histogram ===> ${mip_o_1}/${prefix}.root"
	echo "    Output fit result ===> ${mip_o_2}/${prefix}.root"

	echo "=============================="
done

# if [[ -e ${decode_o}/${date_tag}.root ]]
# then
# 	rm ${decode_o}/${date_tag}.root
# 	hadd ${decode_o}/${date_tag}.root ${decode_o}/*
# 	rm -i ${decode_o}/*_*.root
# else
# 	hadd ${decode_o}/${date_tag}.root ${decode_o}/*
# 	rm -i ${decode_o}/*_*.root
# fi
