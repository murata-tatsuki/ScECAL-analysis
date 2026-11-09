#! /bin/bash
source ./global_config

particle=e+		######<<<<<CHANGE HERE ######
dat=ECAL_Run195_20221027_025512.dat		######<<<<<CHANGE HERE ######

#ln -s /mnt2/BeamData/ECAL/${particle}/${dat} ${raw_path}/

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
#mkdir -p ${decode_o}

#------ Decode -----
if [[ "${energy}" != "all" && "${date_tag}" != "pedestal" ]]
then
	mkdir -p ${decode_o}
	for i in $(ls ${decode_i})
	do
		file_name=`basename $i`
		prefix=${file_name%.dat}
	
		grep -q "${prefix}.root" auto_gain_files.txt
		if [[ $? -eq 0 ]]
		then
			echo "\[${prefix}\] is auto gain data, skip ..."
			continue
		fi
		
		if [[ -e ${decode_o}/${prefix}.root ]]
		then
			echo "${decode_o}/${prefix} has been decoded, skip ..."
			continue
		fi

		echo "Decoding ${prefix}"
		echo " ... ... "
		${bin_dir}/decode ${decode_i}/${file_name} ${decode_o}/${prefix}.root #> ${log_dir}/${prefix}.log
		echo "    Output dir ===> ${decode_o}/${prefix}.root "

		echo "=============================="
	done
fi

if [[ "${date_tag}" == "Pedestal" ]]
then
	decode_o=`dirname ${decode_o}`
	decode_i=`dirname ${decode_i}`
        mkdir -p ${decode_o}
        for i in $(ls ${decode_i})
        do
                file_name=`basename $i`
                prefix=${file_name%.dat}
        
                grep -q "${prefix}.root" auto_gain_files.txt
                if [[ $? -eq 0 ]]
                then
                        echo "\[${prefix}\] is auto gain data, skip ..."
                        continue
                fi
                
                if [[ -e ${decode_o}/${prefix}.root ]]
                then
                        echo "${decode_o}/${prefix} has been decoded, skip ..."
                        continue
                fi

                echo "Decoding ${prefix}"
                echo " ... ... "
                ${bin_dir}/decode ${decode_i}/${file_name} ${decode_o}/${prefix}.root #> ${log_dir}/${prefix}.log
                echo "    Output dir ===> ${decode_o}/${prefix}.root "

                echo "=============================="
        done
fi


ene=0

if [[ "${energy}" == "all" && "${date_tag}" != "pedestal" ]]
then
	decode_i=`dirname ${decode_i}`
	for i in $(ls ${decode_i}/*)
	do
		if [[ ${i} != *".dat"* ]] 
		then
			ene=`basename $i`
			continue
		fi	
		file_name=`basename $i`
		prefix=${file_name%.dat}
		decode_ou=`dirname ${decode_o}`/${ene/%?/}
		decode_in=${decode_i}/${ene/%?/}

		mkdir -p ${decode_ou}

		grep -q "${prefix}.root" auto_gain_files.txt
		if [[ $? -eq 0 ]]
		then
			echo "\[${prefix}\] is auto gain data, skip ..."
			continue
		fi

		if [[ -e ${decode_ou}/${prefix}.root ]]
		then
			echo "${decode_ou}/${prefix} has been decoded, skip ..."
			continue
		fi

		echo "Decoding ${prefix}"
		echo " ... ... "
		${bin_dir}/decode ${decode_in}/${file_name} ${decode_ou}/${prefix}.root #> ${log_dir}/${prefix}.log
		echo "    Output dir ===> ${decode_ou}/${prefix}.root "
		echo "=============================="
	done
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
