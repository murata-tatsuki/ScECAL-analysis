#！/bin/sh

source ./global_config

layer=8

input_file=${simpleCalib_path}
trackfile=${data_dir}/${year}/${beam}/trackFit_layer16/${date_tag}/${energy}


#echo $datafile
#./main ${DIR_IN} ${DIR_OUT1} ${layer} 
mkdir -p ${trackfile}

 ##### for file loop
 num=0
for dat_name in $(ls ${input_file})
  do


  prefix=`basename $dat_name`
  if [[ -e ${trackfile}/${prefix} ]]
  then
                echo "\[${prefix}\] has been decoded, skip ..."
                continue
  fi

	echo ${bin_dir}/track_fit_ ${input_file}/${dat_name} ${trackfile}/ ${layer}
  done

