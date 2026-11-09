#！/bin/sh

source ./global_config

layer=8

input_file=${simpleCalib_path}
# input_file=${calib_path}
trackfile=${trackFit_path}
# trackfile=${data_dir}/${year}/${beam}/trackFit/${date_tag}/${energy}
# trackfile=${data_dir}/${year}/${beam}/trackFit_layer16/${date_tag}/${energy}


#echo $datafile
#./main ${DIR_IN} ${DIR_OUT1} ${layer} 
mkdir -p ${trackfile}

 ##### for file loop
 num=0
for dat_name in $(ls ${input_file})
  do


  prefix=`basename $dat_name`
  # if [[ -e ${trackfile}/${prefix} ]]
  # then
  #   echo "\[${prefix}\] has been decoded, skip ..."
  #   continue
  # fi

  sbatch -o output/test-%A.out --error="output/test-%A.err" execute_analyzer_trackFit.sl ${bin_dir}/track_fit ${input_file}/${dat_name} ${trackfile}/ ${layer}

done

