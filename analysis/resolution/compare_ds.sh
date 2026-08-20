#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

paralell_jobs () {
  Calib_path=$1
  cog_rage=$2

  for txt in $( ls ${Calib_path}/*_${cog_rage}mm.root ); do
    echo -n ${txt} >> ../tmp/compare_resolution.txt
    echo -n " " >> ../tmp/compare_resolution.txt
  done
}


mkdir -p "jobs/compare"

rm ../tmp/compare_resolution.txt
touch ../tmp/compare_resolution.txt


cog_rage=200
thre=threshold
# fig_path=../result/resolution/figures/comparsion/${cog_rage}mm
fig_path=../result/resolution/figures/comparsion/${cog_rage}mm


mkdir -p "${fig_path}"

root_path=../result/resolution/comparsion
mkdir -p ${root_path}

# count_data=$(ls ../result/resolution/data/${cog_rage}mm -U1 | wc -l)
# echo -n ${count_data} >> ../tmp/compare_resolution.txt
# echo -n " " >> ../tmp/compare_resolution.txt

echo -n 2 >> ../tmp/compare_resolution.txt
echo -n " " >> ../tmp/compare_resolution.txt

count_sim=$(ls ../result/resolution/data/${cog_rage}mm -U1 | wc -l)
echo -n ${count_sim} >> ../tmp/compare_resolution.txt
echo -n " " >> ../tmp/compare_resolution.txt



paralell_jobs ../result/resolution/data/${cog_rage}mm ${cog_rage}
paralell_jobs ../result/resolution/simulation/${thre}/${cog_rage}mm ${cog_rage}

# cog_rage=20
# paralell_jobs ../result/resolution/data/${cog_rage}mm ${cog_rage}
# paralell_jobs ../result/resolution/simulation/${thre}/${cog_rage}mm ${cog_rage}

# cog_rage=10
# paralell_jobs ../result/resolution/data/${cog_rage}mm ${cog_rage}
# paralell_jobs ../result/resolution/simulation/${thre}/${cog_rage}mm ${cog_rage}

# cog_rage=5
# paralell_jobs ../result/resolution/data/${cog_rage}mm ${cog_rage}
# paralell_jobs ../result/resolution/simulation/${thre}/${cog_rage}mm ${cog_rage}



var=`cat ../tmp/compare_resolution.txt`
# ./threshold_finding ../result/threshold/analyze_${txt}.root ${var} ${fig_path}
./MultiEnergyAnalysis ${root_path}/comparison_${cog_rage}mm.root ${var} ${fig_path}
# sbatch -o jobs/compare/test-%A.out --error="jobs/compare/test-%A.err" execute_paralell.sl ./MultiEnergyAnalysis ${root_path}/comparison_${cog_rage}mm.root ${var} ${fig_path}
