#!/bin/bash
set -euo pipefail
out=${1:?Usage: summarize.sh OUTPUT_DIRECTORY INPUT_LABEL_DIRECTORY}
base=${2:?Input directory required}
awk -F '\t' '
FNR==1 {next}
$2=="shower_cog20" {split(FILENAME,a,"/");e=a[length(a)-1];sub("GeV","",e);for(i=3;i<=24;i++)v[e,$1,i]=$i;energies[e]=1}
END{print "beam_GeV\tdata_events\tmc_events\tgood_energy_data_over_mc\tgood_energy_data_over_zero\tcentral_energy_data_over_mc\tcentral_energy_data_over_zero\tgood_hits_data_over_mc\tcentral_hits_data_over_mc\tperipheral_energy_data_over_mc\tdata_fallback_energy_fraction\tmc_fallback_energy_fraction";for(e in energies)printf "%g\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\n",e,v[e,"data",3],v[e,"mc",3],v[e,"data",7]/v[e,"mc",7],v[e,"data",7]/v[e,"mc",8],v[e,"data",9]/v[e,"mc",9],v[e,"data",9]/v[e,"mc",10],v[e,"data",13]/v[e,"mc",13],v[e,"data",14]/v[e,"mc",14],(v[e,"data",7]-v[e,"data",9])/(v[e,"mc",7]-v[e,"mc",9]),v[e,"data",17]/v[e,"data",7],v[e,"mc",17]/v[e,"mc",7]}' "$out"/*GeV/event_metrics.tsv > "$out/controlled_response_ratios.tsv"
awk -F '\t' 'FNR==NR {if(FNR>1)ok[$2]=($5==0&&$3>0);next} FNR>1&&ok[$1]&&$8>=100&&$11>=100 {w=($8<$11?$8:$11);n++;sw+=w;d+=w*$9;m+=w*$12;z+=w*$15;dr+=w*$9*$9;mr+=w*$12*$12;zr+=w*$15*$15} END {print "channels\tweighted_data_residual_MeV\tweighted_mc_residual_MeV\tweighted_zero_residual_MeV\tdata_channel_rms\tmc_channel_rms\tzero_channel_rms"; print n"\t"d/sw"\t"m/sw"\t"z/sw"\t"sqrt(dr/sw)"\t"sqrt(mr/sw)"\t"sqrt(zr/sw)}' "$out/100GeV/channels.tsv" "$base/100GeV/channels.tsv" > "$out/gain_consistency_100GeV.tsv"
