#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
label=${1:-20261002}
[[ $label != *[^a-zA-Z0-9_-]* ]] || exit 1
base="$PROJECT/investigations/result/e_ssa_threshold_review_20260930${label:+/$label}"
printf 'energy_GeV\tdata_mean_GeV\tmc_mean_GeV\tdata_over_mc\tdata_RMS_over_mean\tmc_RMS_over_mean\tquadrature_extra_RMS\tcentral_energy_ratio\tperipheral_energy_ratio\tcentral_hits_ratio\tperipheral_hits_ratio\tcentral_delta_over_mc_total\tperipheral_delta_over_mc_total\tdata_noTemp_shift\tmc_digiMip_shift\tdata_digiMip_shift\tboth_digiMip_ratio\tdata_cleanMip_shift\tmc_cleanMip_shift\tboth_cleanMip_ratio\tdata_fallback_energy_fraction\tmc_fallback_energy_fraction\n' > "$base/impact_summary.tsv"
printf 'energy_GeV\tsample\tcommon_channels\tweighted_raw_residual_ADC\tweighted_ELminusEH_MeV\tHG4000\tLG3000\tLG_hits\n' > "$base/gain_summary.tsv"
for energy in 5 20 100; do
 joint="$base/sensitivity/comparison/20mm/${energy}GeV"
 awk -F '\t' -v beam="$energy" 'NR==1{for(i=1;i<=NF;i++)k[$i]=i;next}$2=="exclude"{s=$1;E[s]=$(k["mean_GeV"]);R[s]=$(k["RMS_over_mean"]);C[s]=$(k["central_GeV"]);N[s]=$(k["hits"]);CN[s]=$(k["central_hits"]);T[s]=$(k["temperature20_GeV"]);M[s]=$(k["digiMip_GeV"]);Q[s]=$(k["cleanMip_GeV"]);V[s]=$(k["nonfallback_GeV"])}END{d="data";m="mc";v=R[d]^2-R[m]^2;printf "%.12g",beam;split("",x);x[1]=E[d];x[2]=E[m];x[3]=E[d]/E[m];x[4]=R[d];x[5]=R[m];x[6]=sqrt(v>0?v:0);x[7]=C[d]/C[m];x[8]=(E[d]-C[d])/(E[m]-C[m]);x[9]=CN[d]/CN[m];x[10]=(N[d]-CN[d])/(N[m]-CN[m]);x[11]=(C[d]-C[m])/E[m];x[12]=((E[d]-C[d])-(E[m]-C[m]))/E[m];x[13]=T[d]/E[d]-1;x[14]=M[m]/E[m]-1;x[15]=M[d]/E[d]-1;x[16]=M[d]/M[m];x[17]=Q[d]/E[d]-1;x[18]=Q[m]/E[m]-1;x[19]=Q[d]/Q[m];x[20]=1-V[d]/E[d];x[21]=1-V[m]/E[m];for(i=1;i<=21;i++)printf "\t%.12g",x[i];printf "\n"}' "$joint/summary.tsv" >> "$base/impact_summary.tsv"
 awk -F '\t' -v beam="$energy" 'NR==1{for(i=1;i<=NF;i++)k[$i]=i;next}{s=$1;id=$2;n1=$(k["band1_n"]);n2=$(k["band2_n"]);n=n1+n2;N[id,s]=n;if(n){R[id,s]=(n1*$(k["band1_residual_mean"])+n2*$(k["band2_residual_mean"]))/n;E[id,s]=(n1*$(k["band1_ELminusEH"])+n2*$(k["band2_ELminusEH"]))/n;}ids[id]=1;H[s]+=$(k["HG4000"]);L[s]+=$(k["LG3000"]);LG[s]+=$(k["LG_hits"])}END{for(id in ids)if(N[id,"data"]>=100&&N[id,"mc"]>=100){n++;w=N[id,"data"]<N[id,"mc"]?N[id,"data"]:N[id,"mc"];W+=w;for(j=0;j<2;j++){s=j?"mc":"data";rs[s]+=w*R[id,s];es[s]+=w*E[id,s]}}for(j=0;j<2;j++){s=j?"mc":"data";printf "%s\t%s\t%d\t%.12g\t%.12g\t%d\t%d\t%d\n",beam,s,n,rs[s]/W,es[s]/W,H[s],L[s],LG[s]}}' "$joint/channels.tsv" >> "$base/gain_summary.tsv"
done
printf 'energy_GeV\tsample\tevents\tmedian_GeV\tsigma68_over_median\n' > "$base/robust_width.tsv"
printf 'energy_GeV\tselection\tdiagonal_variance_difference\toffdiagonal_covariance_difference\toffdiagonal_fraction_of_total_difference\n' > "$base/covariance_impact.tsv"
for energy in 5 20 100; do
 joint="$base/sensitivity/comparison/20mm/${energy}GeV"
 for sample in data mc; do
  awk -F '\t' -v s="$sample" -v e="$energy" 'NR>1&&$1==s&&$5==1&&$6==0{a[++n]=$9}function q(p,x,i){x=p*(n-1)+1;i=int(x);return a[i]+(x-i)*(a[i<n?i+1:i]-a[i])}END{asort(a);printf "%s\t%s\t%d\t%.12g\t%.12g\n",e,s,n,q(.5),(q(.84)-q(.16))/2/q(.5)}' "$joint/events.tsv" >> "$base/robust_width.tsv"
 done
 awk -F '\t' -v e="$energy" 'NR>1&&$2!~/^run/{s=$2;d[s,$1]=$6;o[s,$1]=$7;sel[s]=1}END{for(s in sel){dd=d[s,"data"]-d[s,"mc"];oo=o[s,"data"]-o[s,"mc"];printf "%s\t%s\t%.12g\t%.12g\t%.12g\n",e,s,dd,oo,oo/(dd+oo)}}' "$joint/covariance_summary.tsv" >> "$base/covariance_impact.tsv"
done
