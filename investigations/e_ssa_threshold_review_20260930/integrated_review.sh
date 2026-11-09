#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
label=${1:-20261012}
[[ $label != *[^a-zA-Z0-9_-]* && -n $label ]] || { echo 'Invalid label' >&2; exit 1; }
base="$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$label"
out="$base/integrated"
[[ ! -e $out/summary.tsv ]] || { echo "Output exists: $out" >&2; exit 1; }
mkdir -p "$out/inputs"
exec > >(tee "$out/run.log") 2>&1
printf 'Started: %s\nCommand: bash %q %q\n' "$(date --iso-8601=seconds)" "$HERE/integrated_review.sh" "$label"
MC=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
DATA=/home/murata_t/data_beamtest/ECAL_data/analysed/2023
params="$MC/Analysis_edit/share"
args=(--hl "$params/all_hl_electron2023_hlratio_v3.root" --mip "$params/all_auto_muon_v4_trackfit.root" --pedestal "$params/pedestal2023_SPS.root")
make -C "$HERE" PrepareReference IntegratedReview
shopt -s nullglob
for energy in 5 20 100; do
 beam=sps; [[ $energy != 5 ]] || beam=ps
 dir="$base/${energy}GeV/baseline"
 [[ ! -e $dir/study.root ]] || { echo "Baseline exists: $dir" >&2; exit 1; }
 mkdir -p "$dir/inputs"
 files=("$DATA/$beam/decode/e-/${energy}GeV/"*.root)
 ((${#files[@]})) || { echo 'No data input' >&2; exit 1; }
 for raw in "${files[@]}"; do cal="$DATA/$beam/calib/e-/${energy}GeV/${raw##*/}"; [[ -f $cal ]]; printf '%s\t%s\n' "$raw" "$cal"; done > "$dir/inputs/data.tsv"
 files=("$MC/Result_MC/decode/e-/e_ssa/threshold/${energy}GeV/"*.root)
 ((${#files[@]}>=6)) || { echo 'Fewer than 6 MC files' >&2; exit 1; }
 for i in 0 1 2 3 4 5; do idx=$((i*(${#files[@]}-1)/5)); raw=${files[$idx]}; cal="$MC/Result_MC/calib/e-/e_ssa/threshold/${energy}GeV/${raw##*/}"; [[ -f $cal ]]; printf '%s\t%s\n' "$raw" "$cal"; done > "$dir/inputs/mc.tsv"
 while IFS=$'\t' read -r raw cal; do stat -c '%y\t%s\t%n' "$raw" "$cal"; done < <(cat "$dir/inputs/"{data,mc}.tsv) > "$dir/inputs/file_stats_before.tsv"
 "$HERE/PrepareReference" --mc-manifest "$dir/inputs/mc.tsv" "${args[@]}" --energy "$energy" --out "$dir" > "$dir/run.log" 2>&1
 printf '%sGeV reference ready; %s MC files available\n' "$energy" "${#files[@]}"
 sha256sum "$PROJECT/analysis/result/channel_response/hg_lg_tail/threshold/e_ssa/${energy}GeV/"{tail_study.root,inputs/data.tsv} > "$dir/inputs/data_tail_checksums.txt"
 for sample in data simulation/threshold; do for cog in 200 20 10 5; do path="$PROJECT/analysis/result/resolution/$sample/${cog}mm/${energy}GeV_${cog}mm.root"; stat -c '%y\t%s\t%n' "$path"; done; done >> "$out/inputs/production_resolution_stats.tsv"
done
cd "$HERE"
"$HERE/IntegratedReview" "$base" "$out"
for energy in 5 20 100; do
 dir="$base/${energy}GeV/baseline"
 while IFS=$'\t' read -r raw cal; do stat -c '%y\t%s\t%n' "$raw" "$cal"; done < <(cat "$dir/inputs/"{data,mc}.tsv) > "$dir/inputs/file_stats_after.tsv"
 cmp "$dir/inputs/file_stats_before.tsv" "$dir/inputs/file_stats_after.tsv"
done
sha256sum "$HERE/"{IntegratedReview.cc,CurrentParameters.hh,integrated_review.sh,ResolutionImpact.cc} "$PROJECT/analysis/impact_studies/StudyCommon.hh" "$PROJECT/analysis/resolution/HLTailSelection.hh" "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$params/"{all_hl_electron2023_hlratio_v3.root,all_auto_muon_v4_trackfit.root,pedestal2023_SPS.root} > "$out/inputs/checksums.txt"
printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
# Derived comparisons use the fixed nominal shower sample; no event is reselected.
awk 'BEGIN{FS=OFS="\t"} NR>1&&$3=="cog20_energy_window"&&$4=="baseline"{e=$1;s=$2;m[e,s]=$6;c[e,s]=$11;n[e,s]=$12;cn[e,s]=$13;r[e,s]=$7;g[e,s]=$9}END{print "energy_GeV","data_mean_GeV","mc_mean_GeV","data_mc_mean_ratio","central_energy_ratio","peripheral_energy_ratio","central_hit_ratio","peripheral_hit_ratio","data_RMS","mc_RMS","data_gaussian","mc_gaussian";for(e=5;e<=100;e+=(e==5?15:80))print e,m[e,"data"],m[e,"mc"],m[e,"data"]/m[e,"mc"],c[e,"data"]/c[e,"mc"],(m[e,"data"]-c[e,"data"])/(m[e,"mc"]-c[e,"mc"]),cn[e,"data"]/cn[e,"mc"],(n[e,"data"]-cn[e,"data"])/(n[e,"mc"]-cn[e,"mc"]),r[e,"data"],r[e,"mc"],g[e,"data"],g[e,"mc"]}' "$out/summary.tsv" > "$out/response_comparison.tsv"
awk 'BEGIN{FS=OFS="\t"} NR>1{e=$1;s=$2;id=$3;n[e,s,id]=$4;v[e,s,id]=$5;h[e,s]+=$8;l[e,s]+=$9}END{print "energy_GeV","common_channels","data_EL_minus_EH_MeV","mc_EL_minus_EH_MeV","data_HG4000_hits","mc_HG4000_hits","data_LG3000_hits","mc_LG3000_hits";for(k in n){split(k,a,SUBSEP);e=a[1];id=a[3];if(a[2]!="data"||n[e,"data",id]<100||n[e,"mc",id]<100)continue;w=n[e,"data",id]<n[e,"mc",id]?n[e,"data",id]:n[e,"mc",id];cnt[e]++;sw[e]+=w;d[e]+=w*v[e,"data",id];m[e]+=w*v[e,"mc",id]}for(e=5;e<=100;e+=(e==5?15:80))print e,cnt[e],d[e]/sw[e],m[e]/sw[e],h[e,"data"],h[e,"mc"],l[e,"data"],l[e,"mc"]}' "$out/gain_channels.tsv" > "$out/gain_overlap.tsv"
awk 'BEGIN{FS=OFS="\t"} NR>1&&$3=="cog20"{if($4=="baseline")base[$1,$2]=$5;if($4=="back5_le_5percent")cut[$1,$2]=$5} END{print "energy_GeV","data_events","mc_events","data_back5_gt5percent_fraction","mc_back5_gt5percent_fraction";for(e=5;e<=100;e+=(e==5?15:80))print e,base[e,"data"],base[e,"mc"],1-cut[e,"data"]/base[e,"data"],1-cut[e,"mc"]/base[e,"mc"]}' "$out/summary.tsv" > "$out/leak_proxy_rates.tsv"
sha256sum "$MC/Analysis_edit/src/"{MCDigi,Extract}.cxx /home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO/src/{Calibration,EBUdecode,StripSplit}.cxx > "$out/inputs/producer_checksums.txt"
