#!/usr/bin/env bash
set -euo pipefail
PROJECT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../../../../../.." && pwd -P)
MU=$PROJECT/analysis/result/channel_response/muon_mip_response/threshold_MPVwoCT/inputs/calibration_parameters.tsv
REF=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/20261013/raw_adc/inputs/effective_parameters.tsv
awk -F '\t' 'NR==FNR{if(FNR>1){m[$1]=$3;f[$1]=$9}next}FNR>1&&$2=="data"{if(!($1 in m)){print "missing",$1;exit 1}d=$4-m[$1];if(d<0)d=-d;if(d>max)max=d;bad+=($7!=f[$1]);n++}END{print "checked_channels\tmax_MIP_ADC_difference\tfallback_flag_mismatches";print n "\t" max "\t" bad;if(n!=6300||max>1e-6||bad)exit 1}' "$REF" "$MU"
