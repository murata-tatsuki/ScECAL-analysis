#!/bin/bash
set -euo pipefail
out=${1:?Usage: verify_results.sh OUTPUT_DIRECTORY}
[[ $(awk 'END{print NR-1}' "$out/event_populations.tsv") -eq 57 ]]
for energy in 0.5 1 5 30 100 200 250; do
 dir=$out/${energy}GeV
 [[ -s $dir/response.root && -s $dir/figures/energy_at_fixed_nhit.png ]]
 awk -F '\t' 'NR==2{exit !($1<1e-9&&$2<1e-9&&$3==0&&$6==0)}' "$dir/validation.tsv"
 awk -F '\t' 'NR>1{err=$5-$6-$11;if(err<0)err=-err;if(err>1e-8)bad=1;if($1=="mc"&&$2=="all"&&$3!=20000)bad=1} END{exit bad}' "$dir/event_metrics.tsv"
 awk -F '\t' 'FNR==NR{if($2=="shower_cog20"){events[$1]=$3;expected[$1]=$7};next} FNR>1{sums[$1]+=$7} END{for(s in sums){err=sums[s]/events[s]-expected[s];if(err<0)err=-err;if(err>1e-8)exit 1}}' "$dir/event_metrics.tsv" "$dir/channels.tsv"
done
printf '%s\n' 'Verified 19-energy survey, 7 detailed samples, signed/positive bookkeeping, channel sums, and original paired-MC totals.'
