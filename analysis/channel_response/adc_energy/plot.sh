#!/bin/bash
# Stage immutable histogram ROOTs once per node/version; render locally and
# publish only completed images. Does not need original raw/calibrated trees.
set -euo pipefail
HERE=$(dirname -- "$(realpath -e -- "${BASH_SOURCE[0]}")")
declare -A opt=([data]='' [mc]='' [output]='' [figure-dir]='' [channels]=all
 [energy]=100 [layers]=0-29 [save-canvases]=0 [scratch-root]=/tmp)
while (($#));do
 key=${1#--};[[ $1 == --* && -v opt[$key] && $# -ge 2 ]]||{ echo "Invalid plot argument: $1" >&2;exit 1; }
 opt[$key]=$2;shift 2
done
[[ -f ${opt[data]} && -f ${opt[mc]} && -n ${opt[figure-dir]} ]]||{ echo 'Missing plotting inputs/output' >&2;exit 1; }
[[ ${opt[layers]} =~ ^[0-9,-]+$ ]]||exit 1
command -v rsync >/dev/null || { echo 'rsync is required for publishing PNGs' >&2;exit 1; }
mkdir -p -- "${opt[scratch-root]}" "${opt[figure-dir]}"
scratch=$(mktemp -d "${opt[scratch-root]%/}/channel_response_plot.XXXXXX")
lock_key=$(printf '%s' "${opt[layers]}" | sha256sum | cut -d' ' -f1)
plot_lock=${opt[figure-dir]}/.plot_${lock_key}.lock
plot_lock_owned=0
cleanup() {
 local status=$?
 # Remove only our own lock, while its descriptor is still locked.
 if ((plot_lock_owned)) && [[ $plot_lock -ef /proc/$$/fd/8 ]];then
  rm -f -- "$plot_lock" || status=1
 fi
 rm -rf -- "$scratch" || status=1
 return "$status"
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
while :;do
 exec 8>"$plot_lock"
 flock -n 8 || { echo 'Another plot worker is writing these layers' >&2;exit 1; }
 # A finishing worker may have unlinked the file after we opened it.
 # Retry against the current pathname before allowing any output writes.
 if [[ $plot_lock -ef /proc/$$/fd/8 ]];then
  plot_lock_owned=1
  break
 fi
 exec 8>&-
done
cache_dir=${opt[scratch-root]%/}/channel_response_cache_$UID
mkdir -p -- "$cache_dir"
stage_input() {
 local source fingerprint key target
 source=$(realpath -e -- "$1")
 fingerprint=$(stat -Lc '%s:%y:%i' -- "$source")
 key=$(printf '%s\n%s\n' "$source" "$fingerprint" | sha256sum | cut -d' ' -f1)
 target=$cache_dir/$key.root
 (
  flock 9
  if [[ ! -s $target ]];then
   trap 'rm -f -- "$target.tmp.$$"' EXIT
   cp --reflink=auto -- "$source" "$target.tmp.$$"
   [[ $(stat -Lc '%s:%y:%i' -- "$source") == "$fingerprint" ]] || { echo "Input changed while copying: $source" >&2;exit 1; }
   mv -- "$target.tmp.$$" "$target"
   printf 'Cached %s\n' "$source" >&2
  fi
 ) 9>"$cache_dir/$key.lock"
 printf '%s\n' "$target"
}
started=$SECONDS
data_local=$(stage_input "${opt[data]}")
mc_local=$(stage_input "${opt[mc]}")
printf 'Input staging: %s seconds; layers=%s\n' "$((SECONDS-started))" "${opt[layers]}"
started=$SECONDS
"$HERE/PlotChannelResponse" --data "$data_local" --mc "$mc_local" \
 --data-origin "${opt[data]}" --mc-origin "${opt[mc]}" \
 --output "$scratch/comparison.root" --figure-dir "$scratch/figures" \
 --channels "${opt[channels]}" --energy "${opt[energy]}" --layers "${opt[layers]}" \
 --save-canvases "${opt[save-canvases]}"
printf 'Local rendering: %s seconds\n' "$((SECONDS-started))"
started=$SECONDS
if [[ -d $scratch/figures ]];then
 rsync -a --delay-updates -- "$scratch/figures/" "${opt[figure-dir]}/"
fi
if [[ -f $scratch/comparison.root ]];then
 [[ -n ${opt[output]} ]]||{ echo 'Canvas output path missing' >&2;exit 1; }
 mkdir -p -- "$(dirname -- "${opt[output]}")"
 cp -- "$scratch/comparison.root" "${opt[output]}.tmp.$$"
 mv -- "${opt[output]}.tmp.$$" "${opt[output]}"
fi
printf 'Published figures: %s seconds; %s\n' "$((SECONDS-started))" "${opt[figure-dir]}"
