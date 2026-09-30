#!/bin/bash
# Build per-channel response histograms and data/MC comparison plots with ROOT.
set -euo pipefail
export LC_ALL=C
HERE=$(dirname -- "$(realpath -e -- "${BASH_SOURCE[0]}")")
DATA=/megraid01/users/data_beamtest/ECAL_data/analysed/2023
MC=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
PEDESTAL=$MC/Analysis_edit/share/pedestal2023_SPS.root

usage() {
    cat <<EOF
Usage: bash run.sh [options]
  --energy VALUE          Beam energy / default input directory (default: 100)
  --beam NAME             Data beam: auto (default), ps, or sps
                          auto: 0.5,1,2,3,4,5 GeV -> ps; other energies -> sps
  --label NAME            Optional output suffix, e.g. Run105 or smoke
  --mc-tag NAME           MC condition directory (default: threshold)
  --data-decode PATH      Raw_Hit ROOT file or directory
  --data-calib PATH       Calib_Hit ROOT file or directory with matching basenames
  --mc-decode PATH        Raw_Hit ROOT file or directory
  --mc-calib PATH         Calib_Hit ROOT file or directory with matching basenames
  --data-pedestal PATH    Pedestal ROOT containing ChnLevel
  --mc-pedestal PATH      Pedestal ROOT containing ChnLevel
  --channels IDS          Physical CellIDs separated by commas, or all (default)
  --layers RANGES         Layer numbers/ranges (default: 0-29; e.g. 4,9-10)
  --max-events N          Events per sample across files; 0 means all (default)
  --bins N                Bins per axis, 16..1024 (default: 128)
  --hg-max VALUE          Pedestal-subtracted HG ADC upper bound (default: 4200)
  --lg-max VALUE          Pedestal-subtracted LG ADC upper bound (default: 3200)
  --hg-energy-max VALUE   Saved HG energy upper bound [MeV] (default: 40)
  --energy-max VALUE      Saved LG/final energy upper bound [MeV] (default: 400)
  --result-dir PATH       Output directory (default: $HERE/../../result/channel_response/adc_energy)
  --no-plots             Only create histogram ROOT files
  --plots-only           Plot existing histogram ROOT files
  --save-canvases        Also save comparison canvases to ROOT (default: PNG only)
  --scratch-root PATH    Node-local plotting workspace/cache (default: /tmp)
  --dry-run              Validate input paths / print commands; create no outputs
  -h, --help             Show this help

Default data base: $DATA/<beam>/{decode,calib}/e-/<energy>GeV
Default MC base:   $MC/Result_MC/{decode,calib}/e-/e_ssa/<mc-tag>/<energy>GeV
                  Low-energy MC files also live under sps.
Default pedestal: $PEDESTAL
EOF
}
fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
declare -A opt=(
    [energy]=100 [beam]=auto [label]='' [mc-tag]=threshold [channels]=all [layers]=0-29
    [max-events]=0 [bins]=128 [hg-max]=4200 [lg-max]=3200
    [hg-energy-max]=40 [energy-max]=400
    [data-decode]='' [data-calib]='' [mc-decode]='' [mc-calib]=''
    [data-pedestal]="$PEDESTAL" [mc-pedestal]="$PEDESTAL"
    [result-dir]="$HERE/../../result/channel_response/adc_energy"
    [scratch-root]=/tmp
)
no_plots=0 plots_only=0 dry_run=0 save_canvases=0
while (($#)); do
    case $1 in
        -h|--help) usage; exit 0 ;;
        --no-plots) no_plots=1; shift ;;
        --plots-only) plots_only=1; shift ;;
        --save-canvases) save_canvases=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        --*)
            key=${1#--}
            [[ -v opt[$key] ]] || fail "Unknown option: $1"
            (($# >= 2)) || fail "Missing value for $1"
            opt[$key]=$2
            shift 2 ;;
        *) fail "Unexpected argument: $1" ;;
    esac
done
case ${opt[beam]} in
    auto)
        case ${opt[energy]} in
            0.5|1|2|3|4|5) opt[beam]=ps ;;
            *) opt[beam]=sps ;;
        esac ;;
    ps|sps) ;;
    *) fail 'beam must be auto, ps, or sps' ;;
esac
for key in energy mc-tag; do
    [[ ${opt[$key]} =~ ^[A-Za-z0-9_.-]+$ && ${opt[$key]} != . && ${opt[$key]} != .. ]] || fail "Invalid $key"
done
[[ -z ${opt[label]} || ${opt[label]} =~ ^[A-Za-z0-9_-]+$ ]] || fail 'Invalid label'
[[ ${opt[max-events]} =~ ^[0-9]{1,18}$ ]] || fail 'max-events must be a nonnegative integer (up to 18 digits)'
[[ ${opt[bins]} =~ ^[0-9]{1,4}$ ]] || fail 'bins must be 16..1024'
((10#${opt[bins]} >= 16 && 10#${opt[bins]} <= 1024)) || fail 'bins must be 16..1024'
((!no_plots || !plots_only)) || fail '--no-plots and --plots-only cannot be combined'
[[ ${opt[channels]} == all || ${opt[channels]} =~ ^[0-9]+(,[0-9]+)*$ ]] || fail 'channels must be all or comma-separated physical CellIDs'
[[ ${opt[layers]} =~ ^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$ ]] || fail 'Invalid layer selection'
for key in hg-max lg-max hg-energy-max energy-max; do
    [[ ${opt[$key]} =~ ^[+]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]] || fail "Invalid $key"
    minimum=0
    [[ $key != hg-max ]] || minimum=3000
    [[ $key != lg-max ]] || minimum=130
    awk -v value="${opt[$key]}" -v minimum="$minimum" 'BEGIN { exit !(value > minimum && value < 1e308) }' || fail "$key must be finite and > $minimum"
done

tag=${opt[energy]}GeV
[[ -z ${opt[label]} ]] || tag+=_${opt[label]}
out=$(realpath -m -- "${opt[result-dir]}")
opt[result-dir]=$out
for sample in data mc; do
    for stage in decode calib; do
        key=$sample-$stage
        if [[ -z ${opt[$key]} ]]; then
            if [[ $sample == data ]]; then
                opt[$key]=$DATA/${opt[beam]}/$stage/e-/${opt[energy]}GeV
            else
                opt[$key]=$MC/Result_MC/$stage/e-/e_ssa_noHLintercept/${opt[mc-tag]}/${opt[energy]}GeV
                # opt[$key]=$MC/Result_MC/$stage/e-/e_ssa/${opt[mc-tag]}/${opt[energy]}GeV
                # opt[$key]=$MC/Result_MC/$stage/e-/sps/${opt[mc-tag]}/${opt[energy]}GeV
            fi
        fi
    done
done
declare -A outputs=( [data]="$out/data/$tag.root" [mc]="$out/simulation/${opt[mc-tag]}/$tag.root" )
declare -A manifests=( [data]="$out/inputs/${tag}_${opt[mc-tag]}_data.tsv" [mc]="$out/inputs/${tag}_${opt[mc-tag]}_mc.tsv" )
declare -a data_pairs=() mc_pairs=() command=() plot_command=()

collect_pairs() {
    local sample=$1 raw_key=$1-decode cal_key=$1-calib ped_key=$1-pedestal
    local raw=${opt[$raw_key]} cal=${opt[$cal_key]} a b
    local -a files=()
    local -n result=$1_pairs
    [[ -f ${opt[$ped_key]} ]] || fail "Missing pedestal file: ${opt[$ped_key]}"
    opt[$ped_key]=$(realpath -e -- "${opt[$ped_key]}")
    [[ -e $raw && -e $cal ]] || fail "Input does not exist: $raw or $cal"
    if [[ -f $raw ]]; then
        files=("$raw")
    elif [[ -d $raw && -d $cal ]]; then
        shopt -s nullglob
        files=("$raw"/*.root)
        shopt -u nullglob
    else
        fail 'A decode directory requires a calibration directory'
    fi
    ((${#files[@]})) || fail "No ROOT files in $raw"
    for a in "${files[@]}"; do
        [[ -f $a ]] || fail "Not a ROOT file: $a"
        b=$cal
        [[ ! -d $cal ]] || b=$cal/${a##*/}
        [[ -f $b ]] || fail "Missing paired calibration file: $b"
        a=$(realpath -e -- "$a")
        b=$(realpath -e -- "$b")
        [[ $a$b != *$'\t'* && $a$b != *$'\n'* ]] || fail 'Input paths must not contain tabs or newlines'
        result+=("$a"$'\t'"$b")
    done
    printf '%s: %s file pair(s)\n  decode: %s\n  calib:  %s\n' "$sample" "${#result[@]}" "$raw" "$cal"
}
fill_command() {
    local sample=$1 key
    command=("$HERE/ChannelResponse" --manifest "${manifests[$sample]}"
        --pedestal "${opt[$sample-pedestal]}" --output "${outputs[$sample]}" --sample "$sample")
    for key in channels layers max-events bins hg-max lg-max hg-energy-max energy-max; do
        command+=("--$key" "${opt[$key]}")
    done
}
print_command() { printf '%q ' "$@"; printf '\n'; }
for sample in data mc; do
    if ((plots_only)); then
        [[ -f ${outputs[$sample]} ]] || fail "Missing histogram file: ${outputs[$sample]}"
    else
        collect_pairs "$sample"
        fill_command "$sample"
        print_command "${command[@]}"
    fi
done
if ((!no_plots)); then
    canvas_output=$out/comparison_${opt[mc-tag]}_$tag.root
    if [[ ${opt[layers]} =~ ^[0-9]+$ ]];then
        canvas_output=$out/comparisons/${opt[mc-tag]}/$tag/layer${opt[layers]}.root
    fi
    plot_command=(bash "$HERE/plot.sh" --data "${outputs[data]}" --mc "${outputs[mc]}"
        --output "$canvas_output"
        --figure-dir "$out/figures/comparison/${opt[mc-tag]}/$tag"
        --channels "${opt[channels]}" --energy "${opt[energy]}" --layers "${opt[layers]}"
        --save-canvases "$save_canvases" --scratch-root "${opt[scratch-root]}")
    print_command "${plot_command[@]}"
fi
((!dry_run)) || exit 0
build_targets=()
((plots_only)) || build_targets+=(ChannelResponse)
((no_plots)) || build_targets+=(PlotChannelResponse)
make -C "$HERE" "${build_targets[@]}"
if ((!plots_only)); then
    mkdir -p -- "$out/inputs"
    for sample in data mc; do
        declare -n pairs="${sample}_pairs"
        printf '%s\n' "${pairs[@]}" > "${manifests[$sample]}"
        unset -n pairs
        fill_command "$sample"
        "${command[@]}"
    done
    # Shell-quoted values/commands are readable and reproducible without Python.
    {
        printf '# Channel response configuration (shell-quoted values)\n'
        for key in energy beam label mc-tag channels layers max-events bins hg-max lg-max hg-energy-max energy-max result-dir scratch-root data-decode data-calib mc-decode mc-calib data-pedestal mc-pedestal; do
            printf '%s=%q\n' "$key" "${opt[$key]}"
        done
        printf 'no-plots=%s\n' "$no_plots"
        printf 'save-canvases=%s\n' "$save_canvases"
        for sample in data mc; do
            fill_command "$sample"
            print_command "${command[@]}"
        done
        ((no_plots)) || print_command "${plot_command[@]}"
    } > "$out/inputs/${tag}_${opt[mc-tag]}_configuration.txt"
fi
if ((!no_plots)); then
    "${plot_command[@]}"
fi
