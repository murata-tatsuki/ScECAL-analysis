"""Small ROOT plots and channel-balanced summaries; no additional fit cuts."""
import csv
import math
from pathlib import Path

import numpy as np
from mip_tools import write_tsv


CHANNELS_PER_LAYER = 5 * 36 + 30


def channel_index(cellid):
    """Zero-based connected-channel index, ordered by layer/chip/channel."""
    layer, remainder = divmod(int(cellid), 100000)
    chip, channel = divmod(remainder, 10000)
    if not (0 <= layer < 32 and 0 <= chip < 6
            and 0 <= channel < (30 if chip == 5 else 36)):
        raise ValueError(f'Invalid connected physical CellID {cellid}')
    return layer * CHANNELS_PER_LAYER + chip * 36 + channel


def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def report(ROOT, out, channels):
    groups = []
    for p in (10, 15):
        for digi in ('individual', 'fallback', 'runtime_default'):
            for calibration in ('individual', 'fallback'):
                selected = [r for r in channels if r['pitch_um'] == p and r['digi_source'] == digi
                            and r['mc_legacy_source'] == calibration and r['mc_adc_fit_valid']]
                v = [r['mc_adc_mpv_over_input'] for r in selected if math.isfinite(r['mc_adc_mpv_over_input'])]
                paired = [r for r in selected if r['data_fit_valid'] and r['mc_fit_valid']]
                row = dict(pitch_um=p, digi_source=digi, mc_legacy_source=calibration,
                           n_valid_mc_channels=len(v), median_mc_adc_over_input=float(np.median(v)) if v else math.nan,
                           q16=float(np.quantile(v, .16)) if v else math.nan,
                           q84=float(np.quantile(v, .84)) if v else math.nan,
                           n_paired_channels=len(paired))
                for policy in ('legacy', 'clean'):
                    a = [r[policy+'_data_over_mc'] for r in paired if math.isfinite(r[policy+'_data_over_mc'])]
                    row[policy+'_median_data_over_mc'] = float(np.median(a)) if a else math.nan
                groups.append(row)
    write_tsv(out/'group_summary.tsv', groups)
    calibrated = []
    for sample in ('data', 'mc'):
        for policy in ('legacy', 'clean'):
            for p in (10, 15):
                for source in ('individual', 'fallback'):
                    values = [r[f'{sample}_{policy}_mip_mpv'] for r in channels
                              if r['pitch_um'] == p and r[f'{sample}_{policy}_source'] == source
                              and r[sample+'_fit_valid'] and math.isfinite(r[f'{sample}_{policy}_mip_mpv'])]
                    calibrated.append(dict(sample=sample, policy=policy, pitch_um=p, source=source,
                                           n_valid_channels=len(values),
                                           median_mip_mpv=float(np.median(values)) if values else math.nan))
    write_tsv(out/'calibrated_group_summary.tsv', calibrated)
    ROOT.gStyle.SetOptStat(0)
    if channels:
        # Keep detector positions stable even for channels without a valid fit.
        n_layers = max(30, max(int(r['cellid']) // 100000 for r in channels) + 1)
        axis_rows = [dict(channel_index=channel_index(layer*100000 + chip*10000 + channel),
                          cellid=layer*100000 + chip*10000 + channel,
                          layer=layer, chip=chip, channel=channel)
                     for layer in range(n_layers) for chip in range(6)
                     for channel in range(30 if chip == 5 else 36)]
        write_tsv(out/'channel_axis.tsv', axis_rows)
        from visual_report import draw
        draw(ROOT, out, channels)
    # Pre/post policy comparisons use exactly the same event and hit sets.
    energies = rows(out/'energy_summary.tsv')
    beam = rows(out/'beam_truth.tsv')
    coverage = rows(out/'matching_coverage.tsv')
    text = [
        'Muon MIP closure diagnostics',
        'Plots: response distributions, layer/channel map, and labeled paired channels. See visualization.json.',
        'MPV = Landau MPV parameter of a Landau x Gaussian fit; convolution mode is stored separately.',
        'Constants are current-source candidates. saved_calibration_check.tsv tests the saved HG inverse conversion.',
        'Fits use symmetric overlap weights in run/channel/strip position/reconstructed slopes.',
        'Energy totals use unweighted full connected calibrated hit sets of selected muon events, frozen across policies.',
        'An undefined old constant or failed ADC join makes the entire event comparison invalid; subtotals are diagnostic only.',
        'Truth fits are conditional on the same selected detected hits, not an unbiased generated-track dE/dx sample.',
        'The plateau fit uses external threshold constants. MPV below the plateau or failed fits are excluded from summaries.',
        'Input-MIP, temperature law, geometry/track, threshold and fit-window uncertainties are not included in statistical errors.',
        'primary_pe_scale_trial_only = 1/(MC ADC MPV / input MIP): a starting trial, requiring redigitization with crosstalk ON and rechecking MPV.',
        'A high ADC ratio alone does not prove crosstalk double counting: also inspect truth MPV / 0.305 and production provenance.',
        '', 'Saved MC primary direction audit:',
    ]
    for r in beam:
        text.append(f"  Run{r['run']}: {r['n_primaries']} primaries, non-(0,0,1)={r['n_non_axial']}, z=[{r['pz_min']},{r['pz_max']}] mm")
    text += ['', 'Geometric overlap:']+[str(r) for r in coverage]
    text += ['', 'Energy policy changes (all hits of each valid event):']+[str(r) for r in energies]
    text += ['', 'Read STATUS.json and inputs/excluded_runs.tsv for incomplete/recovered input coverage.']
    (out/'report.txt').write_text('\n'.join(text)+'\n')
