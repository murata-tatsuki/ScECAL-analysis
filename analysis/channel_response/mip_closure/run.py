#!/usr/bin/env python3
"""100 GeV muon ADC/MIP closure and a symmetric calibration-policy audit."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
from mip_tools import (write_tsv, pitch, read_mip, read_pedestal, read_threshold,
                       common_weights, load_fit_function, fit_spectrum)

HERE = Path(__file__).resolve().parent
MC = Path('/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new')
DATA = Path('/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps')


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=HERE/'../../result/channel_response/mip_closure/default')
    p.add_argument('--runs', default='all', help='all or comma-separated run numbers')
    p.add_argument('--data-base', type=Path, default=DATA)
    p.add_argument('--mc-base', type=Path, default=MC/'Result_MC')
    for sample in ('data', 'mc'):
        p.add_argument(f'--{sample}-mip', type=Path, default=MC/'Analysis_edit/share/all_auto_muon_v4_trackfit.root')
        p.add_argument(f'--{sample}-pedestal', type=Path, default=MC/'Analysis_edit/share/pedestal2023_SPS.root')
        p.add_argument(f'--{sample}-threshold', type=Path,
                       default=Path('/home/murata_t/ScECAL_BeamTest/analysis/result/threshold/threshold.root'))
    p.add_argument('--digi-mip', type=Path, default=MC/'Analysis_edit/share/all_auto_muon_v4_trackfit.root')
    p.add_argument('--max-tracks', type=int, default=0, help='evenly spaced tracks/file; 0=all')
    p.add_argument('--max-angle', type=float, default=.05, help='absolute projected angle, rad')
    p.add_argument('--max-track-rms', type=float, default=3., help='sqrt(saved objective/NDF), mm, not statistical chi2')
    p.add_argument('--long-bin', type=float, default=5.)
    p.add_argument('--transverse-bin', type=float, default=1.5)
    p.add_argument('--slope-bin', type=float, default=.01)
    p.add_argument('--min-bin', type=int, default=3)
    p.add_argument('--min-fit-hits', type=int, default=100, help='minimum effective entries inside fit interval')
    p.add_argument('--plateau-sigmas', type=float, default=3., help='common lower fit edge above threshold + N sigma')
    p.add_argument('--channels', default='all', help='restrict fitted CellIDs; energy audit still uses all hits')
    p.add_argument('--allow-recovered', action='store_true', help='diagnostic only: include recovered ROOT trees and flag results')
    p.add_argument('--dry-run', action='store_true', help='inspect paths and ROOT metadata only')
    p.add_argument('--no-fit', action='store_true', help='collect and audit without expensive channel fits')
    a = p.parse_args()
    if not all(math.isfinite(v) for v in (a.max_angle, a.max_track_rms, a.long_bin,
                                          a.transverse_bin, a.slope_bin, a.plateau_sigmas)):
        p.error('Selection and bin values must be finite')
    if a.max_tracks < 0 or min(a.max_angle, a.max_track_rms, a.long_bin, a.transverse_bin, a.slope_bin) <= 0:
        p.error('Invalid selection or bin size')
    if a.min_bin < 1 or a.min_fit_hits < 20 or a.plateau_sigmas < 0:
        p.error('Invalid fit/overlap settings')
    return a


def root_info(ROOT, path, tree_name):
    if not path.is_file():
        return dict(path=str(path), tree=tree_name, entries=0, recovered=False, status='missing', size=0, mtime_ns=0)
    f = ROOT.TFile.Open(str(path))
    t = f.Get(tree_name) if f and not f.IsZombie() else None
    row = dict(path=str(path), tree=tree_name, entries=int(t.GetEntries()) if t else 0,
               recovered=bool(f and f.TestBit(ROOT.TFile.kRecovered)), status='ok' if t else 'unreadable',
               size=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns)
    if f:
        f.Close()
    return row


def discover(ROOT, a):
    data_files = sorted((a.data_base/'decode/mu-/100GeV').glob('ECAL_Run*.root'))
    wanted = None if a.runs == 'all' else set(map(int, a.runs.split(',')))
    found, manifest, metadata, excluded = set(), [], [], []
    for raw in data_files:
        run = int(re.search(r'Run(\d+)', raw.name).group(1))
        if wanted is not None and run not in wanted:
            continue
        if run in found:
            raise ValueError(f'Multiple data files for run {run}; provide a directory with one file/run')
        found.add(run)
        paths = [[raw, a.data_base/'calib/mu-/100GeV'/raw.name,
                  a.data_base/'trackFit3D/mu-/100GeV'/raw.name],
                 [a.mc_base/stage/'mu-/mu_track/threshold/100GeV'/f'ECAL_Run{run}.root'
                  for stage in ('decode', 'calib', 'trackFit3D')]]
        paths[1].append(a.mc_base/'generate/mu-/mu_track/100GeV'/f'ECAL_Run{run}.root')
        infos = []
        for sample, files in enumerate(paths):
            for path, tree in zip(files, ('Raw_Hit', 'Calib_Hit', 'T_Event', 'MC_Truth')):
                info = root_info(ROOT, path, tree)
                infos.append(info)
                metadata.append(dict(sample=sample, run=run, **info))
        reasons = sorted({i['status'] for i in infos if i['status'] != 'ok'})
        if not a.allow_recovered and any(i['recovered'] for i in infos):
            reasons.append('recovered_root')
        if reasons:
            excluded.append(dict(run=run, reason=','.join(reasons)))
            continue
        for sample, files in enumerate(paths):
            manifest.append([str(sample), str(run), *map(str, files), *(['-'] if sample == 0 else [])])
    if wanted and wanted-found:
        raise ValueError(f'Data runs not found: {sorted(wanted-found)}')
    return manifest, metadata, excluded


def fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def snapshot(a, out):
    inputs = out/'inputs'; inputs.mkdir()
    (inputs/'configuration.json').write_text(json.dumps(vars(a), default=str, indent=2)+'\n')
    refs = [*HERE.glob('*.py'), *HERE.glob('*.cc'), *HERE.glob('*.sh'), HERE/'Makefile',
            MC/'CEPCScECAL/src/PrimaryGeneratorAction.cc', MC/'CEPCScECAL/generate/generator_sps_mu.sh',
            MC/'Analysis_edit/src/MCDigi.cxx', MC/'Analysis_edit/src/Extract.cxx',
            Path('/megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/src/Calibration.cxx'),
            Path('/megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/src/TrackFit3D.cxx'),
            Path('/home/murata_t/ScECAL_BeamTest/analysis/trackFit/beamsize/beamsize_track.cc')]
    records = []
    for i, path in enumerate(refs):
        if path.is_file():
            shutil.copyfile(path, inputs/f'{i}_{path.name}')
            records.append(dict(kind='current_source_not_production_proof', path=str(path), sha256=fingerprint(path)))
    for name, path in vars(a).items():
        if name.endswith(('_mip', '_pedestal', '_threshold')):
            records.append(dict(kind=name, path=str(path), sha256=fingerprint(path)))
    write_tsv(inputs/'checksums.tsv', records)


def lookup(rows, key, ids):
    table = {r['cellid']: r[key] for r in rows}
    unique, inverse = np.unique(ids, return_inverse=True)
    return np.asarray([table.get(int(c), math.nan) for c in unique])[inverse]


def energy_audit(a, out, h, e, policies, pedestals):
    """Freeze the full connected calibrated-hit set for each selected muon event.

    HG is recalculated from ADC. LG retains the saved numerator, provided the
    saved HG branch validates the hypothesized old MIP constant for that hit.
    Invalid old denominators/join failures invalidate the whole event's comparison.
    """
    size = len(h['cellid'])
    old, new, ped = [np.full(size, np.nan) for _ in range(3)]
    for s in (0, 1):
        m = h['sample'] == s
        old[m] = lookup(policies[s], 'legacy_adc', h['cellid'][m])
        new[m] = lookup(policies[s], 'clean_adc', h['cellid'][m])
        unique, inverse = np.unique(h['cellid'][m], return_inverse=True)
        ped[m] = np.asarray([pedestals[s].get(int(c), np.nan) for c in unique])[inverse]
    layer = h['cellid']//100000
    coefficient = np.where((layer >= 4)&(layer <= 27), 1.6/135, 3.5/230)
    temp_factor = 1-(h['temperature']-20)*coefficient
    adc = h['raw_hg']-ped
    numerator = adc*temp_factor*.305
    with np.errstate(divide='ignore', invalid='ignore'):
        old_hg = numerator/old
        new_hg = numerator/new
        relative = (old_hg-h['saved_hg'])/np.maximum(np.abs(h['saved_hg']), 1.e-8)
    closure = np.isfinite(relative)&(np.abs(relative) < 2.e-4)
    known = (h['raw_ok'] == 1)&np.isfinite(old_hg)&np.isfinite(new_hg)&(new > 0)&(temp_factor > 0)
    known &= np.isfinite(h['saved_energy'])&closure
    high = h['raw_hg'] < 2600
    with np.errstate(divide='ignore', invalid='ignore'):
        baseline = np.where(high, old_hg, h['saved_energy'])
        corrected = np.where(high, new_hg, h['saved_energy']*old/new)
    known &= np.isfinite(baseline)&np.isfinite(corrected)
    event_id = h['event_index'].astype(int)
    n_events = len(e['event_index'])
    if not np.array_equal(e['event_index'], np.arange(n_events)):
        raise ValueError('Noncontiguous event index')
    nh = np.bincount(event_id, minlength=n_events)
    invalid = np.bincount(event_id, weights=~known, minlength=n_events).astype(int)
    mismatch = np.bincount(event_id, weights=~closure, minlength=n_events).astype(int)
    before = np.bincount(event_id, weights=np.where(known, baseline, 0), minlength=n_events)
    after = np.bincount(event_id, weights=np.where(known, corrected, 0), minlength=n_events)
    saved = np.bincount(event_id, weights=np.where(np.isfinite(h['saved_energy']), h['saved_energy'], 0), minlength=n_events)
    if not np.array_equal(nh, e['n_calib']):
        raise ValueError('Frozen hit count does not equal calibrated connected hit count')
    rows = []
    for i in range(n_events):
        valid = invalid[i] == 0 and nh[i] > 0
        rows.append(dict(sample=int(e['sample'][i]), run=int(e['run'][i]), event=int(e['event'][i]),
                         event_time=int(e['event_time'][i]), n_hits=int(nh[i]),
                         n_invalid=int(invalid[i]), n_saved_hg_mismatch=int(mismatch[i]), valid=valid,
                         legacy_energy_mev=float(before[i]) if valid else math.nan,
                         clean_energy_mev=float(after[i]) if valid else math.nan,
                         delta_mev=float(after[i]-before[i]) if valid else math.nan,
                         ratio=float(after[i]/before[i]) if valid and before[i] != 0 else math.nan,
                         valid_hit_subtotal_legacy=float(before[i]), valid_hit_subtotal_clean=float(after[i]),
                         saved_energy_finite_subtotal=float(saved[i])))
    write_tsv(out/'event_energy.tsv', rows)
    summary = []
    for s in (0, 1):
        subset = [r for r in rows if r['sample'] == s]
        valid = [r for r in subset if r['valid']]
        b = sum(r['legacy_energy_mev'] for r in valid)
        n = sum(r['clean_energy_mev'] for r in valid)
        summary.append(dict(sample=s, n_events=len(subset), n_valid_events=len(valid),
                            n_invalid_events=len(subset)-len(valid), sum_legacy_mev=b, sum_clean_mev=n,
                            delta_mev=n-b, ratio=n/b if b else math.nan))
    write_tsv(out/'energy_summary.tsv', summary)
    by_channel = []
    for s in (0, 1):
        table = {r['cellid']: r for r in policies[s]}
        for c in np.unique(h['cellid'][h['sample'] == s]):
            m = (h['sample'] == s)&(h['cellid'] == c)
            usable = m&known
            p = table[int(c)]
            by_channel.append(dict(sample=s, cellid=int(c), pitch_um=pitch(int(c)),
                                   legacy_source=p['legacy_source'], clean_source=p['clean_source'],
                                   legacy_adc=p['legacy_adc'], clean_adc=p['clean_adc'],
                                   n_hits=int(m.sum()), n_invalid=int((m&~known).sum()),
                                   n_negative_legacy=int((usable&(baseline < 0)).sum()),
                                   valid_hit_subtotal_legacy=float(baseline[usable].sum()),
                                   valid_hit_subtotal_clean=float(corrected[usable].sum()),
                                   delta_valid_hit_subtotal=float((corrected[usable]-baseline[usable]).sum())))
    write_tsv(out/'channel_energy.tsv', by_channel)
    closure_rows = []
    for s in (0, 1):
        for c in np.unique(h['cellid'][h['sample'] == s]):
            m = (h['sample'] == s)&(h['cellid'] == c)
            v = relative[m&np.isfinite(relative)]
            closure_rows.append(dict(sample=s, cellid=int(c), n_hits=int(m.sum()),
                                     n_mismatch=int((m&~closure).sum()),
                                     median_relative_residual=float(np.median(v)) if len(v) else math.nan))
    write_tsv(out/'saved_calibration_check.tsv', closure_rows,
              ['sample', 'cellid', 'n_hits', 'n_mismatch', 'median_relative_residual'])
    return adc, temp_factor, old, new


def beam_report(out, b, e):
    rows = []
    for run in np.unique(b['run']):
        m = b['run'] == run
        row = dict(run=int(run), n_primaries=int(m.sum()),
                   n_non_axial=int(np.count_nonzero((np.abs(b['dx'][m])>1e-12)|(np.abs(b['dy'][m])>1e-12)|(np.abs(b['dz'][m]-1)>1e-12))))
        for key in ('px', 'py', 'pz', 'dx', 'dy', 'dz'):
            vals = b[key][m]
            row.update({key+'_mean': float(vals.mean()), key+'_std': float(vals.std()),
                        key+'_min': float(vals.min()), key+'_max': float(vals.max())})
        rows.append(row)
    write_tsv(out/'beam_truth.tsv', rows)
    reco = []
    for s in (0, 1):
        for run in np.unique(e['run'][e['sample'] == s]):
            m = (e['sample'] == s)&(e['run'] == run)
            row = dict(sample=s, run=int(run), n=int(m.sum()))
            for key in ('sx', 'sy', 'x0', 'y0', 'length', 'track_rms'):
                row[key+'_mean'] = float(e[key][m].mean())
                row[key+'_std'] = float(e[key][m].std())
            row['max_saved_length_difference_mm'] = float(np.max(np.abs(e['length'][m]-e['saved_length'][m])))
            reco.append(row)
    write_tsv(out/'beam_reconstructed.tsv', reco)


def analyze(ROOT, a, out, policies, pedestals, thresholds, digi):
    h = dict(ROOT.RDataFrame('Hits', str(out/'flat.root')).AsNumpy())
    e = dict(ROOT.RDataFrame('Events', str(out/'flat.root')).AsNumpy())
    b = dict(ROOT.RDataFrame('Beam', str(out/'flat.root')).AsNumpy())
    if not len(e['event_index']):
        raise ValueError('No selected muon events')
    beam_report(out, b, e)
    del b
    adc, tf, old, new = energy_audit(a, out, h, e, policies, pedestals)
    path_factor = 2/h['length']
    h.update(adc=adc*path_factor, adc20=adc*tf*path_factor,
             saved_mip=h['saved_hg']/.305*path_factor,
             truth=h['truth_edep']*path_factor, tf=tf, path_factor=path_factor)
    mask = (h['mip_selected'] == 1)&np.isfinite(h['adc'])&np.isfinite(h['adc20'])&(tf>0)
    selected = {k: v[mask] for k, v in h.items()}
    if not len(selected['cellid']):
        raise ValueError('No geometrically selected MIP hits')
    w, strata = common_weights(selected, a.long_bin, a.transverse_bin, a.slope_bin, a.min_bin)
    write_tsv(out/'matching_bins.tsv', strata)
    selected['weight'] = w
    ROOT.RDF.FromNumpy({k: np.ascontiguousarray(v) for k, v in selected.items()}).Snapshot('MatchedHits', str(out/'matched_hits.root'))
    coverage = []
    for s in (0, 1):
        m = selected['sample'] == s
        coverage.append(dict(sample=s, geometric_hits=int(m.sum()), overlap_hits=int(np.count_nonzero(w[m])),
                             overlap_sumw=float(w[m].sum()), overlap_fraction=float(np.count_nonzero(w[m])/m.sum()) if m.any() else 0))
    write_tsv(out/'matching_coverage.tsv', coverage)
    diagnostics = []
    for s in (0, 1):
        m = selected['sample'] == s
        for key in ('longitudinal', 'transverse', 'sx', 'sy', 'length', 'temperature'):
            values = selected[key][m]
            weights = w[m]
            mean = np.average(values, weights=weights) if weights.sum() else math.nan
            rms = np.sqrt(np.average((values-mean)**2, weights=weights)) if weights.sum() else math.nan
            diagnostics.append(dict(sample=s, variable=key, n=len(values),
                                    mean_before=float(np.mean(values)) if len(values) else math.nan,
                                    mean_matched=float(mean), std_matched=float(rms)))
    write_tsv(out/'matching_diagnostics.tsv', diagnostics)
    if a.no_fit:
        return []
    if not np.any(w > 0):
        raise ValueError('No common geometry bins; examine matching_bins.tsv or use a larger sample')
    load_fit_function(ROOT)
    rootout = ROOT.TFile.Open(str(out/'spectra.root'), 'CREATE')
    constants = [{r['cellid']: r for r in p} for p in policies]
    generated = {r['cellid']: r for r in digi}
    wanted = None if a.channels == 'all' else set(map(int, a.channels.split(',')))
    fits, comparison = [], []
    for c in np.unique(selected['cellid']):
        c = int(c)
        if wanted is not None and c not in wanted:
            continue
        masks = [(selected['cellid'] == c)&(selected['sample'] == s)&(w > 0) for s in (0, 1)]
        if not all(m.any() for m in masks):
            continue
        scale = 100. if pitch(c) == 10 else 400.
        maxima = []
        minima = {stage: [] for stage in ('adc', 'adc20')}
        for s, m in enumerate(masks):
            if c not in thresholds[s]:
                break
            th, sig, _ = thresholds[s][c]
            factor = selected['path_factor'][m]
            temp = selected['tf'][m]
            floor = np.maximum(th+a.plateau_sigmas*sig, 10/temp)
            minima['adc'].append(float(np.max(floor*factor)))
            minima['adc20'].append(float(np.max(floor*factor*temp)))
            maxima.append(float(np.min((2600-pedestals[s][c])*factor*np.minimum(temp, 1))))
        else:
            upper = min(8*scale, 2000. if pitch(c) == 15 else 800., *maxima)
            channel_fits = {}
            for s, m in enumerate(masks):
                for stage in ('adc', 'adc20', 'saved_mip', 'truth'):
                    if stage == 'truth' and s == 0:
                        continue
                    if stage in minima:
                        fit_scale, low, high = scale, max(minima[stage]), upper
                    elif stage == 'saved_mip':
                        denominator = constants[s][c]['legacy_adc']
                        if not math.isfinite(denominator) or denominator <= 0:
                            continue
                        fit_scale, low, high = scale/denominator, max(minima['adc20'])/denominator, upper/denominator
                    else:
                        fit_scale, low, high = .305, .08, .9
                    f = fit_spectrum(ROOT, rootout, f'{"data" if s==0 else "mc"}_{c}_{stage}',
                                     selected[stage][m], w[m], fit_scale, low, high, a.min_fit_hits)
                    row = dict(sample=s, cellid=c, pitch_um=pitch(c), stage=stage, **f)
                    fits.append(row); channel_fits[s, stage] = row
            row = dict(cellid=c, pitch_um=pitch(c), digi_source=generated[c]['digi_source'],
                       input_digi_adc=generated[c]['digi_adc'])
            for s, label in ((0, 'data'), (1, 'mc')):
                p = constants[s][c]
                for key in ('legacy_source', 'clean_source', 'legacy_adc', 'clean_adc'):
                    row[label+'_'+key] = p[key]
                f = channel_fits.get((s, 'adc20'), {})
                row[label+'_fit_valid'] = f.get('valid', False)
                for key in ('mpv', 'mpv_err'):
                    v = f.get(key, math.nan)
                    row[label+'_adc20_'+key] = v
                    for policy in ('legacy', 'clean'):
                        denom = p[policy+'_adc']
                        row[label+'_'+policy+'_mip_'+key] = v/denom if math.isfinite(denom) and denom > 0 else math.nan
                row[label+'_saved_mip_mpv'] = channel_fits.get((s, 'saved_mip'), {}).get('mpv', math.nan)
            fm = channel_fits.get((1, 'adc'), {})
            ft = channel_fits.get((1, 'truth'), {})
            denom = generated[c]['digi_adc']
            row['mc_adc_fit_valid'] = fm.get('valid', False)
            row['mc_truth_fit_valid'] = ft.get('valid', False)
            row['mc_adc_mpv_over_input'] = fm.get('mpv', math.nan)/denom if math.isfinite(denom) and denom > 0 else math.nan
            row['mc_adc_ratio_stat_error'] = fm.get('mpv_err', math.nan)/denom if math.isfinite(denom) and denom > 0 else math.nan
            row['mc_truth_mpv_mev'] = ft.get('mpv', math.nan)
            row['mc_truth_mpv_over_0p305'] = ft.get('mpv', math.nan)/.305
            row['adc_ratio_div_truth_ratio'] = row['mc_adc_mpv_over_input']/row['mc_truth_mpv_over_0p305'] if row['mc_truth_mpv_over_0p305'] > 0 else math.nan
            r = row['mc_adc_mpv_over_input']
            row['primary_pe_scale_trial_only'] = 1/r if fm.get('valid') and r > 0 else math.nan
            for policy in ('legacy', 'clean'):
                d, m = row[f'data_{policy}_mip_mpv'], row[f'mc_{policy}_mip_mpv']
                row[policy+'_data_over_mc'] = d/m if m > 0 else math.nan
            comparison.append(row)
            if len(comparison) % 25 == 0:
                print(f'Fitted {len(comparison)} common channels', flush=True)
            continue
        # Missing efficiency constants make a plateau claim impossible; keep an explicit record.
        fits.append(dict(sample=-1, cellid=c, pitch_um=pitch(c), stage='all',
                         **dict(n=0, sumw=0, effective_n=0, fit_low=math.nan, fit_high=math.nan,
                                mpv=math.nan, mpv_err=math.nan, convolution_mode=math.nan, width=math.nan,
                                sigma=math.nan, chi2=math.nan, ndf=0, status=-1, covariance=-1,
                                valid=False, reason='missing_threshold', underflow=0, overflow=0, fit_effective_n=0)))
    rootout.Close()
    if fits:
        write_tsv(out/'fits.tsv', fits)
    if comparison:
        write_tsv(out/'channel_closure.tsv', comparison)
    return comparison


def main():
    a = arguments()
    import ROOT
    ROOT.gROOT.SetBatch(True)
    manifest, metadata, excluded = discover(ROOT, a)
    print(f'Usable run pairs: {len(manifest)//2}; excluded: {len(excluded)}', flush=True)
    for row in excluded:
        print('Excluded:', row, flush=True)
    if a.dry_run:
        print(json.dumps(metadata, indent=2)); return
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    snapshot(a, out)
    write_tsv(out/'inputs/files.tsv', metadata)
    write_tsv(out/'inputs/excluded_runs.tsv', excluded, ['run', 'reason'])
    if not manifest:
        raise ValueError('No usable run pairs. See inputs/excluded_runs.tsv. Recovered inputs require --allow-recovered for a diagnostic run.')
    manifest_path = out/'inputs/manifest.tsv'
    manifest_path.write_text(''.join('\t'.join(row)+'\n' for row in manifest))
    for run in sorted({int(row[1]) for row in manifest}):
        p = MC/f'CEPCScECAL/generate/mac/mu-/100GeV/ECAL_Run{run}.mac'
        if p.is_file():
            shutil.copyfile(p, out/f'inputs/ECAL_Run{run}.mac')
    policies, pedestals, thresholds = [], [], []
    for s, label in enumerate(('data', 'mc')):
        rows, means = read_mip(ROOT, getattr(a, label+'_mip'))
        write_tsv(out/f'{label}_constants.tsv', rows)
        write_tsv(out/f'{label}_fallback_means.tsv', means)
        write_tsv(out/f'{label}_changed_channels.tsv', [r for r in rows if r['changed']], list(rows[0]))
        policies.append(rows)
        pedestals.append(read_pedestal(ROOT, getattr(a, label+'_pedestal')))
        thresholds.append(read_threshold(ROOT, getattr(a, label+'_threshold')))
    digi, means = read_mip(ROOT, a.digi_mip)
    write_tsv(out/'digitization_constants.tsv', digi)
    subprocess.run(['make', '-C', str(HERE)], check=True)
    command = [str(HERE/'CollectMuon'), str(manifest_path), str(out/'flat.root'), str(out/'validation.tsv'),
               str(a.max_tracks), str(a.max_angle), str(a.max_track_rms), str(int(a.allow_recovered))]
    (out/'inputs/collector_command.json').write_text(json.dumps(command, indent=2)+'\n')
    subprocess.run(command, check=True)
    comparison = analyze(ROOT, a, out, policies, pedestals, thresholds, digi)
    from report import report
    report(ROOT, out, comparison)
    included = {int(row[1]) for row in manifest}
    diagnostic = any(r['recovered'] and r['run'] in included for r in metadata)
    (out/'STATUS.json').write_text(json.dumps(dict(completed=True, diagnostic_recovered_inputs=diagnostic,
                                                 excluded_runs=excluded, sampled_tracks_per_file=a.max_tracks,
                                                 fitted_channels=len(comparison),
                                                 valid_mc_adc_channels=sum(bool(r['mc_adc_fit_valid']) for r in comparison),
                                                 valid_paired_channels=sum(bool(r['data_fit_valid'] and r['mc_fit_valid']) for r in comparison),
                                                 no_fit=a.no_fit), indent=2)+'\n')
    print(f'Completed: {out}', flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr)
        raise
