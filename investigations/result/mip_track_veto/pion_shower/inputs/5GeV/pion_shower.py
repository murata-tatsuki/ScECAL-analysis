#!/usr/bin/env python3
"""Investigate a MIP-like entrance track followed by shower development."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import warnings

import awkward as ak
import numpy as np
import uproot

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('track_investigation', HERE/'check.py')
track = importlib.util.module_from_spec(spec)
spec.loader.exec_module(track)
helper = track.helper
PROJECT, PREVIOUS = track.PROJECT, track.PREVIOUS
RESULT = track.RESULT/'pion_shower'
MIP = .305
CONFIGS = dict(nominal=8, loose=6, tight=10, interaction_only=8, mip_only=None)
COGS = track.COGS
require, metadata, write_csv = track.require, track.metadata, track.write_csv


def sample_location(out, sample):
    return out/'pion_reference/data' if sample == 'pion_run' else out/sample


def prefix_rms(energy, q, width, end, view):
    columns = np.arange(view, end, 2)
    use = energy[:, columns] > 0
    n = use.sum(axis=1)
    z = track.LAYER_Z[columns][None, :]
    q = np.nan_to_num(q[:, columns])
    w = np.nan_to_num(width[:, columns])
    denominator = np.maximum(n, 1)
    mz, mq = (use*z).sum(axis=1)/denominator, (use*q).sum(axis=1)/denominator
    dz = z-mz[:, None]
    zz = (use*dz*dz).sum(axis=1)
    slope = np.divide((use*dz*(q-mq[:, None])).sum(axis=1), zz,
                      out=np.zeros(len(n)), where=zz > 0)
    residual = q-mq[:, None]-slope[:, None]*dz
    rms = np.sqrt((use*(residual**2+w*w)).sum(axis=1)/denominator)
    rms[(n < 3)|(zz <= 0)] = np.inf
    return rms


def classify(a):
    energy = np.nan_to_num(a['layer_energy_MeV'])/MIP
    hits = np.nan_to_num(a['layer_strips'])
    n = len(energy)
    starts = np.zeros((n, 26), dtype=bool)
    for end in range(6, 26):
        prefix = energy[:, :end]
        active = prefix > 0
        count = active.sum(axis=1)
        enough = (count >= int(np.ceil(.75*end))) & (a['evaluated'])
        candidates = np.flatnonzero(enough)
        if not len(candidates):
            continue
        median = np.zeros(n)
        median[candidates] = np.nanmedian(np.where(active[candidates], prefix[candidates], np.nan), axis=1)
        average = prefix.sum(axis=1)/np.maximum(count, 1)
        occupancy = hits[:, :end].sum(axis=1)/np.maximum(count, 1)
        rx = prefix_rms(energy, a['cog_x_mm'], a['width_x_mm'], end, 1)
        ry = prefix_rms(energy, a['cog_y_mm'], a['width_y_mm'], end, 0)
        before = enough&(median >= .3)&(median <= 2.5)&((prefix > 4).sum(axis=1) <= 1)&\
                 (occupancy <= 2.5)&(rx <= 5.3)&(ry <= 5.3)
        onset = (energy[:, end] >= 6)&(energy[:, end+1] >= 6)&\
                (hits[:, end] >= 4)&(hits[:, end+1] >= 4)
        ratio = energy[:, end:end+4].mean(axis=1)/np.maximum(average, 1e-12)
        starts[:, end] = before&onset&(ratio >= 5)
    for name, minimum in CONFIGS.items():
        if name == 'mip_only':
            continue
        interaction = starts[:, minimum:26].any(axis=1)
        a[name] = interaction if name == 'interaction_only' else interaction|a['mip_only']
    a['shower_start'] = np.where(starts[:, 8:26].any(axis=1), np.argmax(starts[:, 8:26], axis=1)+8, -1)
    a['loose_shower_start'] = np.where(starts[:, 6:26].any(axis=1), np.argmax(starts[:, 6:26], axis=1)+6, -1)
    return a


def input_plan(beam, sample):
    if sample != 'pion_run':
        prior = PREVIOUS/'inputs'/f'{beam}GeV'
        require(json.loads((prior/'complete.json').read_text()).get('current_mc_tail_audited'), 'Prior MC audit missing')
        with np.load(PREVIOUS/sample/'diagnostics'/f'{beam}GeV/events.npz') as saved:
            a = {k: saved[k].copy() for k in ('file_index', 'entry', 'run', 'event_time', 'event',
                                            'tail_rich', 'energy_GeV', 'cog_extent_mm')}
        records = json.loads((prior/sample/'manifest.json').read_text())
        eligible = (~a['tail_rich'].astype(bool))&(a['energy_GeV'] > 0)&(a['cog_extent_mm'] < 200)
    else:
        source = PROJECT.parent/'data_beamtest/ECAL_data/analysed/2023'/('ps' if beam == 5 else 'sps')/'calib_original/pi-'/f'{beam}GeV'
        files = sorted(source.glob('*.root'))
        require(files, 'No pion-run reference inputs')
        counts = []
        for path in files:
            with uproot.open(path) as root:
                counts.append(root['Calib_Hit'].num_entries)
        target = min(6000, sum(counts))
        desired = target*np.array(counts)/sum(counts)
        take = np.floor(desired).astype(int)
        for index in np.argsort(-(desired-take))[:target-int(take.sum())]:
            take[index] += 1
        rng = np.random.default_rng(20261111+beam*10+2)
        records, file_indices, entries = [], [], []
        for index, (path, count, count_take) in enumerate(zip(files, counts, take)):
            spans = helper.sample_blocks(count, int(count_take), rng)
            records.append(dict(file_index=index, entries=count, sampled=int(count_take), spans=spans, **metadata(path)))
            for lo, hi in spans:
                file_indices.extend([index]*(hi-lo))
                entries.extend(range(lo, hi))
        a = dict(file_index=np.array(file_indices), entry=np.array(entries))
        for key in ('run', 'event_time', 'event', 'energy_GeV', 'cog_extent_mm'):
            a[key] = np.full(target, np.nan)
        a['tail_rich'] = np.zeros(target, dtype=bool)
        eligible = np.ones(target, dtype=bool)
    return a, records, eligible


def collect(beam, sample, out):
    a, records, eligible = input_plan(beam, sample)
    mask_source = PREVIOUS/'inputs'/f'{beam}GeV/resolution_source.cc'
    text = mask_source.read_text().split('const unordered_set<int> bad_channels = {', 1)[1].split('};', 1)[0]
    bad = np.array(sorted(set(map(int, re.findall(r'\d+', re.sub(r'//[^\n]*', '', text))))))
    n = len(a['entry'])
    for key in ('layer_energy_MeV', 'layer_strips', 'cog_x_mm', 'cog_y_mm', 'width_x_mm', 'width_y_mm'):
        a[key] = np.full((n, 30), np.nan)
    a['mip_only'] = np.zeros(n, dtype=bool)
    a['evaluated'] = eligible.copy()
    read_count, max_difference = 0, 0.
    for index, record in enumerate(records):
        path = Path(record['path'])
        original = {key: record[key] for key in ('path', 'bytes', 'mtime_ns')}
        require(metadata(path) == original, 'Input changed: '+str(path))
        selected = np.flatnonzero(eligible&(a['file_index'] == index))
        if not len(selected):
            continue
        with uproot.open(path, array_cache=None) as root:
            tree = root['Calib_Hit']
            require(tree.num_entries == record['entries'], 'Entry count changed')
            for lo, hi in record['spans']:
                rows = selected[(a['entry'][selected] >= lo)&(a['entry'][selected] < hi)]
                if not len(rows):
                    continue
                start, stop = int(a['entry'][rows].min()), int(a['entry'][rows].max())+1
                arrays = tree.arrays(['Run_Num', 'Event_Time', 'Event_Num', 'CellID', 'Hit_Energy',
                                     'Hit_X', 'Hit_Y', 'Hit_Z'], entry_start=start, entry_stop=stop, library='ak')
                for row in rows:
                    j = int(a['entry'][row])-start
                    keys = tuple(int(arrays[k][j]) for k in ('Run_Num', 'Event_Time', 'Event_Num'))
                    if sample == 'pion_run':
                        a['run'][row], a['event_time'][row], a['event'][row] = keys
                    else:
                        require(keys == tuple(int(a[k][row]) for k in ('run', 'event_time', 'event')), 'Event identity differs')
                    ids, en, x, y, z = [ak.to_numpy(arrays[k][j]) for k in ('CellID', 'Hit_Energy', 'Hit_X', 'Hit_Y', 'Hit_Z')]
                    require(len({len(v) for v in (ids, en, x, y, z)}) == 1, 'Hit vector length mismatch')
                    values = track.moments(ids, en, x, y, z, bad)
                    recomputed = values['recomputed_energy_GeV']
                    if sample == 'pion_run':
                        a['energy_GeV'][row] = recomputed
                        coords = np.array([values['cog_x_mm'][9], values['cog_y_mm'][9],
                                           values['cog_x_mm'][10], values['cog_y_mm'][10]])
                        a['cog_extent_mm'][row] = np.max(np.abs(coords)) if np.isfinite(coords).all() else np.nan
                    else:
                        difference = abs(recomputed-a['energy_GeV'][row])
                        max_difference = max(max_difference, difference)
                        require(difference < 1e-10, 'Signal differs from previous baseline')
                    for key in ('layer_energy_MeV', 'cog_x_mm', 'cog_y_mm', 'width_x_mm', 'width_y_mm'):
                        a[key][row] = values[key]
                    layer, chip, channel = ids//100000, ids%100000//10000, ids%100
                    good = (ids >= 0)&(layer < 30)&(chip < 6)&(channel < 36)&~((chip == 5)&(channel > 29))&\
                           (en > 0)&~np.isin(ids, bad)
                    cells = np.unique(helper.physical(ids[good]))
                    a['layer_strips'][row] = np.bincount(cells//100000, minlength=30)
                    a['mip_only'][row] = 0 < recomputed <= .03 and values['hits'] <= 60 and values['active'] >= 15 and\
                        values['front'] >= 2 and values['back'] >= 2 and values['occupancy'] <= 2 and\
                        values['rms_x_mm'] <= 5.3 and values['rms_y_mm'] <= 5.3
                    read_count += 1
        require(metadata(path) == original, 'Input changed during reading')
        if sample != 'mc' or (index+1) % 20 == 0:
            print(f'{beam}GeV {sample}: {index+1}/{len(records)} files; {read_count}/{eligible.sum()} evaluated', flush=True)
    require(read_count == int(eligible.sum()), 'Incomplete moment coverage')
    if sample != 'pion_run':
        with np.load(track.RESULT/sample/'diagnostics'/f'{beam}GeV/events.npz') as previous_track:
            require(np.array_equal(a['mip_only'], previous_track['nominal']), 'Previous MIP veto not reproduced')
    a = classify(a)
    directory = sample_location(out, sample)/'diagnostics'/f'{beam}GeV'
    directory.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(directory/'events.npz', **a)
    record_out = out/'inputs'/f'{beam}GeV'/sample
    record_out.mkdir(parents=True, exist_ok=True)
    (record_out/'manifest.json').write_text(json.dumps(records, indent=2)+'\n')
    (record_out/'validation.json').write_text(json.dumps(dict(sampled=n, evaluated=read_count,
        max_abs_signal_difference_GeV=max_difference, mask_sha256=hashlib.sha256(mask_source.read_bytes()).hexdigest(),
        previous_baseline=str(PREVIOUS), input_identity='All input size/mtime and evaluated Run/Time/Event checked',
        pion_reference_caveat='Pion-setting runs: PID, matched tail flags and calibration-production equivalence unavailable; shape reference only' if sample == 'pion_run' else None), indent=2)+'\n')
    return a


def reference_rates(beam, a, out):
    rows = []
    for cog in COGS:
        use = (a['energy_GeV'] > 0)&(a['cog_extent_mm'] < cog)
        for name in CONFIGS:
            n, tagged = int(use.sum()), int(a[name][use].sum())
            rows.append(dict(energy_GeV=beam, sample='pion_setting_reference', cog_mm=cog,
                condition=name, before=n, removed=tagged, removed_fraction=tagged/n if n else np.nan,
                caveat='Not pure-pion efficiency; older calibration; no matched HG-LG tail selection'))
    write_csv(out/'comparison'/f'{beam}GeV/pion_reference_rates.csv', rows)


def plots(beam, sample, a, out):
    import matplotlib.pyplot as plt
    location = sample_location(out, sample)/'diagnostics'/f'{beam}GeV'
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for axis, cog in zip(axes.flat, COGS):
        use = (~a['tail_rich'].astype(bool))&np.isfinite(a['energy_GeV'])&(a['cog_extent_mm'] < cog)
        values = a['energy_GeV'][use]
        if not len(values):
            continue
        edges = np.linspace(0, max(np.quantile(values, .999)*1.1, .001), 151)
        for v, label, color, style in ((values, 'Before', 'black', '-'),
            (values[~a['nominal'][use]], 'After combined veto', 'tab:blue', '--'),
            (values[a['interaction_only'][use]], 'Track then shower', 'tab:orange', '-'),
            (values[a['mip_only'][use]], 'Through-going MIP', 'tab:red', ':')):
            axis.hist(v, bins=edges, histtype='step', color=color, linestyle=style, label=f'{label} (n={len(v)})')
        axis.set(yscale='log', ylim=(.7, None), xlabel='Calibrated signal sum [GeV]', ylabel='Events / bin',
                 title=('nocut' if cog == 200 else f'CoG {cog} mm'))
        axis.legend(fontsize=7)
    fig.suptitle(f'{beam} GeV {sample}: MIP track followed by shower; fixed cuts')
    fig.tight_layout()
    fig.savefig(location/'energy_before_after.png', dpi=140)
    plt.close(fig)
    candidates = np.flatnonzero(a['interaction_only']&(a['cog_extent_mm'] < 20))[:3]
    for i, row in enumerate(candidates):
        start = int(a['shower_start'][row])
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        for axis, view, coordinate in ((axes[0], 1, 'x'), (axes[1], 0, 'y')):
            layers = np.arange(view, 30, 2)
            layers = layers[a['layer_energy_MeV'][row, layers] > 0]
            axis.errorbar(track.LAYER_Z[layers], a['cog_'+coordinate+'_mm'][row, layers],
                          yerr=a['width_'+coordinate+'_mm'][row, layers], fmt='o', markersize=3)
            axis.axvline(track.LAYER_Z[start], color='tab:red', linestyle='--')
            axis.set(xlabel='Layer z [mm]', ylabel=coordinate+' [mm]')
        axes[2].plot(range(30), a['layer_energy_MeV'][row]/MIP, 'o-', label='Signal [MIP]')
        axes[2].axvline(start, color='tab:red', linestyle='--')
        axes[2].axhline(6, color='0.5', linestyle=':')
        axes[2].set(yscale='log', ylim=(.1, None), xlabel='Layer', ylabel='Signal [MIP]')
        other = axes[2].twinx()
        other.plot(range(30), a['layer_strips'][row], color='tab:green', linestyle=':', label='Unique strips')
        other.set_ylabel('Positive unique strips', color='tab:green')
        fig.suptitle(f'{beam} GeV {sample}; run {a["run"][row]:.0f}, entry {a["entry"][row]}; '
                     f'start layer {start}, signal {a["energy_GeV"][row]:.4g} GeV')
        fig.tight_layout()
        fig.savefig(location/f'candidate_{i}.png', dpi=140)
        plt.close(fig)


def combine(out):
    import matplotlib.pyplot as plt
    tables = {}
    for name in ('metrics', 'tag_rates', 'paired_changes', 'pion_reference_rates'):
        rows = []
        for beam in (5, 40, 100):
            require(json.loads((out/'inputs'/f'{beam}GeV/complete.json').read_text())['complete'], 'Incomplete energy')
            with (out/'comparison'/f'{beam}GeV'/(name+'.csv')).open() as stream:
                rows.extend(csv.DictReader(stream))
        tables[name] = rows
        write_csv(out/'comparison'/(name+'.csv'), rows)
    summary = []
    for beam in (5, 40, 100):
        ms = {r['stage']: r for r in tables['metrics'] if int(r['energy_GeV']) == beam and r['sample'] == 'data' and r['cog_mm'] == '20' and r['condition'] == 'nominal'}
        rates = {(r['sample'], r['condition']): r for r in tables['tag_rates'] if int(r['energy_GeV']) == beam and r['cog_mm'] == '20'}
        row = dict(energy_GeV=beam, data_before=rates['data','nominal']['before'],
            data_new_shower_removed=rates['data','interaction_only']['removed'], data_mip_removed=rates['data','mip_only']['removed'],
            data_combined_removed=rates['data','nominal']['removed'], data_combined_removed_percent=100*float(rates['data','nominal']['removed_fraction']),
            electron_mc_before=rates['mc','nominal']['before'], electron_mc_removed=rates['mc','nominal']['removed'],
            electron_mc_removed_percent=100*float(rates['mc','nominal']['removed_fraction']))
        for key in ('rms_mean', 'w68_median', 'gaussian_core'):
            for stage in ('before', 'after'):
                row[key+'_'+stage+'_percent'] = 100*float(ms[stage][key])
        summary.append(row)
    write_csv(out/'comparison/summary20mm.csv', summary)
    for measure, label in (('gaussian_core', 'Gaussian core sigma / mean [%]'), ('rms_mean', 'Full RMS / mean [%]')):
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        for axis, beam in zip(axes, (5, 40, 100)):
            for sample, color in (('data', 'black'), ('mc', 'tab:red')):
                for stage, style, marker in (('before', '-', 'o'), ('after', '--', 's')):
                    rows = [r for r in tables['metrics'] if int(r['energy_GeV']) == beam and r['sample'] == sample and r['condition'] == 'nominal' and r['stage'] == stage]
                    rows.sort(key=lambda r: COGS.index(int(r['cog_mm'])))
                    axis.plot(range(4), [100*float(r[measure]) for r in rows], color=color,
                              linestyle=style, marker=marker, label=f'{sample} {stage}')
            axis.set(xticks=range(4), xticklabels=['nocut', '20', '10', '5'], xlabel='CoG half-width [mm]', ylabel=label, title=f'{beam} GeV')
            axis.legend(fontsize=8)
            axis.grid(alpha=.2)
        fig.tight_layout()
        fig.savefig(out/'comparison'/(measure+'_summary.png'), dpi=160)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energy', type=int, choices=(5, 40, 100))
    parser.add_argument('--combine', action='store_true')
    parser.add_argument('--summarize-only', action='store_true')
    parser.add_argument('--out', type=Path, default=RESULT)
    parser.add_argument('--bootstrap', type=int, default=400)
    args = parser.parse_args()
    require(args.combine or args.energy is not None, 'Specify energy or combine')
    out = args.out.resolve()
    record = out/'inputs'/('combined' if args.combine else f'{args.energy}GeV')
    record.mkdir(parents=True, exist_ok=True)
    require(not (record/'complete.json').exists() or args.summarize_only, 'Complete output exists')
    for source in (Path(__file__), HERE/'check.py', track.HELPER_PATH):
        destination = record/('width_helper.py' if source == track.HELPER_PATH else 'track_helper.py' if source == HERE/'check.py' else source.name)
        if not args.summarize_only:
            shutil.copy2(source, destination)
    (record/('summary_configuration.json' if args.summarize_only else 'configuration.json')).write_text(json.dumps(dict(
        command=sys.argv, configs=CONFIGS, prefix_active_fraction=.75, prefix_median_MIP=[.3,2.5], prefix_max4MIP_outliers=1,
        prefix_strips_per_active_layer=2.5, prefix_rms_mm=5.3, prefix_min_layers_per_projection=3,
        onset_layers=2, onset_min_MIP=6, onset_min_unique_positive_strips=4, next4_average_ratio_minimum=5,
        ratio_denominator='Prefix active-layer average', total_energy_window=None,
        note='Nominal is previous through-going MIP veto OR MIP-prefix then shower; not pure-pion PID',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')
    with tempfile.TemporaryDirectory(prefix='scecal-pion-veto-') as cache:
        fonts = Path(importlib.util.find_spec('matplotlib').origin).parent/'mpl-data/fonts/ttf'
        fc = Path(cache)/'fonts.conf'
        fc.write_text(f'<fontconfig><dir>{fonts}</dir><cachedir>{cache}</cachedir></fontconfig>')
        os.environ['FONTCONFIG_FILE'], os.environ['MPLCONFIGDIR'] = str(fc), cache
        import matplotlib
        matplotlib.use('Agg')
        if args.combine:
            combine(out)
        else:
            tables = [[], [], []]
            track.CONFIGS = CONFIGS
            for sample in ('data', 'mc', 'pion_run'):
                if args.summarize_only:
                    with np.load(sample_location(out, sample)/'diagnostics'/f'{args.energy}GeV/events.npz') as saved:
                        a = classify({key: saved[key] for key in saved.files})
                else:
                    a = collect(args.energy, sample, out)
                if sample == 'pion_run':
                    reference_rates(args.energy, a, out)
                else:
                    for rows, more in zip(tables, track.summarize(args.energy, sample, a, out, args.bootstrap)):
                        rows.extend(more)
                plots(args.energy, sample, a, out)
            for name, rows in zip(('metrics', 'tag_rates', 'paired_changes'), tables):
                write_csv(out/'comparison'/f'{args.energy}GeV'/(name+'.csv'), rows)
            for row in tables[1]:
                if row['cog_mm'] == 20:
                    print(json.dumps(row), flush=True)
    (record/'complete.json').write_text(json.dumps(dict(complete=True))+'\n')


if __name__ == '__main__':
    main()
