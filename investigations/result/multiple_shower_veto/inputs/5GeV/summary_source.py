#!/usr/bin/env python3
"""Small, reproducible investigation of persistent separated entrance components."""
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
import time
import warnings

import awkward as ak
import numpy as np
import uproot
from scipy.optimize import curve_fit, OptimizeWarning

PROJECT = Path(__file__).resolve().parents[2]
DATA = PROJECT.parent / 'data_beamtest'
RESULT = PROJECT / 'investigations/result/multiple_shower_veto'
SOURCE = PROJECT / 'analysis/resolution/SingleEnergyAnalysis.cc'
PITCH = 5.3
MIP = .305
COGS = (200, 20, 10, 5)
# Fixed before examining the results; no optimization on resolution.
CONFIGS = {
    'nominal': (3, .15, 3),
    'two_layers': (3, .15, 2),
    'wide_separation': (4, .15, 3),
    'fraction10': (3, .10, 3),
    'fraction20': (3, .20, 3),
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metadata(path):
    stat = path.stat()
    return dict(path=str(path.resolve()), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def physical(ids):
    return ids // 100000 * 100000 + ids % 100000 // 10000 * 10000 + ids % 100


def relocated(text):
    old = '/net/megsrv01/export/megraid01/users/data_beamtest/'
    if text.startswith(old):
        text = str(DATA) + '/' + text[len(old):]
    path = Path(text).resolve(strict=True)
    require(path.is_relative_to(DATA), 'Unexpected input location: ' + str(path))
    return path


def projections(ids, energy, x, y, good):
    """Use short coordinates only; combine memory cells and long-axis segments."""
    result = []
    layers = ids // 100000
    for layer in range(8):
        use = (layers == layer) & good & (energy > 0)
        q = (x if layer % 2 else y)[use]
        e = energy[use]
        if not len(e):
            result.append((np.empty(0), np.empty(0)))
            continue
        q, reverse = np.unique(np.round(q, 6), return_inverse=True)
        e = np.bincount(reverse, weights=e)
        result.append((q, e))
    return result


def layer_pairs(projection, separation, fraction):
    q, energy = projection
    total = energy.sum()
    keep = energy >= .5 * MIP
    q, energy = q[keep], energy[keep]
    if len(q) < 2 or total <= 0:
        return []
    groups = np.split(np.arange(len(q)), np.flatnonzero(np.diff(q) > PITCH + .01) + 1)
    clusters = []
    for group in groups:
        value = energy[group].sum()
        if value >= MIP and value >= fraction * total:
            clusters.append((float(np.average(q[group], weights=energy[group])), value))
    return [(a[0], b[0]) for i, a in enumerate(clusters) for b in clusters[i+1:]
            if b[0] - a[0] >= separation * PITCH - .01]


def persistent(profiles, config, view):
    separation, fraction, needed = config
    pairs = [layer_pairs(profiles[layer], separation, fraction) for layer in range(view, 8, 2)]
    found = []
    for choices in pairs:
        for a, b in choices:
            support = sum(any(abs(c-a) <= PITCH+.01 and abs(d-b) <= PITCH+.01
                              for c, d in other) for other in pairs)
            if support >= needed:
                found.append((a, b, support))
    return found


def tag(profiles, unmasked_profiles, config):
    """Require the same two axes with and without the bad-channel mask."""
    flags, axes = [], []
    for view in (0, 1):
        masked = persistent(profiles, config, view)
        unmasked = persistent(unmasked_profiles, config, view)
        common = [a for a in masked if any(abs(a[0]-b[0]) <= PITCH+.01 and
                  abs(a[1]-b[1]) <= PITCH+.01 for b in unmasked)]
        flags.append(bool(common))
        axes.append(max(common, key=lambda a: (a[2], a[1]-a[0])) if common else (np.nan, np.nan, 0))
    return flags, axes


def sample_blocks(count, take, rng):
    if take == count:
        return [(start, min(start+1000, count)) for start in range(0, count, 1000)]
    blocks = min(12, max(1, take // 40))
    output = []
    for i in range(blocks):
        lo, hi = count*i//blocks, count*(i+1)//blocks
        size = take//blocks + (i < take % blocks)
        require(hi-lo >= size, 'Invalid sample allocation')
        start = int(rng.integers(lo, hi-size+1))
        output.append((start, start+size))
    return output


def collect(beam, sample, target, out):
    source = PROJECT / f'analysis/result/channel_response/hg_lg_tail/threshold/e_ssa/{beam}GeV'
    text = SOURCE.read_text().split('const unordered_set<int> bad_channels = {', 1)[1].split('};', 1)[0]
    text = re.sub(r'//[^\n]*', '', text)
    bad = np.array(sorted(set(map(int, re.findall(r'\d+', text)))), dtype=np.int64)
    manifest = source / 'inputs' / f'{sample}.tsv'
    files = [relocated(line.split('\t')[1]) for line in manifest.read_text().splitlines() if line]
    with uproot.open(source / 'tail_study.root') as root:
        require(root['min_tail_channels'].member('fVal') == 3 and root['nsigma'].member('fVal') == 5,
                'Unexpected HG-LG tail definition')
        tail = root['events'].arrays(['sample', 'file_index', 'entry', 'Run_Num', 'Event_Time',
                                     'TriggerID', 'tail_rich', 'low_tail_channels'], library='np')
    sample_id = 0 if sample == 'data' else 1
    tail = {key: value[tail['sample'] == sample_id] for key, value in tail.items()}
    counts, lookup = [], []
    for index, path in enumerate(files):
        selection = np.flatnonzero(tail['file_index'] == index)
        selection = selection[np.argsort(tail['entry'][selection])]
        require(np.array_equal(tail['entry'][selection], np.arange(len(selection))), 'Incomplete tail coverage')
        require(np.array_equal(tail['tail_rich'][selection], tail['low_tail_channels'][selection] >= 3),
                'Tail flags disagree with channel multiplicity')
        with uproot.open(path) as root:
            count = root['Calib_Hit'].num_entries
        require(count == len(selection), 'Input and tail entry counts differ: ' + str(path))
        counts.append(count)
        lookup.append(selection)
    total = sum(counts)
    target = min(total, target) if target else total
    desired = target * np.array(counts) / total
    takes = np.floor(desired).astype(int)
    for index in np.argsort(-(desired-takes))[:target-int(takes.sum())]:
        takes[index] += 1
    rng = np.random.default_rng(20261111 + beam*10 + sample_id)
    rows, records, displays = [], [], []
    input_dir = out / 'inputs' / f'{beam}GeV' / sample
    input_dir.mkdir(parents=True, exist_ok=True)
    branches = ['Run_Num', 'Event_Time', 'Event_Num', 'DetectorID', 'CellID', 'Hit_Energy', 'Hit_X', 'Hit_Y']
    for index, (path, count, take) in enumerate(zip(files, counts, takes)):
        spans = sample_blocks(count, int(take), rng)
        before = metadata(path)
        records.append(dict(file_index=index, entries=count, sampled=int(take), spans=spans, **before))
        with uproot.open(path, array_cache=None) as root:
            tree = root['Calib_Hit']
            for start, stop in spans:
                arrays = tree.arrays(branches, entry_start=start, entry_stop=stop, library='ak')
                for j, entry in enumerate(range(start, stop)):
                    k = lookup[index][entry]
                    scalar = lambda field: int(arrays[field][j])
                    require((scalar('Run_Num'), scalar('Event_Time'), scalar('Event_Num')) ==
                            (int(tail['Run_Num'][k]), int(tail['Event_Time'][k]), int(tail['TriggerID'][k])),
                            f'Event identity mismatch: {path}:{entry}')
                    row = dict(file_index=index, entry=entry, run=scalar('Run_Num'),
                               event_time=scalar('Event_Time'), event=scalar('Event_Num'),
                               tail_rich=bool(tail['tail_rich'][k]))
                    if row['tail_rich']:
                        rows.append(row)
                        continue
                    ids, en, x, y = (ak.to_numpy(arrays[field][j]) for field in
                                     ['CellID', 'Hit_Energy', 'Hit_X', 'Hit_Y'])
                    require(len(ids) == len(en) == len(x) == len(y), 'Mismatched hit vector lengths')
                    require(np.isfinite(en).all() and np.isfinite(x).all() and np.isfinite(y).all(),
                            'Non-finite calibrated hit')
                    layer, chip, channel = ids//100000, ids%100000//10000, ids%100
                    connected = (ids >= 0) & (layer < 30) & (chip < 6) & (channel < 36) & ~((chip == 5) & (channel > 29))
                    # Keep the exact current resolution mask lookup, including original CellID.
                    good = connected & (en >= 0) & ~np.isin(ids, bad)
                    row['energy_GeV'] = float(en[good].sum()/1000)
                    cog = []
                    for ilayer in (9, 10):
                        use = good & (layer == ilayer)
                        value = en[use].sum()
                        cog.extend([float(np.dot(x[use], en[use])/value),
                                    float(np.dot(y[use], en[use])/value)] if value > 0 else [np.nan, np.nan])
                    row['cog_extent_mm'] = float(max(map(abs, cog))) if np.isfinite(cog).all() else np.nan
                    profiles = projections(ids, en, x, y, good)
                    raw_profiles = projections(ids, en, x, y, connected & (en >= 0))
                    row['front_active_x'] = sum(profiles[l][1].sum() >= MIP for l in (1, 3, 5, 7))
                    row['front_active_y'] = sum(profiles[l][1].sum() >= MIP for l in (0, 2, 4, 6))
                    nominal_axes = None
                    for name, config in CONFIGS.items():
                        flags, axes = tag(profiles, raw_profiles, config)
                        row[name] = any(flags)
                        if name == 'nominal':
                            nominal_axes = axes
                            row['nominal_y'], row['nominal_x'] = flags
                            for view, axis in zip(('y', 'x'), axes):
                                row[view+'_axis1_mm'], row[view+'_axis2_mm'], row[view+'_support'] = axis
                    if row['nominal'] and len(displays) < 12:
                        displays.append((dict(row), profiles, raw_profiles))
                    rows.append(row)
        require(metadata(path) == before, 'Input file changed during investigation')
        if sample == 'data' or (index+1) % 20 == 0:
            print(f'{beam}GeV {sample}: file {index+1}/{len(files)}, sampled {len(rows)}', flush=True)
    # Uniform columns, also for events removed by the pre-existing tail veto.
    fields = sorted(set().union(*(r.keys() for r in rows)))
    rows = [{key: row.get(key, np.nan) for key in fields} for row in rows]
    path = out / sample / 'diagnostics' / f'{beam}GeV'
    path.mkdir(parents=True, exist_ok=True)
    columns = {key: np.array([row[key] for row in rows]) for key in fields}
    np.savez_compressed(path/'events.npz', **columns)
    (input_dir/'manifest.json').write_text(json.dumps(records, indent=2)+'\n')
    shutil.copy2(manifest, input_dir/'original_tail_manifest.tsv')
    (input_dir/'tail_source.json').write_text(json.dumps(metadata(source/'tail_study.root'), indent=2)+'\n')
    return columns, displays


def widths(values):
    if len(values) < 2 or np.mean(values) <= 0:
        return np.array([np.nan, np.nan])
    q16, median, q84 = np.quantile(values, [.16, .5, .84])
    return np.array([np.std(values)/np.mean(values), (q84-q16)/(2*median) if median > 0 else np.nan])


def gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-.5*((x-mean)/sigma)**2)


def fit_core(values, edges, window=None):
    if len(values) < 100:
        return dict(sigma_mean=np.nan, error=np.nan, lo=np.nan, hi=np.nan)
    counts, _ = np.histogram(values, edges)
    x = .5*(edges[1:]+edges[:-1])
    refine = window is None
    if refine:
        search = x >= .5*np.median(values)
        peak = x[search][np.argmax(counts[search])]
        window = (.75*peak, 1.25*peak)
    try:
        for iteration in range(2 if refine else 1):
            use = (x >= window[0]) & (x <= window[1])
            require(use.sum() >= 6, 'Too few bins in Gaussian core')
            guess = [max(counts[use]), np.mean(window), (window[1]-window[0])/5]
            with warnings.catch_warnings():
                warnings.simplefilter('error', OptimizeWarning)
                parameters, cov = curve_fit(gaussian, x[use], counts[use], p0=guess,
                    sigma=np.sqrt(np.maximum(counts[use], 1)), absolute_sigma=True,
                    bounds=([0, window[0], (edges[1]-edges[0])/2],
                            [np.inf, window[1], window[1]-window[0]]), maxfev=3000)
            if refine and iteration == 0:
                window = (parameters[1]-1.5*parameters[2], parameters[1]+1.5*parameters[2])
        _, mean, sigma = parameters
        gradient = np.array([0, -sigma/mean**2, 1/mean])
        return dict(sigma_mean=float(sigma/mean), error=float(np.sqrt(max(0, gradient@cov@gradient))),
                    lo=float(window[0]), hi=float(window[1]))
    except (RuntimeError, ValueError, OptimizeWarning):
        return dict(sigma_mean=np.nan, error=np.nan, lo=np.nan, hi=np.nan)


def summarize(beam, sample, columns, out, bootstrap):
    rows, changes, rates, run_rates = [], [], [], []
    energy = columns['energy_GeV']
    base = (~columns['tail_rich'].astype(bool)) & np.isfinite(energy)
    for cog in COGS:
        selected = base & (columns['cog_extent_mm'] < cog)
        before = energy[selected]
        if not len(before):
            continue
        edges = np.linspace(0, max(np.quantile(before, .995)*1.1, .001), 201)
        fit_before = fit_core(before, edges)
        fixed_window = (fit_before['lo'], fit_before['hi'])
        for name in CONFIGS:
            tag_mask = columns[name][selected].astype(bool)
            kept = before[~tag_mask]
            n, tagged = len(before), int(tag_mask.sum())
            rate = tagged/n
            rates.append(dict(energy_GeV=beam, sample=sample, cog_mm=cog, condition=name,
                              before=n, removed=tagged, removed_fraction=rate,
                              binomial_se=np.sqrt(rate*(1-rate)/n),
                              evaluable=int(((columns['front_active_x'][selected] >= 3) |
                                             (columns['front_active_y'][selected] >= 3)).sum()),
                              both_views=int((columns['nominal_x'][selected].astype(bool) &
                                              columns['nominal_y'][selected].astype(bool)).sum())))
            for stage, mask in [('before', np.ones(n, dtype=bool)), ('after', ~tag_mask), ('removed', tag_mask)]:
                if stage == 'before' and name != 'nominal':
                    continue
                values = before[mask]
                if not len(values):
                    rms, w68, mean, median = [np.nan]*4
                else:
                    rms, w68 = widths(values)
                    mean, median = np.mean(values), np.median(values)
                fit = fit_before if stage == 'before' else fit_core(values, edges, fixed_window)
                row = dict(energy_GeV=beam, sample=sample, cog_mm=cog, condition=name, stage=stage,
                           n=len(values), mean_GeV=mean, median_GeV=median,
                           rms_mean=rms, w68_median=w68,
                           gaussian_core=fit['sigma_mean'], gaussian_error=fit['error'],
                           fit_lo_GeV=fit['lo'], fit_hi_GeV=fit['hi'])
                rows.append(row)
                location = out / sample / (stage if name == 'nominal' else name+'_'+stage) / f'{cog}mm' / f'{beam}GeV'
                write_csv(location/'metrics.csv', [row])
                counts, _ = np.histogram(values, edges)
                np.savez_compressed(location/'energy_histogram.npz', edges_GeV=edges, counts=counts,
                                    underflow=(values < edges[0]).sum(), overflow=(values >= edges[-1]).sum())
            if name != 'nominal':
                continue
            delta = widths(kept)-widths(before)
            rng = np.random.default_rng(beam*100+cog+(sample == 'mc'))
            differences = []
            for _ in range(bootstrap):
                indices = rng.integers(0, n, n)
                values, flags = before[indices], tag_mask[indices]
                differences.append(widths(values[~flags])-widths(values))
            interval = np.nanquantile(differences, [.025, .975], axis=0)
            changes.append(dict(energy_GeV=beam, sample=sample, cog_mm=cog,
                delta_rms_mean=delta[0], delta_rms_ci_low=interval[0, 0], delta_rms_ci_high=interval[1, 0],
                delta_w68=delta[1], delta_w68_ci_low=interval[0, 1], delta_w68_ci_high=interval[1, 1],
                bootstrap=bootstrap))
            for run in np.unique(columns['run'][selected]):
                use = selected & (columns['run'] == run)
                run_rates.append(dict(energy_GeV=beam, sample=sample, cog_mm=cog, run=int(run),
                    n=int(use.sum()), tagged=int(columns['nominal'][use].astype(bool).sum())))
    return rows, rates, changes, run_rates


def plots(beam, sample, columns, displays, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    energy = columns['energy_GeV']
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for axis, cog in zip(axes.flat, COGS):
        selected = (~columns['tail_rich'].astype(bool)) & (columns['cog_extent_mm'] < cog)
        before = energy[selected]
        tagged = columns['nominal'][selected].astype(bool)
        if not len(before):
            continue
        edges = np.linspace(0, np.quantile(before, .999)*1.15, 141)
        for values, label, color, linestyle in [(before, 'Before', 'black', '-'),
                (before[~tagged], 'After veto', 'tab:blue', '--'),
                (before[tagged], 'Removed', 'tab:red', '-')]:
            axis.hist(values, bins=edges, histtype='step', label=f'{label} (n={len(values)})', color=color, linestyle=linestyle)
        axis.set_yscale('log')
        axis.set_ylim(bottom=.7)
        axis.set_title(('nocut' if cog == 200 else f'CoG {cog} mm') + f'; tag {tagged.mean():.2%}')
        axis.set_xlabel('Calibrated signal sum [GeV]')
        axis.set_ylabel('Events / bin')
        axis.legend(fontsize=8)
    fig.suptitle(f'{beam} GeV {sample}: persistent separated entrance components; existing HG-LG tail veto')
    fig.tight_layout()
    location = out / sample / 'diagnostics' / f'{beam}GeV'
    fig.savefig(location/'energy_before_after.png', dpi=150)
    plt.close(fig)
    for index, (row, profiles, raw) in enumerate(displays[:6]):
        fig, axes = plt.subplots(4, 2, figsize=(11, 9), sharex=True)
        for layer, axis in enumerate(axes.flat):
            axis.bar(*profiles[layer], width=4.7, color='tab:blue')
            axis.bar(*raw[layer], width=4.7, fill=False, edgecolor='0.5', linewidth=.8, label='Without mask')
            axis.set_title(f'Layer {layer}: '+('x' if layer % 2 else 'y'))
            axis.set_xlim(-120, 120)
            axis.set_ylabel('Signal [MeV]')
            view = 'x' if layer % 2 else 'y'
            for key in (view+'_axis1_mm', view+'_axis2_mm'):
                if np.isfinite(row[key]):
                    axis.axvline(row[key], color='tab:red', linestyle='--')
        for axis in axes[-1]:
            axis.set_xlabel('Strip short-coordinate [mm]')
        fig.suptitle(f'{beam} GeV {sample}; file {row["file_index"]}, entry {row["entry"]}, run {row["run"]}; '
                     f'signal {row["energy_GeV"]:.4g} GeV')
        fig.tight_layout()
        fig.savefig(location/f'candidate_{index}.png', dpi=130)
        plt.close(fig)


def load_displays(beam, sample, columns, out):
    """Re-read a few candidate entries for plots; never rerun event selection."""
    records = json.loads((out/'inputs'/f'{beam}GeV'/sample/'manifest.json').read_text())
    text = SOURCE.read_text().split('const unordered_set<int> bad_channels = {', 1)[1].split('};', 1)[0]
    bad = np.array(sorted(set(map(int, re.findall(r'\d+', re.sub(r'//[^\n]*', '', text))))))
    candidates = np.flatnonzero((~columns['tail_rich'].astype(bool)) & columns['nominal'].astype(bool))[:6]
    output = []
    for index in candidates:
        row = {key: value[index].item() for key, value in columns.items()}
        path = Path(records[int(row['file_index'])]['path'])
        entry = int(row['entry'])
        with uproot.open(path) as root:
            a = root['Calib_Hit'].arrays(['CellID', 'Hit_Energy', 'Hit_X', 'Hit_Y'],
                                       entry_start=entry, entry_stop=entry+1, library='ak')
        ids, en, x, y = (ak.to_numpy(a[key][0]) for key in a.fields)
        layer, chip, channel = ids//100000, ids%100000//10000, ids%100
        connected = (ids >= 0) & (layer < 30) & (chip < 6) & (channel < 36) & ~((chip == 5) & (channel > 29))
        good = connected & (en >= 0) & ~np.isin(ids, bad)
        output.append((row, projections(ids, en, x, y, good),
                       projections(ids, en, x, y, connected & (en >= 0))))
    return output


def apply_mc_audit(beam, columns, out):
    path = out/'inputs'/f'{beam}GeV'/'current_mc_tail_audit.csv'
    if not path.exists():
        return False
    with path.open() as source:
        rows = list(csv.DictReader(source))
    lookup = {(int(row['file_index']), int(row['entry'])): row for row in rows}
    require(len(lookup) == len(columns['entry']), 'MC tail audit does not cover the sampled events')
    flags = []
    for index, key in enumerate(zip(columns['file_index'], columns['entry'])):
        row = lookup[key]
        require(int(row['run']) == int(columns['run'][index]) and
                int(row['event_time']) == int(columns['event_time'][index]) and
                int(row['event']) == int(columns['event'][index]), 'MC audit identity mismatch')
        flags.append(bool(int(row['tail_rich'])))
    columns['tail_rich'] = np.array(flags, dtype=bool)
    require(np.isfinite(columns['energy_GeV']).all(), 'MC feature coverage incomplete')
    np.savez_compressed(out/'mc/diagnostics'/f'{beam}GeV'/'events.npz', **columns)
    return True


def combine(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for beam in (5, 40, 100):
        status = json.loads((out/'inputs'/f'{beam}GeV'/'complete.json').read_text())
        require(status.get('current_mc_tail_audited'), f'{beam}GeV current MC tail audit is pending')
    tables = {}
    for name in ('metrics', 'tag_rates', 'paired_changes', 'run_rates'):
        rows = []
        for beam in (5, 40, 100):
            with (out/'comparison'/f'{beam}GeV'/(name+'.csv')).open() as source:
                rows.extend(csv.DictReader(source))
        tables[name] = rows
        write_csv(out/'comparison'/(name+'.csv'), rows)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for axis, beam in zip(axes, (5, 40, 100)):
        for sample, color in (('data', 'black'), ('mc', 'tab:red')):
            for stage, marker, linestyle in (('before', 'o', '-'), ('after', 's', '--')):
                rows = [r for r in tables['metrics'] if int(r['energy_GeV']) == beam and
                        r['sample'] == sample and r['condition'] == 'nominal' and r['stage'] == stage]
                rows.sort(key=lambda r: COGS.index(int(r['cog_mm'])))
                axis.plot(range(4), [float(r['w68_median'])*100 for r in rows],
                          color=color, marker=marker, linestyle=linestyle, label=f'{sample} {stage}')
        axis.set_xticks(range(4), ['nocut', '20', '10', '5'])
        axis.set_title(f'{beam} GeV')
        axis.set_xlabel('CoG half-width [mm]')
        axis.set_ylabel('Central 68% half-width / median [%]')
        axis.legend(fontsize=8)
        axis.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(out/'comparison'/'width_summary.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energy', type=int, choices=(5, 40, 100))
    parser.add_argument('--events', type=int, default=20000, help='per sample before existing selections; zero=all')
    parser.add_argument('--bootstrap', type=int, default=400)
    parser.add_argument('--out', type=Path, default=RESULT)
    parser.add_argument('--summarize-only', action='store_true')
    parser.add_argument('--combine', action='store_true')
    args = parser.parse_args()
    require(args.events >= 0 and args.bootstrap >= 50, 'Invalid sampling/bootstrap options')
    out = args.out.resolve()
    # Avoid fontconfig scanning remote software trees during this small investigation.
    font_cache = tempfile.TemporaryDirectory(prefix='scecal-veto-mpl-')
    import atexit
    atexit.register(font_cache.cleanup)
    fonts = Path(importlib.util.find_spec('matplotlib').origin).parent/'mpl-data/fonts/ttf'
    fontconfig = Path(font_cache.name)/'fonts.conf'
    fontconfig.write_text(f'<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd">'
                         f'<fontconfig><dir>{fonts}</dir><cachedir>{font_cache.name}</cachedir></fontconfig>')
    os.environ['FONTCONFIG_FILE'] = str(fontconfig)
    os.environ['MPLCONFIGDIR'] = font_cache.name
    if args.combine:
        combine(out)
        return
    require(args.energy is not None, 'Specify --energy or --combine')
    record = out/'inputs'/f'{args.energy}GeV'
    record.mkdir(parents=True, exist_ok=True)
    if args.summarize_only:
        shutil.copy2(Path(__file__), record/'summary_source.py')
        (record/'summary_configuration.json').write_text(json.dumps(dict(
            command=sys.argv, bootstrap=args.bootstrap,
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')
    if not args.summarize_only:
        require(not (record/'complete.json').exists(), 'Results already complete; use summarize-only or another output base')
        shutil.copy2(Path(__file__), record/'check.py')
        shutil.copy2(SOURCE, record/'resolution_source.cc')
        config = dict(energy=args.energy, target_per_sample=args.events, bootstrap=args.bootstrap,
            seed=20261111, pitch_mm=PITCH, MIP_MeV=MIP, configs=CONFIGS,
            layer_range=[0, 7], projected_bin_threshold_MIP=.5, component_minimum_MIP=1,
            matching_tolerance_mm=PITCH, mask_confirmation='same two axes with and without bad-channel mask',
            existing_selection='HG-LG tail-rich excluded; nonnegative good-hit sum; layer9/10 square CoG',
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            source_mask_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            python=sys.version, numpy=np.__version__, uproot=uproot.__version__,
            caveat='Topology candidates, not identified multiple particles; current source cannot prove historical production')
        (record/'configuration.json').write_text(json.dumps(config, indent=2)+'\n')
        (record/'command.txt').write_text(' '.join(sys.argv)+'\n')
    combined = [[], [], [], []]
    with tempfile.TemporaryDirectory(prefix='scecal-veto-mpl-') as cache:
        os.environ['MPLCONFIGDIR'] = cache
        for sample in ('data', 'mc'):
            if args.summarize_only:
                with np.load(out/sample/'diagnostics'/f'{args.energy}GeV'/'events.npz') as saved:
                    columns = {key: saved[key] for key in saved.files}
                displays = load_displays(args.energy, sample, columns, out)
            else:
                columns, displays = collect(args.energy, sample, args.events, out)
            if sample == 'mc':
                audit_complete = apply_mc_audit(args.energy, columns, out)
            summaries = summarize(args.energy, sample, columns, out, args.bootstrap)
            for rows, more in zip(combined, summaries):
                rows.extend(more)
            plots(args.energy, sample, columns, displays, out)
    for name, rows in zip(('metrics', 'tag_rates', 'paired_changes', 'run_rates'), combined):
        write_csv(out/'comparison'/f'{args.energy}GeV'/(name+'.csv'), rows)
    (record/'complete.json').write_text(json.dumps(dict(completed=True, current_mc_tail_audited=audit_complete,
                                                      unix_time=time.time()))+'\n')
    for row in combined[1]:
        if row['condition'] == 'nominal':
            print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
