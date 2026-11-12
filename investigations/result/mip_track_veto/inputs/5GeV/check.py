#!/usr/bin/env python3
"""Estimate a sparse through-going track veto on the previous investigation's events."""
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

import awkward as ak
import numpy as np
import uproot

sys.dont_write_bytecode = True
PROJECT = Path(__file__).resolve().parents[2]
PREVIOUS = PROJECT/'investigations/result/multiple_shower_veto'
RESULT = PROJECT/'investigations/result/mip_track_veto'
HELPER_PATH = PROJECT/'investigations/multiple_shower_veto/check.py'
spec = importlib.util.spec_from_file_location('previous_veto_investigation', HELPER_PATH)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
require, metadata, write_csv = helper.require, helper.metadata, helper.write_csv
COGS = (200, 20, 10, 5)
CONFIGS = {'nominal': (15, 5.3), 'layers20': (20, 5.3), 'rms10p6': (15, 10.6)}
LAYER_Z = np.array([1+(l//2)*19.9 if l % 2 == 0 else 12.2+(l//2)*19.9 for l in range(30)])


def straightness(energy, q, width, view):
    layers = np.arange(view, 30, 2)
    layers = layers[energy[layers] > 0]
    if len(layers) < 4:
        return np.nan, np.nan, np.nan
    z, position = LAYER_Z[layers], q[layers]
    dz = z-z.mean()
    slope = np.dot(dz, position-position.mean())/np.dot(dz, dz)
    intercept = position.mean()-slope*z.mean()
    residual = position-(intercept+slope*z)
    return np.sqrt(np.mean(residual**2+width[layers]**2)), intercept, slope


def moments(ids, energy, x, y, z, bad):
    layer, chip, channel = ids//100000, ids%100000//10000, ids%100
    connected = (ids >= 0)&(layer < 30)&(chip < 6)&(channel < 36)&~((chip == 5)&(channel > 29))
    good = connected&(energy >= 0)&~np.isin(ids, bad)
    require(np.isfinite(energy).all() and np.isfinite(x).all() and np.isfinite(y).all(), 'Nonfinite hit')
    require(np.allclose(z[connected], LAYER_Z[layer[connected]], atol=.01, rtol=0), 'Layer geometry differs')
    l, e = layer[good], energy[good]
    sums = np.bincount(l, weights=e, minlength=30)
    divide = lambda values: np.divide(values, sums, out=np.full(30, np.nan), where=sums > 0)
    qx = divide(np.bincount(l, weights=e*x[good], minlength=30))
    qy = divide(np.bincount(l, weights=e*y[good], minlength=30))
    wx = np.sqrt(np.maximum(0, divide(np.bincount(l, weights=e*x[good]**2, minlength=30))-qx*qx))
    wy = np.sqrt(np.maximum(0, divide(np.bincount(l, weights=e*y[good]**2, minlength=30))-qy*qy))
    rx, ax, bx = straightness(sums, qx, wx, 1)
    ry, ay, by = straightness(sums, qy, wy, 0)
    active = sums > 0
    physical = helper.physical(ids)
    physical_good = connected&(energy >= 0)&~np.isin(physical, bad)
    return dict(hits=int(good.sum()), active=int(active.sum()), front=int(active[:5].sum()),
        back=int(active[25:].sum()), occupancy=float(good.sum()/max(active.sum(), 1)),
        rms_x_mm=rx, rms_y_mm=ry, intercept_x_mm=ax, intercept_y_mm=ay, slope_x=bx, slope_y=by,
        layer_energy_MeV=sums, cog_x_mm=qx, cog_y_mm=qy, width_x_mm=wx, width_y_mm=wy,
        recomputed_energy_GeV=float(e.sum()/1000),
        physical_mask_energy_GeV=float(energy[physical_good].sum()/1000),
        memory_hits=int((physical[connected] != ids[connected]).sum()))


def collect(beam, sample, out):
    prior_record = PREVIOUS/'inputs'/f'{beam}GeV'
    require(json.loads((prior_record/'complete.json').read_text()).get('current_mc_tail_audited'), 'Prior MC audit incomplete')
    baseline_path = PREVIOUS/sample/'diagnostics'/f'{beam}GeV'/'events.npz'
    with np.load(baseline_path) as saved:
        a = {k: saved[k].copy() for k in ('file_index', 'entry', 'run', 'event_time', 'event',
                                        'tail_rich', 'energy_GeV', 'cog_extent_mm')}
    records = json.loads((prior_record/sample/'manifest.json').read_text())
    source = prior_record/'resolution_source.cc'
    text = source.read_text().split('const unordered_set<int> bad_channels = {', 1)[1].split('};', 1)[0]
    bad = np.array(sorted(set(map(int, re.findall(r'\d+', re.sub(r'//[^\n]*', '', text))))))
    n = len(a['entry'])
    # Other events cannot satisfy this veto and do not need a ROOT re-read.
    eligible = (~a['tail_rich'].astype(bool))&(a['energy_GeV'] > 0)&(a['energy_GeV'] <= .03)&(a['cog_extent_mm'] < 200)
    a['evaluated'] = eligible
    for key in ('hits', 'active', 'front', 'back', 'occupancy', 'rms_x_mm', 'rms_y_mm',
                'intercept_x_mm', 'intercept_y_mm', 'slope_x', 'slope_y', 'physical_mask_energy_GeV', 'memory_hits'):
        a[key] = np.full(n, np.nan)
    for key in ('layer_energy_MeV', 'cog_x_mm', 'cog_y_mm', 'width_x_mm', 'width_y_mm'):
        a[key] = np.full((n, 30), np.nan)
    max_energy_difference = 0.
    read_count = 0
    for index, record in enumerate(records):
        path = Path(record['path'])
        require(metadata(path) == {k: record[k] for k in ('path', 'bytes', 'mtime_ns')}, 'Input changed since previous investigation')
        selected = np.flatnonzero(eligible&(a['file_index'] == index))
        if not len(selected):
            continue
        with uproot.open(path, array_cache=None) as root:
            tree = root['Calib_Hit']
            require(tree.num_entries == record['entries'], 'Entry count changed')
            for start, stop in record['spans']:
                selected_rows = selected[(a['entry'][selected] >= start)&(a['entry'][selected] < stop)]
                if not len(selected_rows):
                    continue
                lo, hi = int(a['entry'][selected_rows].min()), int(a['entry'][selected_rows].max())+1
                arrays = tree.arrays(['Run_Num', 'Event_Time', 'Event_Num', 'CellID', 'Hit_Energy',
                                     'Hit_X', 'Hit_Y', 'Hit_Z'], entry_start=lo, entry_stop=hi, library='ak')
                for row in selected_rows:
                    j = int(a['entry'][row])-lo
                    require(tuple(int(arrays[k][j]) for k in ('Run_Num', 'Event_Time', 'Event_Num')) ==
                            tuple(int(a[k][row]) for k in ('run', 'event_time', 'event')), 'Event identity differs')
                    inputs = [ak.to_numpy(arrays[k][j]) for k in ('CellID', 'Hit_Energy', 'Hit_X', 'Hit_Y', 'Hit_Z')]
                    require(len(set(map(len, inputs))) == 1, 'Hit vector lengths differ')
                    values = moments(*inputs, bad)
                    difference = abs(values.pop('recomputed_energy_GeV')-a['energy_GeV'][row])
                    max_energy_difference = max(max_energy_difference, difference)
                    require(difference < 1e-10, 'Current signal differs from previous baseline')
                    for key, value in values.items():
                        a[key][row] = value
                    read_count += 1
        require(metadata(path) == {k: record[k] for k in ('path', 'bytes', 'mtime_ns')}, 'Input changed during investigation')
        print(f'{beam}GeV {sample}: {read_count}/{eligible.sum()} possible track entries checked', flush=True)
    require(read_count == int(eligible.sum()), 'Incomplete topology coverage')
    for name, (layers, rms) in CONFIGS.items():
        a[name] = eligible&(a['hits'] <= 60)&(a['active'] >= layers)&(a['front'] >= 2)&(a['back'] >= 2)&\
                  (a['occupancy'] <= 2)&(a['rms_x_mm'] <= rms)&(a['rms_y_mm'] <= rms)
    directory = out/sample/'diagnostics'/f'{beam}GeV'
    directory.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(directory/'events.npz', **a)
    record_out = out/'inputs'/f'{beam}GeV'/sample
    record_out.mkdir(parents=True, exist_ok=True)
    shutil.copy2(prior_record/sample/'manifest.json', record_out/'manifest.json')
    (record_out/'validation.json').write_text(json.dumps(dict(
        sampled=n, possible_entries=int(eligible.sum()), reread=read_count,
        max_abs_energy_difference_GeV=max_energy_difference,
        previous_features=metadata(baseline_path), previous_features_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        mask_source=metadata(source), mask_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        identity='file/entry and Run/Time/Event checked for every potentially vetoed entry; unchanged input size/mtime for all files',
        mask_definition='Original full CellID lookup, same as previous investigation/current resolution; not the physical-ID lookup of shower_leak Aggregate.cc',
        track_geometry='TrackSelection.hh v1 geometry; z checked against current EBUdecode; accepted zero-energy and repeated memory hits counted separately'), indent=2)+'\n')
    return a


def summarize(beam, sample, a, out, bootstrap):
    metrics, rates, changes = [], [], []
    for cog in COGS:
        use = (~a['tail_rich'].astype(bool))&np.isfinite(a['energy_GeV'])&(a['cog_extent_mm'] < cog)
        before = a['energy_GeV'][use]
        edges = np.linspace(0, max(np.quantile(before, .995)*1.1, .001), 201)
        core = helper.fit_core(before, edges)
        window = (core['lo'], core['hi'])
        for name in CONFIGS:
            rejected = a[name][use]
            n, removed = len(before), int(rejected.sum())
            rates.append(dict(energy_GeV=beam, sample=sample, cog_mm=cog, condition=name,
                before=n, removed=removed, removed_fraction=removed/n,
                low_signal_before=int((before <= .03).sum()), low_signal_after=int(((before <= .03)&~rejected).sum())))
            for stage, keep in (('before', np.ones(n, dtype=bool)), ('after', ~rejected), ('removed', rejected)):
                if stage == 'before' and name != 'nominal':
                    continue
                v = before[keep]
                rms, w68 = helper.widths(v)
                fit = core if stage == 'before' else helper.fit_core(v, edges, window)
                counts, _ = np.histogram(v, edges)
                row = dict(energy_GeV=beam, sample=sample, cog_mm=cog, condition=name, stage=stage,
                    n=len(v), mean_GeV=float(v.mean()) if len(v) else np.nan,
                    rms_mean=rms, w68_median=w68, gaussian_core=fit['sigma_mean'],
                    gaussian_error=fit['error'], fit_lo_GeV=fit['lo'], fit_hi_GeV=fit['hi'])
                metrics.append(row)
                location = out/sample/(stage if name == 'nominal' else name+'_'+stage)/f'{cog}mm'/f'{beam}GeV'
                write_csv(location/'metrics.csv', [row])
                np.savez_compressed(location/'energy_histogram.npz', edges_GeV=edges, counts=counts,
                                    underflow=(v < edges[0]).sum(), overflow=(v > edges[-1]).sum())
            if name == 'nominal':
                rng = np.random.default_rng(beam*100+cog+(sample == 'mc'))
                diffs = []
                for _ in range(bootstrap):
                    indices = rng.integers(0, n, n)
                    values, flags = before[indices], rejected[indices]
                    diffs.append(helper.widths(values[~flags])-helper.widths(values))
                interval = np.nanquantile(diffs, [.025, .975], axis=0)
                delta = helper.widths(before[~rejected])-helper.widths(before)
                changes.append(dict(energy_GeV=beam, sample=sample, cog_mm=cog,
                    delta_rms_mean=delta[0], delta_rms_low=interval[0, 0], delta_rms_high=interval[1, 0],
                    delta_w68=delta[1], delta_w68_low=interval[0, 1], delta_w68_high=interval[1, 1], bootstrap=bootstrap))
    return metrics, rates, changes


def plots(beam, sample, a, out):
    import matplotlib.pyplot as plt
    location = out/sample/'diagnostics'/f'{beam}GeV'
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for axis, cog in zip(axes.flat, COGS):
        selected = (~a['tail_rich'].astype(bool))&(a['cog_extent_mm'] < cog)
        values, rejected = a['energy_GeV'][selected], a['nominal'][selected]
        edges = np.linspace(0, max(np.quantile(values, .999)*1.15, .031), 151)
        for v, label, color, style in ((values, 'Before', 'black', '-'),
                (values[~rejected], 'After veto', 'tab:blue', '--'), (values[rejected], 'Removed', 'tab:red', '-')):
            axis.hist(v, bins=edges, histtype='step', color=color, linestyle=style, label=f'{label} (n={len(v)})')
        axis.set(yscale='log', ylim=(.7, None), xlabel='Calibrated signal sum [GeV]', ylabel='Events / bin',
                 title=('nocut' if cog == 200 else f'CoG {cog} mm')+f'; removed {rejected.mean():.2%}')
        axis.legend(fontsize=8)
    fig.suptitle(f'{beam} GeV {sample}: sparse through-going MIP-like veto')
    fig.tight_layout()
    fig.savefig(location/'energy_before_after.png', dpi=140)
    plt.close(fig)
    candidates = np.flatnonzero(a['nominal']&(a['cog_extent_mm'] < 20))[:3]
    for i, row in enumerate(candidates):
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        for axis, view, coordinate in ((axes[0], 1, 'x'), (axes[1], 0, 'y')):
            layers = np.arange(view, 30, 2)
            layers = layers[a['layer_energy_MeV'][row, layers] > 0]
            axis.errorbar(LAYER_Z[layers], a['cog_'+coordinate+'_mm'][row, layers],
                          yerr=a['width_'+coordinate+'_mm'][row, layers], fmt='o', markersize=4)
            axis.plot(LAYER_Z, a['intercept_'+coordinate+'_mm'][row]+a['slope_'+coordinate][row]*LAYER_Z, '--')
            axis.set(xlabel='Layer z [mm]', ylabel=coordinate+' [mm]', title=f'{coordinate}: RMS {a["rms_"+coordinate+"_mm"][row]:.2f} mm')
        axes[2].bar(range(30), a['layer_energy_MeV'][row]/.305)
        axes[2].set(xlabel='Layer', ylabel='Signal [MIP]', title=f'{a["active"][row]:.0f} active layers; {a["hits"][row]:.0f} hits')
        fig.suptitle(f'{beam} GeV {sample}: run {a["run"][row]}, entry {a["entry"][row]}; signal {a["energy_GeV"][row]*1000:.2f} MeV')
        fig.tight_layout()
        fig.savefig(location/f'candidate_{i}.png', dpi=140)
        plt.close(fig)


def combine(out):
    tables = {}
    for name in ('metrics', 'tag_rates', 'paired_changes'):
        rows = []
        for beam in (5, 40, 100):
            require(json.loads((out/'inputs'/f'{beam}GeV'/'complete.json').read_text())['complete'], 'Incomplete energy')
            with (out/'comparison'/f'{beam}GeV'/(name+'.csv')).open() as source:
                rows.extend(csv.DictReader(source))
        tables[name] = rows
        write_csv(out/'comparison'/(name+'.csv'), rows)
    summary = []
    for beam in (5, 40, 100):
        m = {r['stage']: r for r in tables['metrics'] if int(r['energy_GeV']) == beam and r['sample'] == 'data' and r['cog_mm'] == '20' and r['condition'] == 'nominal'}
        rates = {r['sample']: r for r in tables['tag_rates'] if int(r['energy_GeV']) == beam and r['cog_mm'] == '20' and r['condition'] == 'nominal'}
        row = dict(energy_GeV=beam, data_before=rates['data']['before'], data_removed=rates['data']['removed'],
                   data_removed_percent=float(rates['data']['removed_fraction'])*100,
                   mc_before=rates['mc']['before'], mc_removed=rates['mc']['removed'])
        for key in ('rms_mean', 'w68_median', 'gaussian_core'):
            for stage in ('before', 'after'):
                row[key+'_'+stage+'_percent'] = float(m[stage][key])*100
        summary.append(row)
    write_csv(out/'comparison/summary20mm.csv', summary)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for axis, beam in zip(axes, (5, 40, 100)):
        for sample, color in (('data', 'black'), ('mc', 'tab:red')):
            for stage, style, marker in (('before', '-', 'o'), ('after', '--', 's')):
                rows = [r for r in tables['metrics'] if int(r['energy_GeV']) == beam and r['sample'] == sample and r['condition'] == 'nominal' and r['stage'] == stage]
                rows.sort(key=lambda r: COGS.index(int(r['cog_mm'])))
                axis.plot(range(4), [float(r['rms_mean'])*100 for r in rows], color=color,
                          linestyle=style, marker=marker, label=f'{sample} {stage}')
        axis.set(xticks=range(4), xticklabels=['nocut', '20', '10', '5'], xlabel='CoG half-width [mm]',
                 ylabel='Full distribution RMS / mean [%]', title=f'{beam} GeV')
        axis.legend(fontsize=8)
        axis.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(out/'comparison/rms_summary.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energy', type=int, choices=(5, 40, 100))
    parser.add_argument('--combine', action='store_true')
    parser.add_argument('--out', type=Path, default=RESULT)
    parser.add_argument('--bootstrap', type=int, default=400)
    args = parser.parse_args()
    out = args.out.resolve()
    record = out/'inputs'/('combined' if args.combine else f'{args.energy}GeV')
    record.mkdir(parents=True, exist_ok=True)
    require(not (record/'complete.json').exists(), 'Completed output exists; choose another --out')
    require(args.combine or args.energy is not None, 'Specify --energy or --combine')
    require(args.bootstrap >= 50, 'Too few bootstrap samples')
    for source in (Path(__file__), HELPER_PATH, PROJECT/'analysis/shower_leak/TrackSelection.hh', PROJECT/'analysis/calibration_parameters/EBUdecode.cxx'):
        destination = record/('width_helper.py' if source == HELPER_PATH else source.name)
        shutil.copy2(source, destination)
    (record/'configuration.json').write_text(json.dumps(dict(command=sys.argv,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        helper_sha256=hashlib.sha256(HELPER_PATH.read_bytes()).hexdigest(), configs=CONFIGS,
        max_signal_GeV=.03, max_hits=60, max_hits_per_active_layer=2,
        min_active_front5=2, min_active_back5=2, previous_result=str(PREVIOUS),
        mask='Same original full CellID mask and nonnegative signal as previous investigation',
        note='Topology tag, not muon/pion identification; energy ceiling affects full width but need not affect the peak core'), indent=2)+'\n')
    with tempfile.TemporaryDirectory(prefix='scecal-mip-veto-') as cache:
        fonts = Path(importlib.util.find_spec('matplotlib').origin).parent/'mpl-data/fonts/ttf'
        config = Path(cache)/'fonts.conf'
        config.write_text(f'<fontconfig><dir>{fonts}</dir><cachedir>{cache}</cachedir></fontconfig>')
        os.environ['FONTCONFIG_FILE'], os.environ['MPLCONFIGDIR'] = str(config), cache
        import matplotlib
        matplotlib.use('Agg')
        if args.combine:
            combine(out)
        else:
            tables = [[], [], []]
            for sample in ('data', 'mc'):
                a = collect(args.energy, sample, out)
                for rows, more in zip(tables, summarize(args.energy, sample, a, out, args.bootstrap)):
                    rows.extend(more)
                plots(args.energy, sample, a, out)
            for name, rows in zip(('metrics', 'tag_rates', 'paired_changes'), tables):
                write_csv(out/'comparison'/f'{args.energy}GeV'/(name+'.csv'), rows)
            for row in tables[1]:
                if row['condition'] == 'nominal':
                    print(json.dumps(row), flush=True)
    (record/'complete.json').write_text(json.dumps(dict(complete=True))+'\n')


if __name__ == '__main__':
    main()
