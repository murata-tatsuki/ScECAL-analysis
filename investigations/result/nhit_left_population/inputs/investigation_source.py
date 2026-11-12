#!/usr/bin/env python3
"""Small data-only investigation of the low-signal side of E versus Nhit."""
import importlib.util
import json
import os
import re
import sys
import tempfile
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

sys.dont_write_bytecode = True
PROJECT = Path(__file__).resolve().parents[2]
OUT = PROJECT/'investigations/result/nhit_left_population'
spec = importlib.util.spec_from_file_location('veto_helper', PROJECT/'investigations/multiple_shower_veto/check.py')
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
SOURCE = PROJECT/'analysis/resolution/SingleEnergyAnalysis.cc'
text = SOURCE.read_text().split('const unordered_set<int> bad_channels = {', 1)[1].split('};', 1)[0]
BAD = np.array(sorted(set(map(int, re.findall(r'\d+', re.sub(r'//[^\n]*', '', text))))))


def collect(beam, target=18000):
    source = PROJECT/f'analysis/result/channel_response/hg_lg_tail/threshold/e_ssa/{beam}GeV'
    with uproot.open(source/'tail_study.root') as root:
        a = root['events'].arrays(['sample', 'file_index', 'entry', 'Run_Num', 'Event_Time',
            'TriggerID', 'tail_rich', 'low_tail_channels', 'energy_GeV', 'good_hits'], library='np')
    a = {k: v[a['sample'] == 0] for k, v in a.items()}
    files = [helper.relocated(line.split('\t')[1]) for line in (source/'inputs/data.tsv').read_text().splitlines() if line]
    counts = [int((a['file_index'] == i).sum()) for i in range(len(files))]
    desired = target*np.array(counts)/sum(counts)
    takes = np.floor(desired).astype(int)
    for i in np.argsort(-(desired-takes))[:target-int(takes.sum())]:
        takes[i] += 1
    rng = np.random.default_rng(20261111+beam*10)
    records, rows = [], []
    max_energy_diff = 0.
    for i, (path, count, take) in enumerate(zip(files, counts, takes)):
        lookup = np.flatnonzero(a['file_index'] == i)
        lookup = lookup[np.argsort(a['entry'][lookup])]
        helper.require(np.array_equal(a['entry'][lookup], np.arange(count)), 'Incomplete tail flags')
        spans = helper.sample_blocks(count, int(take), rng)
        meta = helper.metadata(path)
        records.append(dict(file_index=i, entries=count, spans=spans, **meta))
        with uproot.open(path, array_cache=None) as root:
            tree = root['Calib_Hit']
            helper.require(tree.num_entries == count, 'Entry count differs from tail tree')
            for lo, hi in spans:
                hit = tree.arrays(['Run_Num', 'Event_Time', 'Event_Num', 'CellID',
                    'Hit_Energy', 'Hit_X', 'Hit_Y'], entry_start=lo, entry_stop=hi, library='ak')
                for j, entry in enumerate(range(lo, hi)):
                    k = lookup[entry]
                    ident = tuple(int(hit[key][j]) for key in ('Run_Num', 'Event_Time', 'Event_Num'))
                    helper.require(ident == tuple(int(a[key][k]) for key in ('Run_Num', 'Event_Time', 'TriggerID')), 'Event identity mismatch')
                    ids, en, x, y = [ak.to_numpy(hit[key][j]) for key in ('CellID', 'Hit_Energy', 'Hit_X', 'Hit_Y')]
                    layer = ids//100000
                    helper.require(np.all(ids >= 0) and np.isfinite(en).all(), 'Invalid calibrated hit')
                    good = (layer < 30)&(en >= 0)&~np.isin(ids, BAD)
                    l, e = layer[good], en[good]
                    q = np.where(l % 2, x[good], y[good])
                    energy = e.sum()
                    le = np.bincount(l, weights=e, minlength=30)
                    nh = np.bincount(l, minlength=30)
                    div = lambda v: np.divide(v, le, out=np.full(30, np.nan), where=le > 0)
                    cx = div(np.bincount(l, weights=e*x[good], minlength=30))
                    cy = div(np.bincount(l, weights=e*y[good], minlength=30))
                    cq = div(np.bincount(l, weights=e*q, minlength=30))
                    w2 = np.maximum(0, div(np.bincount(l, weights=e*q*q, minlength=30))-cq*cq)
                    cog = np.array([cx[9], cy[9], cx[10], cy[10]])
                    cells = np.unique(helper.physical(ids[good&(en > 0)]))
                    strips = np.bincount(cells//100000, minlength=30)
                    dense = (le >= 6*.305)&(strips >= 4)
                    pairs = dense[:-1]&dense[1:]
                    cached = float(a['energy_GeV'][k])*1000
                    max_energy_diff = max(max_energy_diff, abs(energy-cached))
                    helper.require(abs(energy-cached) < 1e-5 and good.sum() == a['good_hits'][k], 'Current data differ from tail cache')
                    row = dict(file_index=i, entry=entry, run=ident[0], event_time=ident[1], event=ident[2],
                        energy_GeV=energy/1000, hits=int(good.sum()), tail_rich=bool(a['tail_rich'][k]),
                        low_tail_channels=int(a['low_tail_channels'][k]),
                        cog_extent_mm=np.max(np.abs(cog)) if np.isfinite(cog).all() else np.nan,
                        shower_start=int(np.argmax(pairs)) if pairs.any() else 99,
                        depth=float(np.dot(le, np.arange(30))/energy) if energy > 0 else np.nan,
                        back_fraction=le[25:].sum()/energy if energy > 0 else np.nan,
                        front_fraction=le[:5].sum()/energy if energy > 0 else np.nan,
                        width_mm=np.sqrt(np.nansum(le*w2)/energy) if energy > 0 else np.nan,
                        edge_fraction=e[np.abs(q) > 70].sum()/energy if energy > 0 else np.nan,
                        zero_hits=int((e == 0).sum()), unique_strips=len(cells),
                        memory_fraction=float((helper.physical(ids[good]) != ids[good]).mean()) if good.any() else 0,
                        low_hit_fraction=float((e < .5*.305).mean()) if len(e) else np.nan,
                        layer_energy_MeV=le, layer_hits=nh, cog_x_mm=cx, cog_y_mm=cy, layer_width_mm=np.sqrt(w2))
                    rows.append(row)
        helper.require(helper.metadata(path) == meta, 'Input changed while reading')
        print(f'{beam} GeV file {i+1}/{len(files)}: {len(rows)} events; max signal difference {max_energy_diff:g} MeV', flush=True)
    location = OUT/'data/diagnostics'/f'{beam}GeV'
    location.mkdir(parents=True, exist_ok=True)
    columns = {key: np.array([r[key] for r in rows]) for key in rows[0]}
    np.savez_compressed(location/'events.npz', **columns)
    inp = OUT/'inputs'/f'{beam}GeV'
    inp.mkdir(parents=True, exist_ok=True)
    (inp/'manifest.json').write_text(json.dumps(records, indent=2)+'\n')
    (inp/'validation.json').write_text(json.dumps(dict(sampled=len(rows), max_signal_difference_MeV=max_energy_diff,
        identity='All sampled Run/Time/Event, energy, Nhit match tail cache; full file size/mtime checked',
        tail_source=helper.metadata(source/'tail_study.root'), mask_source=helper.metadata(SOURCE)), indent=2)+'\n')
    return columns


def summarize():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    metrics, features, rates, position = [], [], [], []
    for beam in (20, 30, 40):
        a = dict(np.load(OUT/'data/diagnostics'/f'{beam}GeV/events.npz'))
        e, n = a['energy_GeV'], a['hits']
        base = (~a['tail_rich'])&(a['cog_extent_mm'] < 200)
        central = base&(a['cog_extent_mm'] < 20)
        eref, nref = np.median(e[central]), np.median(n[central])
        ratio = e/np.maximum(n, 1)
        low = (n > .7*nref)&(n < 1.35*nref)&(ratio < .8*eref/nref)
        main = (n > .7*nref)&(n < 1.35*nref)&(ratio > .9*eref/nref)&(ratio < 1.1*eref/nref)
        # These exploratory cuts are diagnostics, not a validated pion classifier.
        cuts = {'left_band': low, 'width14': a['width_mm'] > 14,
                'late_shower': a['shower_start'] > 6}
        for cog in (200, 20, 10, 5):
            use = base&(a['cog_extent_mm'] < cog)
            values = e[use]
            edges = np.linspace(0, np.quantile(values, .995)*1.1, 201)
            initial = helper.fit_core(values, edges)
            window = (initial['lo'], initial['hi'])
            for condition, rejected in {'before': np.zeros(len(e), dtype=bool), **cuts}.items():
                keep = use&~rejected
                fit = initial if condition == 'before' else helper.fit_core(e[keep], edges, window)
                rms, w68 = helper.widths(e[keep])
                row = dict(energy_GeV=beam, cog_mm=cog, condition=condition,
                    n=int(keep.sum()), removed=int((use&rejected).sum()),
                    removed_fraction=float((use&rejected).sum()/use.sum()),
                    mean_GeV=float(e[keep].mean()), rms_mean=float(rms), w68_median=float(w68),
                    gaussian_core=fit['sigma_mean'], gaussian_error=fit['error'],
                    fit_lo_GeV=fit['lo'], fit_hi_GeV=fit['hi'])
                metrics.append(row)
                path = OUT/'data'/condition/f'{cog}mm'/f'{beam}GeV'
                helper.write_csv(path/'metrics.csv', [row])
                np.savez_compressed(path/'energy_histogram.npz', edges_GeV=edges,
                    counts=np.histogram(e[keep], edges)[0],
                    underflow=int((e[keep] < edges[0]).sum()),
                    overflow=int((e[keep] > edges[-1]).sum()))
            for label, region in [('left', low), ('main', main)]:
                mask = use&region
                rates.append(dict(energy_GeV=beam, cog_mm=cog, group=label, n=int(mask.sum()),
                    fraction=float(mask.sum()/use.sum()),
                    early_fraction=float(np.mean(a['shower_start'][mask] <= 6)),
                    outside20_fraction=float(np.mean(a['cog_extent_mm'][mask] >= 20)),
                    energy_per_hit_limit_MeV=float(.8*1000*eref/nref)))
                for key in ('energy_GeV', 'hits', 'cog_extent_mm', 'shower_start', 'depth',
                            'back_fraction', 'front_fraction', 'width_mm', 'edge_fraction',
                            'zero_hits', 'unique_strips', 'memory_fraction', 'low_hit_fraction'):
                    q = np.quantile(a[key][mask], [.1, .5, .9])
                    features.append(dict(energy_GeV=beam, cog_mm=cog, group=label,
                        feature=key, q10=q[0], median=q[1], q90=q[2]))
        # Position comparison at ordinary Nhit, independent of energy selection.
        ordinary = base&(n > .7*nref)&(n < 1.35*nref)
        x, y = a['cog_x_mm'][:, 9], a['cog_y_mm'][:, 10]
        for axis, coords in [('x9', x), ('y10', y)]:
            for lo, hi in zip(range(-40, 40, 5), range(-35, 45, 5)):
                mask = ordinary&(coords >= lo)&(coords < hi)
                if mask.sum() < 20:
                    continue
                position.append(dict(energy_GeV=beam, axis=axis, lo_mm=lo, hi_mm=hi,
                    n=int(mask.sum()), median_energy_GeV=np.median(e[mask]),
                    median_hits=np.median(n[mask]), median_MeV_per_hit=1000*np.median(ratio[mask]),
                    median_width_mm=np.median(a['width_mm'][mask]),
                    median_start=np.median(a['shower_start'][mask])))
        fig, axs = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
        for ax, cog in zip(axs[0, :2], (200, 20)):
            use = base&(a['cog_extent_mm'] < cog)
            ax.hist2d(e[use], n[use], bins=100, range=((0, 1.5*eref), (0, 1.7*nref)),
                      norm=matplotlib.colors.LogNorm(), cmap='viridis')
            grid = np.linspace(.7*nref, 1.35*nref, 100)
            ax.plot(.8*eref/nref*grid, grid, 'r--', lw=1, label='left-band boundary')
            ax.set(xlabel='Calibrated signal sum [GeV]', ylabel='Nhit',
                   title=f'{beam} GeV, '+('nocut' if cog == 200 else 'CoG 20 mm'))
            ax.legend(fontsize=8)
        ax = axs[0, 2]
        bins = np.linspace(-50, 50, 26)
        count = np.histogram2d(x[ordinary], y[ordinary], bins=(bins, bins))[0]
        sums = np.histogram2d(x[ordinary], y[ordinary], bins=(bins, bins), weights=1000*ratio[ordinary])[0]
        mean = np.divide(sums, count, out=np.full(count.shape, np.nan), where=count >= 10)
        mesh = ax.pcolormesh(bins, bins, mean.T, shading='auto', cmap='coolwarm',
                             vmin=.65*1000*eref/nref, vmax=1.15*1000*eref/nref)
        fig.colorbar(mesh, ax=ax, label='Mean signal / hit [MeV]')
        ax.set(xlabel='Layer 9 CoG X [mm]', ylabel='Layer 10 CoG Y [mm]', title='Ordinary Nhit; no E selection')
        ax = axs[1, 0]
        for label, region, color in [('left', low, 'tab:red'), ('main', main, 'tab:blue')]:
            mask = base&region
            le = a['layer_energy_MeV'][mask]
            profile = (le/le.sum(axis=1)[:, None]).mean(axis=0)
            ax.plot(np.arange(30), profile, label=label, color=color)
        ax.set(xlabel='Layer (0 based)', ylabel='Mean fraction of signal', title='Longitudinal shape; nocut')
        ax.legend()
        ax = axs[1, 1]
        bins = np.linspace(7, 30, 70)
        for label, region, color in [('left', low, 'tab:red'), ('main', main, 'tab:blue')]:
            ax.hist(a['width_mm'][base&region], bins=bins, density=True, histtype='step', label=label, color=color)
        ax.set(xlabel='Energy-weighted within-layer transverse RMS [mm]', ylabel='Density')
        ax.legend()
        ax = axs[1, 2]
        bins = np.linspace(.4*eref, 1.3*eref, 90)
        for label, rejected, color in [('before', np.zeros(len(e), dtype=bool), 'black'),
                                       ('remove left band', low, 'tab:red'),
                                       ('width <= 14 mm', cuts['width14'], 'tab:green')]:
            ax.hist(e[central&~rejected], bins=bins, histtype='step', label=label, color=color)
        ax.set(xlabel='Calibrated signal sum [GeV]', ylabel='Events', title='CoG 20 mm')
        ax.legend(fontsize=8)
        fig.savefig(OUT/'data/diagnostics'/f'{beam}GeV/populations.png', dpi=150)
        plt.close(fig)
        print(f'{beam} GeV: summary completed', flush=True)
    helper.write_csv(OUT/'comparison/metrics.csv', metrics)
    helper.write_csv(OUT/'comparison/features.csv', features)
    helper.write_csv(OUT/'comparison/rates.csv', rates)
    helper.write_csv(OUT/'comparison/position.csv', position)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    for energy in (20, 30, 40):
        location = OUT/'data/diagnostics'/f'{energy}GeV/events.npz'
        if location.exists():
            continue
        collect(energy)
    if '--summarize' in sys.argv:
        with tempfile.TemporaryDirectory(prefix='nhit_plot_') as temporary:
            import matplotlib
            os.environ['MPLCONFIGDIR'] = temporary
            fontconfig = Path(temporary)/'fonts.conf'
            fontconfig.write_text('<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd"><fontconfig>'
                f'<dir>{matplotlib.get_data_path()}/fonts/ttf</dir><cachedir>{temporary}/fontcache</cachedir></fontconfig>')
            os.environ['FONTCONFIG_FILE'] = str(fontconfig)
            summarize()
