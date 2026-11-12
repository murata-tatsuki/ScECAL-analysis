#!/usr/bin/env python3
"""Compare shower shape and instrumented-edge signals on the saved event sample."""
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

sys.dont_write_bytecode = True
PROJECT = Path(__file__).resolve().parents[2]
BASE = PROJECT/'investigations/result/nhit_left_population'
OUT = BASE/'shower_leakage'
MASK = BASE/'inputs/resolution_source.cc'
text = MASK.read_text().split('const unordered_set<int> bad_channels = {', 1)[1].split('};', 1)[0]
BAD = np.array(sorted(set(map(int, re.findall(r'\d+', re.sub(r'//[^\n]*', '', text))))))


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def groups(a):
    e, n = a['energy_GeV'], a['hits']
    base = (~a['tail_rich'])&(a['cog_extent_mm'] < 200)
    central = base&(a['cog_extent_mm'] < 20)
    er, nr = np.median(e[central]), np.median(n[central])
    ordinary = base&(n > .7*nr)&(n < 1.35*nr)
    ratio = e/np.maximum(n, 1)/(er/nr)
    return ordinary&(ratio < .8), ordinary&(ratio > .9)&(ratio < 1.1), nr


def weighted_median(v, w):
    order = np.argsort(v)
    return np.interp(.5*w.sum(), np.cumsum(w[order])-.5*w[order], v[order])


def compare_cached():
    rows, matched = [], []
    for beam in (20, 30, 40):
        a = dict(np.load(BASE/'data/diagnostics'/f'{beam}GeV/events.npz'))
        left, main, nr = groups(a)
        le = a['layer_energy_MeV']
        features = dict(width_mm=a['width_mm'], depth_layer=a['depth'],
            back5_fraction=a['back_fraction'], back5_MeV=le[:, 25:].sum(axis=1),
            back3_fraction=le[:, 27:].sum(axis=1)/np.maximum(le.sum(axis=1), 1e-15),
            back3_MeV=le[:, 27:].sum(axis=1), first5_fraction=a['front_fraction'],
            first5_MeV=le[:, :5].sum(axis=1), first_layer_hits=a['layer_hits'][:, 0],
            first_layer_MeV=le[:, 0], late=a['shower_start'] > 6)
        for cog in (200, 20):
            for label, use in [('left', left), ('main', main)]:
                use = use&(a['cog_extent_mm'] < cog)
                for name, value in features.items():
                    value = np.asarray(value, dtype=float)
                    rows.append(dict(energy_GeV=beam, cog_mm=cog, group=label, n=int(use.sum()),
                        feature=name, mean=float(value[use].mean()), median=float(np.median(value[use])),
                        q90=float(np.quantile(value[use], .9))))
        # Overlap weighting in reconstructed position bins; no energy matching.
        x, y = a['cog_x_mm'][:, 9], a['cog_y_mm'][:, 10]
        for include_n in (False, True):
            strata = np.column_stack([np.floor(x/5), np.floor(y/5), a['run']])
            if include_n:
                strata = np.column_stack([strata, np.floor(a['hits']/(.1*nr))])
            both = np.flatnonzero(left|main)
            _, inverse = np.unique(strata[both], axis=0, return_inverse=True)
            nl = np.bincount(inverse, weights=left[both])
            nm = np.bincount(inverse, weights=main[both])
            target = np.minimum(nl, nm)
            for label, mask, counts in [('left', left, nl), ('main', main, nm)]:
                use = mask[both]&(target[inverse] > 0)
                ids = both[use]
                weights = target[inverse[use]]/counts[inverse[use]]
                for name, value in features.items():
                    v = np.asarray(value[ids], dtype=float)
                    matched.append(dict(energy_GeV=beam,
                        matching='xy_run_nhit' if include_n else 'xy_run', group=label,
                        n_raw=len(ids), n_weighted=float(weights.sum()),
                        n_effective=float(weights.sum()**2/np.dot(weights, weights)), feature=name,
                        mean=float(np.average(v, weights=weights)), median=float(weighted_median(v, weights))))
    write(OUT/'comparison/shape.csv', rows)
    write(OUT/'comparison/shape_position_matched.csv', matched)


def edge_subset(beam):
    a = dict(np.load(BASE/'data/diagnostics'/f'{beam}GeV/events.npz'))
    selected = np.arange(0, len(a['entry']), 3)
    records = json.loads((BASE/'inputs'/f'{beam}GeV/manifest.json').read_text())
    result = {k: a[k][selected].copy() for k in ('entry', 'file_index', 'run', 'event_time', 'event',
                                               'energy_GeV', 'tail_rich', 'cog_extent_mm')}
    for key in ('fine_edge_MeV', 'long_edge_MeV', 'outer_union_MeV', 'hit_width_mm'):
        result[key] = np.full(len(selected), np.nan)
    checked = 0
    for record in records:
        path = Path(record['path'])
        before = path.stat()
        assert before.st_size == record['bytes'] and before.st_mtime_ns == record['mtime_ns']
        with uproot.open(path, array_cache=None) as root:
            tree = root['Calib_Hit']
            assert tree.num_entries == record['entries']
            for lo, hi in record['spans']:
                rr = np.flatnonzero((result['file_index'] == record['file_index'])&
                    (result['entry'] >= lo)&(result['entry'] < hi))
                if not len(rr):
                    continue
                arrays = tree.arrays(['Run_Num', 'Event_Time', 'Event_Num', 'CellID',
                    'Hit_Energy', 'Hit_X', 'Hit_Y'], entry_start=lo, entry_stop=hi, library='ak')
                for row in rr:
                    j = int(result['entry'][row])-lo
                    assert tuple(int(arrays[k][j]) for k in ('Run_Num', 'Event_Time', 'Event_Num')) == \
                           tuple(int(result[k][row]) for k in ('run', 'event_time', 'event'))
                    ids, e, x, y = [ak.to_numpy(arrays[k][j]) for k in ('CellID', 'Hit_Energy', 'Hit_X', 'Hit_Y')]
                    layer = ids//100000
                    good = (layer < 30)&(e >= 0)&~np.isin(ids, BAD)
                    l, e = layer[good], e[good]
                    fine = np.where(l % 2, x[good], y[good])
                    long = np.where(l % 2, y[good], x[good])
                    assert abs(e.sum()/1000-result['energy_GeV'][row]) < 1e-12
                    fine_edge = np.abs(fine)+2.5 > 111.15-10
                    long_edge = np.abs(long)+22.5 > 113.3-10
                    result['fine_edge_MeV'][row] = e[fine_edge].sum()
                    result['long_edge_MeV'][row] = e[long_edge].sum()
                    result['outer_union_MeV'][row] = e[fine_edge|long_edge].sum()
                    counts = np.bincount(l, minlength=30)
                    sums = np.bincount(l, weights=fine, minlength=30)
                    sums2 = np.bincount(l, weights=fine*fine, minlength=30)
                    centered = sums2-np.divide(sums*sums, counts, out=np.zeros(30), where=counts > 0)
                    result['hit_width_mm'][row] = np.sqrt(max(0, centered.sum())/max(len(e), 1))
                    checked += 1
        after = path.stat()
        assert before.st_size == after.st_size and before.st_mtime_ns == after.st_mtime_ns
    assert checked == len(selected)
    dest = OUT/'data/diagnostics'/f'{beam}GeV'
    dest.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest/'edge_events.npz', source_row=selected, **result)
    left, main, nr = groups(a)
    geometry = []
    x, y = a['cog_x_mm'][selected, 9], a['cog_y_mm'][selected, 10]
    strata = np.column_stack([np.floor(x/5), np.floor(y/5), a['run'][selected],
                              np.floor(a['hits'][selected]/(.1*nr))])
    both = np.flatnonzero((left|main)[selected])
    _, inv = np.unique(strata[both], axis=0, return_inverse=True)
    nl = np.bincount(inv, weights=left[selected[both]])
    nm = np.bincount(inv, weights=main[selected[both]])
    target = np.minimum(nl, nm)
    for name, group, counts in [('left', left, nl), ('main', main, nm)]:
        mask = group[selected]
        value = result['hit_width_mm'][mask]
        geometry.append(dict(energy_GeV=beam, condition='unmatched', group=name,
            n_raw=len(value), n_effective=len(value), mean_mm=value.mean(), median_mm=np.median(value)))
        use = mask[both]&(target[inv] > 0)
        weights = target[inv[use]]/counts[inv[use]]
        value = result['hit_width_mm'][both[use]]
        geometry.append(dict(energy_GeV=beam, condition='xy_run_nhit_matched', group=name,
            n_raw=len(value), n_effective=weights.sum()**2/np.dot(weights, weights),
            mean_mm=np.average(value, weights=weights), median_mm=weighted_median(value, weights)))
    write(dest/'hit_geometry.csv', geometry)
    rows = []
    for cog in (200, 20):
        for name, group in [('left', left), ('main', main)]:
            use = group[selected]&(result['cog_extent_mm'] < cog)
            for feature in ('fine_edge_MeV', 'long_edge_MeV', 'outer_union_MeV'):
                value = result[feature][use]
                fraction = value/(result['energy_GeV'][use]*1000)
                rows.append(dict(energy_GeV=beam, cog_mm=cog, group=name, n=int(use.sum()),
                    feature=feature, median_MeV=float(np.median(value)), mean_MeV=float(value.mean()),
                    median_fraction=float(np.median(fraction)), mean_fraction=float(fraction.mean()),
                    fraction_above_3pct=float((fraction > .03).mean()),
                    fraction_above_5pct=float((fraction > .05).mean())))
    return rows


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    compare_cached()
    rows = []
    for beam in (20, 30, 40):
        rows.extend(edge_subset(beam))
        print(f'{beam} GeV: 6000 sampled events rechecked; edges completed', flush=True)
    write(OUT/'comparison/edge.csv', rows)
    (OUT/'run.json').write_text(json.dumps(dict(command='python -u investigations/nhit_left_population/shower_leakage.py',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        mask_sha256=hashlib.sha256(MASK.read_bytes()).hexdigest(),
        source_sample=str(BASE/'data/diagnostics'), energies_GeV=[20, 30, 40],
        selection='Same left/main definitions as check.py; HG-LG tail excluded; valid CoG<200 or20',
        matching='5mm x9,y10 bins and run, optionally 0.1*Nref Nhit bins; each group weighted to min(group counts) in common cells',
        outer='Entire strip counted if it intersects physical outer 10mm band: fine abs(q)+2.5>111.15-10, long abs(q)+22.5>113.3-10',
        limitation='Observed boundary-strip and final-layer signals are proxies, not escaped energy. Energy-derived position/width/fractions can correlate with response without proving causality.',
        hit_geometry='Hit-weighted within-layer fine-coordinate RMS uses no hit-energy weights; matching uses the same 5mm CoG/run/0.1*Nref bins on the 6000-event subset',
        validation='All reread IDs and energy sums match, all original input sizes/mtime unchanged'), indent=2)+'\n')
