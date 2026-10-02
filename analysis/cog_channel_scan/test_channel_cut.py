#!/usr/bin/env python3
"""Run the real ROOT executable and verify independent CoG/channel cuts (requires PyROOT)."""
from array import array
from collections import Counter
from pathlib import Path
import math
import random
import re
import subprocess

import ROOT

ROOT.gROOT.SetBatch(True)
base = Path(__file__).resolve().parent
out = base.parent / 'result' / 'cog_channel_scan' / 'test_channel_cut'
out.mkdir(parents=True, exist_ok=True)
source = (base / 'SingleEnergyAnalysis.cc').read_text()
bad_block = source.split('const unordered_set<int> bad_channels = {', 1)[1].split('};', 1)[0]
bad = set(map(int, re.findall(r'\b\d+\b', bad_block)))
good = {
    layer: [layer * 100000 + chip * 10000 + channel
            for chip in range(6) for channel in range(36)
            if layer * 100000 + chip * 10000 + channel not in bad]
    for layer in range(30)
}
rng = random.Random(20260929)
# Channel acceptance and CoG outcomes are specified independently of the C++ logic.
# Each hit: cell ID, energy, x, y, accepted at 5 / 200 / 22.5 mm.
events = [[], [(good[9][0], 0.75, 40., 40., False, True, False)]]
cog_pass = [(False, False), (False, False)]  # CoG half-widths 2.5 / 5 mm.
bad9 = next(cid for cid in sorted(bad) if cid // 100000 == 9)
for event_id in range(600):
    hits = []
    scale = rng.uniform(0.8, 1.2)
    for layer in range(30):
        if layer in (9, 10):
            continue
        for idx in range(rng.randint(1, 5)):
            x, y = ((4.999, -4.999) if idx == 0 else (-3., 4.))
            hits.append((good[layer][idx], rng.uniform(0.5, 1.5) * scale, x, y, True, True, True))
        # Both signs and both coordinates, including the strict boundary.
        for idx, (x, y) in enumerate(((5., 0.), (-5., 0.), (0., 5.), (0., -5.), (7., 0.), (0., -7.)), 5):
            hits.append((good[layer][idx], 0.125, x, y, False, True, True))
        hits.append((good[layer][11], -2., 0., 0., False, False, False))
    hits.extend(((1, 9., 0., 0., False, False, False),
                 (3000000, 8., 0., 0., False, False, False),
                 (3100000, 8., 0., 0., False, False, False),
                 (good[0][12], 2., 22.5, 0., False, True, False),
                 (good[0][13], 2., -22.499, 0., False, True, True)))
    center9 = (good[9][0], 1., 0., 0., True, True, True)
    center10 = (good[10][0], 1., 0., 0., True, True, True)
    # Outside-channel energy must affect CoG; missing/zero-energy layers fail;
    # the cut is energy-weighted, rectangular, strict, and applied to BOTH layers.
    scenarios = [
        ([center9, center10], (True, True)),
        ([center9, (good[9][1], 1., 40., 0., False, True, False), center10], (False, False)),
        ([center9, center10, (good[10][1], 1., 0., -40., False, True, False)], (False, False)),
        ([center9, (good[9][1], 1., -40., 0., False, True, False),
          (good[9][2], 1., 40., 0., False, True, False), center10], (True, True)),
        ([(good[9][0], 1., 2.5, 0., True, True, True), center10], (False, True)),
        ([center9, (good[10][0], 1., 0., 5., False, True, True)], (False, False)),
        ([center10], (False, False)),
        ([(good[9][0], 0., 0., 0., True, True, True), center10], (False, False)),
        ([center9, (good[10][0], 1., 0., -2.5, True, True, True)], (False, True)),
        ([(good[9][0], 1., -5., 0., False, True, True), center10], (False, False)),
        ([center9, center10, (bad9, 100., 40., 0., False, False, False),
          (good[9][1], -100., 40., 0., False, False, False)], (True, True)),
        ([(good[9][0], 3., 0., 0., True, True, True),
          (good[9][1], 1., 8., 0., False, True, True), center10], (True, True)),
        ([center9, (good[10][0], 1., 0., 2.499, True, True, True)], (True, True)),
        ([(good[9][0], 1., 4., 0., True, True, True), center10], (False, True)),
        ([(good[9][0], 1., 2., 2., True, True, True),
          (good[10][0], 1., -2., -2., True, True, True)], (True, True)),
    ]
    cog_hits, expected_cog = scenarios[event_id % len(scenarios)]
    hits.extend(cog_hits)
    events.append(hits)
    cog_pass.append(expected_cog)

fixture = out / 'calibration.root'
f = ROOT.TFile(str(fixture), 'RECREATE')
tree = ROOT.TTree('Calib_Hit', 'Synthetic channel-cut regression input')
cell = ROOT.std.vector('int')()
values = {name: ROOT.std.vector('double')() for name in
          ('Hit_Energy', 'Hit_X', 'Hit_Y', 'Hit_Z', 'NewTemperature')}
total = array('d', [0.])
tree.Branch('CellID', cell)
for name, value in values.items():
    tree.Branch(name, value)
tree.Branch('TotalEnergyDep', total, 'TotalEnergyDep/D')
for hits in events:
    cell.clear()
    for value in values.values():
        value.clear()
    total[0] = sum(hit[1] for hit in hits)
    for cid, energy, x, y, *_ in hits:
        cell.push_back(cid)
        for name, value in zip(values, (energy, x, y, 0., 23.)):
            values[name].push_back(value)
    tree.Fill()
tree.Write()
f.Close()


def check_hist(hist, values):
    assert hist, 'Missing histogram'
    counts = Counter(hist.FindBin(value) for value in values)
    assert int(hist.GetEntries()) == len(values), (hist.GetName(), hist.GetEntries(), len(values))
    for bin_id in range(hist.GetNbinsX() + 2):
        assert hist.GetBinContent(bin_id) == counts[bin_id], (hist.GetName(), bin_id)


for cog_range, aperture, flag, cog_flag in ((None, 5, 4, None), (0, 200, 5, None),
                                             (2.5, 22.5, 6, 0), (5, 22.5, 6, 1)):
    condition = f'cog{cog_range or 0}mm_channel{aperture}mm'
    figures = out / condition / 'figures'
    (figures / '1GeV' / 'raw').mkdir(parents=True, exist_ok=True)
    result = out / f'1GeV_{condition}.root'
    print(f'Running {condition} integration test...', flush=True)
    with (out / f'{condition}.log').open('w') as log:
        subprocess.run([str(base / 'SingleEnergyAnalysis'), str(result), '1', '1', '1',
                        str(fixture), *([] if cog_range is None else [str(cog_range)]),
                        str(aperture), '1', str(figures)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    result_file = ROOT.TFile.Open(str(result))
    assert result_file and not result_file.IsZombie()
    assert result_file.Get('cog_range_mm').GetVal() == (cog_range or 0)
    assert result_file.Get('channel_range_mm').GetVal() == aperture
    energy_fit = result_file.Get('beforeEventCut/fit_gaus')
    hit_fit = result_file.Get('beforeEventCut/fit_gaus_nhit')
    assert energy_fit and hit_fit
    assert math.isfinite(energy_fit.GetParameter(1)) and energy_fit.GetParameter(2) > 0
    assert math.isfinite(hit_fit.GetParameter(1)) and hit_fit.GetParameter(2) > 0
    selected = [[hit for hit in hits if hit[flag]] for hits in events]
    sums = [sum(hit[1] for hit in hits) for hits in selected]
    counts = list(map(len, selected))
    for directory in ('beforeEventCut', 'afterEventCut'):
        indices = [i for i in range(len(events))
                   if (cog_flag is None or cog_pass[i][cog_flag])
                   and (directory == 'beforeEventCut' or
                   ((sums[i] - energy_fit.GetParameter(1)) / (3 * energy_fit.GetParameter(2))) ** 2 +
                   ((counts[i] - hit_fit.GetParameter(1)) / (3 * hit_fit.GetParameter(2))) ** 2 <= 1)]
        check_hist(result_file.Get(f'{directory}/edep_1GeV'), [sums[i] for i in indices])
        check_hist(result_file.Get(f'{directory}/nhits_1GeV'), [counts[i] for i in indices])
        joint = result_file.Get(f'{directory}/hit_vs_e_1GeV')
        expected_joint = Counter(joint.FindBin(sums[i], counts[i]) for i in indices)
        assert joint.GetEntries() == len(indices)
        for bin_id, count in expected_joint.items():
            assert joint.GetBinContent(bin_id) == count, (directory, bin_id)
        for layer in range(30):
            layer_hits = [[hit for hit in selected[i] if hit[0] // 100000 == layer] for i in indices]
            check_hist(result_file.Get(f'{directory}/edep_1GeV_layer{layer}'),
                       [sum(hit[1] for hit in hits) for hits in layer_hits])
            check_hist(result_file.Get(f'{directory}/nhits_1GeV_layer{layer}'), list(map(len, layer_hits)))
        for sipm, layers in (('10um', range(4, 28)), ('15um', (0, 1, 2, 3, 28, 29))):
            check_hist(result_file.Get(f'{directory}/edep_1hit_{sipm}'),
                       [hit[1] for i in indices for hit in selected[i] if hit[0] // 100000 in layers])
        print(f'PASS {condition} {directory}: {len(indices)} events; energy, hit counts, joint E/Nhit, 30 layers and SiPM spectra', flush=True)
    result_file.Close()
print(f'All CoG/channel-cut integration checks passed. Output: {out}', flush=True)
