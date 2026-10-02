#!/usr/bin/env python3
"""End-to-end event veto and selection provenance checks (requires PyROOT)."""
from array import array
from pathlib import Path
import argparse
import math
import random
import subprocess
import tempfile
import ROOT

ROOT.gROOT.SetBatch(True)
SOURCE = Path(__file__).resolve().parents[1]

def run(args, expected=None, log=None):
    if log:
        with log.open('w') as stream:
            result = subprocess.run(list(map(str, args)), text=True, stdout=stream, stderr=subprocess.STDOUT, timeout=600)
        result.stdout = log.read_text()
    else:
        result = subprocess.run(list(map(str, args)), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
    if expected:
        assert result.returncode != 0 and expected in result.stdout, result.stdout[-5000:]
    else:
        assert result.returncode == 0, result.stdout[-5000:]
    return result.stdout


def fixtures(work, omit=False, mismatch=False):
    calib = work / 'calib' / '10GeV' / 'events.root'
    calib.parent.mkdir(parents=True, exist_ok=True)
    output = ROOT.TFile(str(calib), 'RECREATE')
    tree = ROOT.TTree('Calib_Hit', '')
    scalars = {name: array('i', [0]) for name in ('Run_Num', 'Event_Time', 'Event_Num')}
    for name, value in scalars.items():
        tree.Branch(name, value, name + '/I')
    total = array('d', [0.])
    tree.Branch('TotalEnergyDep', total, 'TotalEnergyDep/D')
    ids = ROOT.std.vector('int')()
    tree.Branch('CellID', ids)
    vectors = {name: ROOT.std.vector('double')() for name in ('Hit_Energy', 'Hit_X', 'Hit_Y', 'Hit_Z', 'NewTemperature')}
    for name, value in vectors.items():
        tree.Branch(name, value)
    rng = random.Random(14)
    for entry in range(1200):
        scalars['Run_Num'][0] = 42
        scalars['Event_Time'][0] = entry + 100
        scalars['Event_Num'][0] = entry % 100  # Repeated trigger IDs must not collide.
        rich = entry % 5 == 0
        total[0] = rng.gauss(360 if rich else 420, 28 if rich else 14)
        ids.clear()
        for value in vectors.values():
            value.clear()
        for layer in range(30):
            # Vary hit counts so the existing nhit Gaussian fit is defined.
            for channel in range(2 + rng.randrange(4)):
                ids.push_back(layer * 100000 + channel)
                vectors['Hit_Energy'].push_back(total[0] / 100)
                vectors['Hit_X'].push_back(0.)
                vectors['Hit_Y'].push_back(0.)
                vectors['Hit_Z'].push_back(layer * 10.)
                vectors['NewTemperature'].push_back(20.)
        tree.Fill()
    tree.Write()
    output.Close()
    directory = work / 'tails' / '10GeV'
    (directory / 'inputs').mkdir(parents=True, exist_ok=True)
    (directory / 'inputs' / 'data.tsv').write_text(f'/unused/raw.root\t{calib}\n')
    (directory / 'inputs' / 'mc.tsv').write_text('')
    output = ROOT.TFile(str(directory / 'tail_study.root'), 'RECREATE')
    tree = ROOT.TTree('events', '')
    fields = {name: array('i', [0]) for name in ('sample', 'file_index', 'Run_Num', 'Event_Time', 'TriggerID', 'low_tail_channels')}
    for name, value in fields.items():
        tree.Branch(name, value, name + '/I')
    index, rich, selected = array('q', [0]), array('b', [0]), array('b', [0])
    tree.Branch('entry', index, 'entry/L')
    tree.Branch('tail_rich', rich, 'tail_rich/O')
    tree.Branch('selected', selected, 'selected/O')
    # Reverse order proves matching is by explicit entry, not tail-tree row order.
    for entry in reversed(range(1200)):
        if omit and entry == 777:
            continue
        index[0] = entry
        fields['Run_Num'][0] = 42
        fields['Event_Time'][0] = entry + 100
        fields['TriggerID'][0] = entry % 100 + (1 if mismatch and entry == 444 else 0)
        fields['low_tail_channels'][0] = 3 if entry % 5 == 0 else entry % 3
        rich[0] = entry % 5 == 0
        selected[0] = 0  # Must still veto tail-rich entries outside the study's selection.
        tree.Fill()
    tree.Write()
    ROOT.TParameter('int')('min_tail_channels', 3).Write()
    ROOT.TParameter('double')('nsigma', 5.).Write()
    ROOT.TNamed('tail_definition', 'synthetic low tail: z<-5; >=3 distinct channels').Write()
    output.Close()
    return calib, work / 'tails'


def checks(work):
    calib, tails = fixtures(work)
    single = SOURCE / 'SingleEnergyAnalysis'
    def command(mode):
        out = work / mode / '10GeV_20mm.root'
        figures = work / mode / 'figures'
        (figures / '10GeV' / 'raw').mkdir(parents=True, exist_ok=True)
        return [single, out, 1, '10GeV', 1, calib, 20, 0, figures, '--tail-results', tails]
    fixtures(work, omit=True)
    run(command('missing'), 'incomplete tail coverage')
    assert not (work / 'missing/10GeV_20mm.root').exists()
    fixtures(work, mismatch=True)
    run(command('mismatch'), 'event identity mismatch')
    fixtures(work)
    unknown = work / 'unknown.root'
    unknown.write_bytes(calib.read_bytes())
    args = command('unknown'); args[5] = unknown
    run(args, 'input is absent from tail manifest')
    print('PASS: incomplete coverage, mismatched event identity and unregistered input rejected before output', flush=True)
    hist_entries = {}
    fits = {}
    for mode, option in [('exclude', []), ('keep', ['--keep-tail-events'])]:
        run(command(mode) + option, log=work / f'{mode}.log')
        output = ROOT.TFile(str(work / mode / '10GeV_20mm.root'))
        hist_entries[mode] = output.Get('beforeEventCut/edep_10GeV').GetEntries()
        fit = output.Get('beforeEventCut/fit_gaus')
        fits[mode] = fit.GetParameter(2) / fit.GetParameter(1)
        assert math.isfinite(fits[mode]) and fits[mode] > 0
        assert output.Get('hl_tail_excluded').GetVal() == (mode == 'exclude')
        assert output.Get('hl_tail_status').GetTitle() == 'complete'
        if mode == 'exclude':
            assert output.Get('hl_tail_checked_events').GetVal() == 1200
            assert output.Get('hl_tail_rejected_events').GetVal() == 240
        assert output.Get('hl_tail_events_after_cog').GetVal() == hist_entries[mode]
        output.Close()
    assert hist_entries == {'exclude': 960, 'keep': 1200}, hist_entries
    print('PASS: default excludes exactly 240/1200 events, keep retains all, repeated triggers and selected=0 handled', flush=True)
    multi = SOURCE / 'MultiEnergyAnalysis'
    excluded = work / 'exclude/10GeV_20mm.root'
    kept = work / 'keep/10GeV_20mm.root'
    def comparison(first, second, mode):
        return [multi, work / f'comparison_{mode}.root', 2, 1, first, second, work / f'multi_{mode}', '--labels', 'A', 'B', '--skip-channel-plots']
    run(comparison(kept, kept, 'wrong'), 'tail selection mismatch')
    run(comparison(excluded, kept, 'mixed'), 'tail selection mismatch')
    run(comparison(excluded, excluded, 'wrongkeep') + ['--keep-tail-events'], 'tail selection mismatch')
    for mode, path, options in [('exclude', excluded, []), ('keep', kept, ['--keep-tail-events'])]:
        run(comparison(path, path, mode) + options, log=work / f'multi_{mode}.log')
        output = ROOT.TFile(str(work / f'comparison_{mode}.root'))
        canvas = output.Get('resolution')
        graphs = [o for o in canvas.GetListOfPrimitives() if o.InheritsFrom('TGraphErrors')]
        assert len(graphs) == 2
        assert all(abs(g.GetPointY(0) - fits[mode]) < 1e-12 for g in graphs)
        output.Close()
    print('PASS: Multi rejects mixed/wrong selections and preserves SingleEnergy Gaussian sigma/mu', flush=True)
    print('Fit resolutions:', fits, flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--work-dir', type=Path)
    options = parser.parse_args()
    if options.work_dir:
        options.work_dir.mkdir(parents=True, exist_ok=True)
        checks(options.work_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='resolution-hl-tail-') as directory:
            checks(Path(directory))
