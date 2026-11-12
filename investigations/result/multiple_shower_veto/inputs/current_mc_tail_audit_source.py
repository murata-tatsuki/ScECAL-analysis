#!/usr/bin/env python3
"""Apply the saved HG–LG residual model to the current sampled MC decode."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

PROJECT = Path(__file__).resolve().parents[2]
DATA = PROJECT.parent / 'data_beamtest'
RESULT = PROJECT / 'investigations/result/multiple_shower_veto'
OLD = '/net/megsrv01/export/megraid01/users/data_beamtest/'
PED_SHA = '9310639a8738c1e3bba07288bf4cffea07df2826e1f042049e5a2ff12adc5b4c'


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def relocated(text):
    require(text.startswith(OLD), 'Unexpected manifest path: ' + text)
    return (DATA / text[len(OLD):]).resolve(strict=True)


def metadata(path):
    stat = path.stat()
    return dict(path=str(path.resolve()), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def physical(ids):
    return ids // 100000 * 100000 + ids % 100000 // 10000 * 10000 + ids % 100


def sample_blocks(count, take, rng):
    # Identical to check.py; the completed manifest is asserted separately.
    if take == count:
        return [(start, min(start+1000, count)) for start in range(0, count, 1000)]
    blocks = min(12, max(1, take // 40))
    spans = []
    for i in range(blocks):
        lo, hi = count*i//blocks, count*(i+1)//blocks
        size = take//blocks + (i < take % blocks)
        require(hi-lo >= size, 'Invalid sample allocation')
        start = int(rng.integers(lo, hi-size+1))
        spans.append((start, start+size))
    return spans


def verify(beam, required=False):
    directory = RESULT / 'inputs' / f'{beam}GeV'
    path = directory / 'mc/manifest.json'
    audit = directory / 'current_mc_tail_audit.json'
    state = json.loads(audit.read_text())
    if not path.exists():
        require(not required, 'Main MC manifest not yet available: ' + str(path))
        return False
    require(json.loads(path.read_text()) == state['sampling_records'],
            'Main collection and tail audit sampling manifests differ')
    state['main_manifest_verified'] = True
    state['main_manifest_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    audit.write_text(json.dumps(state, indent=2)+'\n')
    print(f'{beam} GeV: main MC sampling manifest matches', flush=True)
    return True


def audit(beam, target):
    source = PROJECT / f'analysis/result/channel_response/hg_lg_tail/threshold/e_ssa/{beam}GeV'
    out = RESULT / 'inputs' / f'{beam}GeV'
    out.mkdir(parents=True, exist_ok=True)
    model_path = source / 'channels.tsv'
    model = list(csv.DictReader(model_path.open(), delimiter='\t'))
    common = {int(r['cell']): r for r in model if int(r['common']) == 1}
    require(common and all(int(r['good']) == 1 for r in common.values()), 'Invalid common mask')
    common_ids = np.array(sorted(common), dtype=np.int64)
    # Dense arrays avoid per-hit Python model lookups.
    size = max(int(r['cell']) for r in model) + 1
    gain, offset, sigma, median, ph, pl = [np.zeros(size) for _ in range(6)]
    for cell, row in common.items():
        gain[cell], offset[cell] = float(row['HL_A']), float(row['HL_B'])
        sigma[cell], median[cell] = float(row['noise_sigma_LG_ADC']), float(row['mc_median'])
    require(np.all(gain[common_ids] > 0) and np.all(sigma[common_ids] > 0), 'Invalid residual model')
    ped = DATA / 'simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/pedestal2023_SPS.root'
    ped_hash = hashlib.sha256(ped.read_bytes()).hexdigest()
    require(ped_hash == PED_SHA, 'Pedestal differs from saved model checksum')
    with uproot.open(ped) as root:
        arrays = root['ChnLevel'].arrays(['CellID', 'PedHighMean', 'PedLowMean'], library='ak')
    for ids, high, low in zip(arrays['CellID'], arrays['PedHighMean'], arrays['PedLowMean']):
        for cell, h, l in zip(physical(ak.to_numpy(ids)), ak.to_numpy(high), ak.to_numpy(low)):
            if int(cell) in common:
                ph[cell], pl[cell] = float(np.float32(h)), float(np.float32(l))
    pairs = [tuple(relocated(p) for p in line.split('\t')) for line in
             (source/'inputs/mc.tsv').read_text().splitlines() if line]
    count_rows = list(csv.DictReader((source/'inputs/sampled_files.tsv').open(), delimiter='\t'))
    counts = np.array([int(r['entries']) for r in count_rows if r['sample'] == 'mc'])
    require(len(pairs) == len(counts), 'Invalid source file counts')
    target = min(int(counts.sum()), target)
    desired = target * counts / counts.sum()
    takes = np.floor(desired).astype(int)
    for i in np.argsort(-(desired-takes))[:target-int(takes.sum())]:
        takes[i] += 1
    rng = np.random.default_rng(20261111 + beam*10 + 1)
    with uproot.open(source/'tail_study.root') as root:
        old = root['events'].arrays(['sample', 'file_index', 'entry', 'Run_Num', 'Event_Time',
                                     'TriggerID', 'tail_rich'], library='np')
    old = {k: a[old['sample'] == 1] for k, a in old.items()}
    records, input_records, rows = [], [], []
    with (out/'current_mc_tail_audit.log').open('w') as log:
        for index, ((raw_path, cal_path), count, take) in enumerate(zip(pairs, counts, takes)):
            spans = sample_blocks(int(count), int(take), rng)
            before_raw, before_cal = metadata(raw_path), metadata(cal_path)
            records.append(dict(file_index=index, entries=int(count), sampled=int(take), spans=spans, **before_cal))
            input_records.append(dict(file_index=index, raw=before_raw, calib=before_cal))
            oi = np.flatnonzero(old['file_index'] == index)
            oi = oi[np.argsort(old['entry'][oi])]
            require(np.array_equal(old['entry'][oi], np.arange(count)), 'Incomplete old event keys')
            with uproot.open(raw_path, array_cache=None) as raw, uproot.open(cal_path, array_cache=None) as cal:
                rtree, ctree = raw['Raw_Hit'], cal['Calib_Hit']
                require(rtree.num_entries == ctree.num_entries == count, 'Entry counts changed')
                for start, stop in spans:
                    r = rtree.arrays(['Run_Num', 'Event_Time', 'TriggerID', 'CellID', 'HitTag',
                                     'HG_Charge', 'LG_Charge'], entry_start=start, entry_stop=stop, library='ak')
                    c = ctree.arrays(['Run_Num', 'Event_Time', 'Event_Num', 'CellID'],
                                    entry_start=start, entry_stop=stop, library='ak')
                    for j, entry in enumerate(range(start, stop)):
                        keys = tuple(int(r[k][j]) for k in ('Run_Num', 'Event_Time', 'TriggerID'))
                        require(keys == tuple(int(c[k][j]) for k in ('Run_Num', 'Event_Time', 'Event_Num')),
                                'Current raw/calib event key mismatch')
                        k = oi[entry]
                        require(keys == tuple(int(old[f][k]) for f in ('Run_Num', 'Event_Time', 'TriggerID')),
                                'Saved event key mismatch')
                        ids, tags, h, l = (ak.to_numpy(r[f][j]) for f in
                                           ('CellID', 'HitTag', 'HG_Charge', 'LG_Charge'))
                        cell = physical(ids)
                        current_cells = np.unique(physical(ak.to_numpy(c['CellID'][j])))
                        use = ((tags == 1) & (h >= 800) & (h < 2200) &
                               np.isin(cell, common_ids) & np.isin(cell, current_cells))
                        cells, high, low = cell[use], h[use], l[use]
                        z = ((low-pl[cells]-(high-ph[cells])/gain[cells] +
                              offset[cells]/gain[cells]-median[cells])/sigma[cells])
                        low_cells = np.unique(cells[z < -5])
                        rows.append(dict(file_index=index, entry=entry, run=keys[0], event_time=keys[1],
                                         event=keys[2], old_tail_rich=int(old['tail_rich'][k]),
                                         tail_rich=int(len(low_cells) >= 3),
                                         current_tail_rich=int(len(low_cells) >= 3),
                                         low_tail_channels=len(low_cells), overlap_common_hits=int(use.sum()),
                                         min_z=float(z.min()) if len(z) else '',
                                         low_tail_cellids=';'.join(map(str, low_cells))))
            require(metadata(raw_path) == before_raw and metadata(cal_path) == before_cal, 'Input changed during audit')
            if (index+1) % 20 == 0:
                message = f'{beam} GeV: {index+1}/{len(pairs)} files, {len(rows)} sampled events'
                print(message, flush=True)
                print(message, file=log, flush=True)
        require(len(rows) == target, 'Incorrect sampled count')
        with (out/'current_mc_tail_audit.csv').open('w', newline='') as output:
            writer = csv.DictWriter(output, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        summary = dict(energy_GeV=beam, sampled=len(rows), current_tail_rich=sum(r['current_tail_rich'] for r in rows),
                       old_tail_rich=sum(r['old_tail_rich'] for r in rows),
                       events_with_any_low_tail=sum(r['low_tail_channels'] > 0 for r in rows),
                       maximum_low_tail_channels=max(r['low_tail_channels'] for r in rows))
        state = dict(**summary, model_condition='Saved common mask, MC median, sigma and HL factors fixed; current decode reclassified. Not a full refit of the current MC residual population.',
                     formula='z=((LG-pedLG)-(HG-pedHG)/A+B/A-mc_median)/noise_sigma_LG_ADC; low iff z<-5; event iff >=3 distinct physical cells',
                     hit_selection='HitTag=1; 800<=HG_raw<2200; saved common=1 (connected good channels); physical CellID present in current Calib_Hit',
                     pedestal_sha256=ped_hash, saved_pedestal_sha256=PED_SHA, pedestal_float32_cast=True,
                     model=dict(**metadata(model_path), sha256=hashlib.sha256(model_path.read_bytes()).hexdigest()),
                     old_tail_source=metadata(source/'tail_study.root'), sampling_records=records,
                     current_inputs=input_records, seed=20261111+beam*10+1, main_manifest_verified=False,
                     content_change_evidence=dict(method='Same saved event IDs, old tail event energy/good_hits compared with current calibration using saved physical good mask. Raw ADC checked against saved tail_hits.',
                       first_ten_energy_ratio_ranges={'5':[.960456123,1.039675499], '40':[.980501607,1.016010295], '100':[.988528836,1.002344633]},
                       raw_examples=[dict(energy=40,file_index=33,entry=910,cell=1200024,old_HG=1190,new_HG=1280,old_LG=455,new_LG=474),
                                     dict(energy=100,file_index=12,entry=766,cell=1200027,old_HG=994,new_HG=935,old_LG=447,new_LG=468)]),
                     source_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        (out/'current_mc_tail_audit.json').write_text(json.dumps(state, indent=2)+'\n')
        print(json.dumps(summary), file=log)
        print(json.dumps(summary), flush=True)
    verify(beam)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energy', type=int, nargs='+', default=[5, 40, 100])
    parser.add_argument('--target', type=int, default=20000)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    for energy in args.energy:
        verify(energy, required=True) if args.verify_only else audit(energy, args.target)
