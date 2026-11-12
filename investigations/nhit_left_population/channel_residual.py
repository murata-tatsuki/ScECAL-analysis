#!/usr/bin/env python3
"""Read-only 30 GeV check of channel and gain response in the left population."""
import csv
import hashlib
import json
import re
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

PROJECT = Path(__file__).resolve().parents[2]
BASE = PROJECT / 'investigations/result/nhit_left_population'
OUT = BASE / 'channel_residual'
PARAM = Path('/gpfs/group/ilc/users/murata/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share')


def physical(ids):
    return ids // 100000 * 100000 + ids % 100000 // 10000 * 10000 + ids % 100


def metadata(p):
    s = p.stat()
    return dict(path=str(p), bytes=s.st_size, mtime_ns=s.st_mtime_ns)


def csvwrite(p, rows):
    with p.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = dict(np.load(BASE / 'data/diagnostics/30GeV/events.npz'))
    manifest = json.loads((BASE / 'inputs/30GeV/manifest.json').read_text())
    baseline = (~a['tail_rich']) & (a['cog_extent_mm'] < 200)
    center = baseline & (a['cog_extent_mm'] < 20)
    er, nr = np.median(a['energy_GeV'][center]), np.median(a['hits'][center])
    ordinary = (a['hits'] > .7 * nr) & (a['hits'] < 1.35 * nr)
    ratio = a['energy_GeV'] / np.maximum(a['hits'], 1)
    left = ordinary & (ratio < .8 * er / nr)
    main_pop = ordinary & (ratio > .9 * er / nr) & (ratio < 1.1 * er / nr)
    source = PROJECT / 'analysis/resolution/SingleEnergyAnalysis.cc'
    s = source.read_text().split('const unordered_set<int> bad_channels = {', 1)[1].split('};', 1)[0]
    bad = np.array(sorted(set(map(int, re.findall(r'\d+', re.sub(r'//[^\n]*', '', s))))))
    cells = np.array([l * 100000 + chip * 10000 + c for l in range(30)
                      for chip in range(6) for c in range(30 if chip == 5 else 36)])
    idx = {int(c): i for i, c in enumerate(cells)}
    group = ((cells // 100000 > 3) & (cells // 100000 < 28)).astype(int)
    mip = np.full(len(cells), np.nan)
    rawmpv = np.full(len(cells), np.nan)
    pars = [PARAM / 'all_auto_muon_v4_trackfit.root', PARAM / 'all_hl_electron2023_hlratio_v3.root', PARAM / 'pedestal2023_SPS.root']
    with uproot.open(pars[0]) as f:
        t = f['MIP_Fit'].arrays(library='np')
    for cid, mpv, chi, ndf in zip(t['CellID'], t['LandauMPV'], t['ChiSquare'], t['NDF']):
        j = idx.get(int(physical(cid)))
        if j is None:
            continue
        rawmpv[j] = mpv
        if np.isfinite(mpv) and mpv > 0 and ndf > 0 and np.isfinite(chi) and chi >= 0 and chi / ndf <= (2 if group[j] else 1.7):
            mip[j] = mpv
    mip_fallback = ~np.isfinite(mip)
    for g in (0, 1):
        mip[(group == g) & mip_fallback] = np.mean(mip[(group == g) & ~mip_fallback], dtype=np.longdouble)
    hl_valid = np.zeros(len(cells), bool)
    gain, offset = np.full(len(cells), np.nan), np.full(len(cells), np.nan)
    with uproot.open(pars[1]) as f:
        t = f['InterCalib'].arrays(library='np')
    for cid, slope, intercept in zip(t['CellID'], t['Slope'], t['Intercept']):
        j = idx.get(int(physical(cid)))
        if j is not None and np.isfinite(slope) and .02 <= slope <= .05 and np.isfinite(intercept):
            hl_valid[j] = True
            gain[j], offset[j] = 1 / slope, -intercept / slope
    for g in (0, 1):
        for v in (gain, offset):
            v[(group == g) & ~hl_valid] = np.mean(v[(group == g) & hl_valid], dtype=np.longdouble)
    ph, pl = np.full(len(cells), np.nan), np.full(len(cells), np.nan)
    with uproot.open(pars[2]) as f:
        t = f['ChnLevel'].arrays(['CellID', 'PedHighMean', 'PedLowMean'])
    for cid, vh, vl in zip(ak.to_numpy(ak.flatten(t['CellID'])), ak.to_numpy(ak.flatten(t['PedHighMean'])), ak.to_numpy(ak.flatten(t['PedLowMean']))):
        j = idx.get(int(physical(cid)))
        if j is not None:
            ph[j], pl[j] = float(np.float32(vh)), float(np.float32(vl))
    channel_energy = np.zeros((2, len(cells)))
    channel_hits = np.zeros((2, len(cells)), dtype=np.int64)
    rows, provenance = [], []
    for rec in manifest:
        cal = Path(rec['path'])
        raw = Path(str(cal).replace('/calib/', '/decode/'))
        assert metadata(cal) == {k: rec[k] for k in ('path', 'bytes', 'mtime_ns')}
        before = [metadata(p) for p in (cal, raw)]
        provenance.extend(before)
        with uproot.open(cal, array_cache=None) as fc, uproot.open(raw, array_cache=None) as fr:
            assert fc['Calib_Hit'].num_entries == fr['Raw_Hit'].num_entries
            for lo, hi in rec['spans']:
                hi = min(hi, lo + 250)
                ch = fc['Calib_Hit'].arrays(['Run_Num', 'Event_Time', 'Event_Num', 'CellID',
                    'Hit_Energy', 'Hit_HG_Energy', 'Hit_LG_Energy', 'NewTemperature'], entry_start=lo, entry_stop=hi)
                rh = fr['Raw_Hit'].arrays(['Run_Num', 'Event_Time', 'TriggerID', 'CellID',
                    'HitTag', 'HG_Charge', 'LG_Charge'], entry_start=lo, entry_stop=hi)
                for j, entry in enumerate(range(lo, hi)):
                    k = np.flatnonzero((a['file_index'] == rec['file_index']) & (a['entry'] == entry))
                    assert len(k) == 1
                    k = k[0]
                    for ck, rk, akkey in [('Run_Num', 'Run_Num', 'run'), ('Event_Time', 'Event_Time', 'event_time'), ('Event_Num', 'TriggerID', 'event')]:
                        assert int(ch[ck][j]) == int(rh[rk][j]) == int(a[akkey][k])
                    ids, en, eh, el, temp = [ak.to_numpy(ch[key][j]) for key in ('CellID', 'Hit_Energy', 'Hit_HG_Energy', 'Hit_LG_Energy', 'NewTemperature')]
                    rid, tag, hg, lg = [ak.to_numpy(rh[key][j]) for key in ('CellID', 'HitTag', 'HG_Charge', 'LG_Charge')]
                    assert len(np.unique(rid)) == len(rid)
                    mapping = {int(cid): i for i, cid in enumerate(rid)}
                    order = np.array([mapping[int(cid)] for cid in ids], dtype=int)
                    assert np.all(tag[order] == 1)
                    hg, lg = hg[order], lg[order]
                    good0 = (ids // 100000 < 30) & ~np.isin(ids, bad)
                    good = good0 & (en >= 0)
                    assert good.sum() == a['hits'][k]
                    assert abs(en[good].sum() / 1000 - a['energy_GeV'][k]) < 1e-8
                    pindex = np.array([idx[int(cid)] for cid in physical(ids[good])], dtype=int)
                    e, h, l = en[good], eh[good], el[good]
                    if not len(e):
                        # Empty events are outside the valid CoG baseline.
                        assert not baseline[k]
                        continue
                    hraw, lraw = hg[good], lg[good]
                    overlap = (hraw >= 800) & (hraw <= 2200)
                    delta = l[overlap] - h[overlap]
                    lowgain = hraw >= 2600
                    assert np.max(np.abs(e - np.where(lowgain, l, h)), initial=0) < 1e-8
                    factor = 1 - (temp[good] - 20) * np.where(group[pindex], 1.6 / 135, 3.5 / 230)
                    hcalc = (hraw - ph[pindex]) * factor / mip[pindex] * .305
                    lcalc = ((lraw - pl[pindex]) * factor * gain[pindex] + offset[pindex]) / mip[pindex] * .305
                    closure_h = float(np.max(np.abs(hcalc - h), initial=0))
                    closure_l = float(np.max(np.abs(lcalc - l), initial=0))
                    assert max(closure_h, closure_l) < 1e-8
                    row = dict(cache_index=int(k), file_index=int(rec['file_index']), entry=entry,
                        run=int(a['run'][k]), base=bool(baseline[k]), left=bool(left[k]), main=bool(main_pop[k]),
                        cog_extent_mm=float(a['cog_extent_mm'][k]), energy_GeV=float(a['energy_GeV'][k]), hits=int(good.sum()),
                        mip_fallback_hit_fraction=float(mip_fallback[pindex].mean()),
                        mip_fallback_energy_fraction=float(e[mip_fallback[pindex]].sum() / e.sum()),
                        hl_fallback_energy_fraction=float(e[~hl_valid[pindex]].sum() / e.sum()),
                        lg_energy_fraction=float(e[lowgain].sum() / e.sum()),
                        overlap_hits=int(overlap.sum()), overlap_delta_sum_MeV=float(delta.sum()),
                        overlap_delta2_sum_MeV2=float(np.dot(delta, delta)),
                        overlap_delta_median_MeV=float(np.median(delta)) if len(delta) else np.nan,
                        negative_selected_MeV=float(-en[good0 & (en < 0)].sum()),
                        lg3000_hits=int((lraw >= 3000).sum()), temperature_C=float(np.mean(temp[good])),
                        closure_h_absmax_MeV=closure_h, closure_l_absmax_MeV=closure_l)
                    rows.append(row)
                    if baseline[k] and (left[k] or main_pop[k]):
                        pop = 0 if left[k] else 1
                        channel_energy[pop] += np.bincount(pindex, weights=e, minlength=len(cells))
                        channel_hits[pop] += np.bincount(pindex, minlength=len(cells))
                print(f'run {int(a["run"][k])} entry {lo}:{hi}, total {len(rows)}', flush=True)
        assert before == [metadata(p) for p in (cal, raw)]
    csvwrite(OUT / 'events.csv', rows)
    event = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    np.savez_compressed(OUT / 'events.npz', **event)
    stats = []
    for cog in (200, 20):
        for pop in ('left', 'main'):
            use = event['base'] & event[pop] & (event['cog_extent_mm'] < cog)
            n = int(use.sum())
            out = dict(energy_GeV=30, cog_mm=cog, population=pop, events=n)
            for key in ('mip_fallback_hit_fraction', 'mip_fallback_energy_fraction', 'hl_fallback_energy_fraction',
                        'lg_energy_fraction', 'negative_selected_MeV', 'temperature_C'):
                out[key + '_mean'] = float(np.mean(event[key][use]))
                out[key + '_median'] = float(np.median(event[key][use]))
            out['overlap_hits'] = int(event['overlap_hits'][use].sum())
            out['overlap_delta_mean_MeV'] = float(event['overlap_delta_sum_MeV'][use].sum() / out['overlap_hits'])
            out['overlap_delta_RMS_MeV'] = float(np.sqrt(event['overlap_delta2_sum_MeV2'][use].sum() / out['overlap_hits'] - out['overlap_delta_mean_MeV'] ** 2))
            out['lg3000_hits'] = int(event['lg3000_hits'][use].sum())
            stats.append(out)
    csvwrite(OUT / 'summary.csv', stats)
    csvwrite(OUT / 'channels.csv', [dict(cellid=int(cid), mip_adopted_ADC=mip[i], mip_raw_ADC=rawmpv[i], mip_fallback=int(mip_fallback[i]),
        hl_fallback=int(~hl_valid[i]), left_energy_MeV=channel_energy[0, i], main_energy_MeV=channel_energy[1, i],
        left_hits=int(channel_hits[0, i]), main_hits=int(channel_hits[1, i])) for i, cid in enumerate(cells)])
    (OUT / 'run.json').write_text(json.dumps(dict(command='python investigations/nhit_left_population/channel_residual.py',
        events=len(rows), sample='First 250 events of each of the original 24 stratified blocks; both 30 GeV runs',
        eref_GeV=er, nref=nr, mip_fallback=int(mip_fallback.sum()), hl_fallback=int((~hl_valid).sum()),
        closure_h_absmax_MeV=float(event['closure_h_absmax_MeV'].max()),
        closure_l_absmax_MeV=float(event['closure_l_absmax_MeV'].max()),
        inputs=provenance, parameters=[dict(**metadata(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in pars],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        validation='All event keys, full CellID mapping, HitTag, selected gain, cache energy/hits and current-parameter HG/LG formula match; input size/mtime unchanged'), indent=2) + '\n')
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == '__main__':
    main()
