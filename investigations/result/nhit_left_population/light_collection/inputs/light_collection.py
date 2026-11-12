#!/usr/bin/env python3
"""Read-only summary of existing independent muon strip-position aggregates.

These are mean ADC ratios of acquired hits, not fitted MIP MPVs or optical
attenuation measurements. No correction is applied to electron events.
"""
import csv
import hashlib
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parents[2]
SOURCE = PROJECT / 'analysis/result/channel_response/strip_position_response/muon_allruns_100GeV'
OUT = PROJECT / 'investigations/result/nhit_left_population/light_collection'


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fingerprint(path):
    stat = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                sha256=digest.hexdigest())


def audit_inputs():
    counts, differences = {}, []
    data_base = PROJECT.parent / 'data_beamtest'
    for line in (SOURCE / 'inputs/input_stats_before.tsv').read_text().splitlines():
        original, size, mtime = line.split('\t')
        suffix = original.split('/data_beamtest/', 1)[1]
        current = data_base / suffix
        sample = 'data' if suffix.startswith('ECAL_data/') else 'mc'
        kind = 'raw' if '/decode/' in suffix else 'calib' if '/calib/' in suffix else 'truth'
        stat = current.stat() if current.exists() else None
        match = stat is not None and stat.st_size == int(size) and int(stat.st_mtime) == int(mtime)
        row = counts.setdefault(f'{sample}_{kind}', dict(all=0, matching_size_mtime=0, missing=0))
        row['all'] += 1
        row['matching_size_mtime'] += int(match)
        row['missing'] += int(stat is None)
        if not match:
            differences.append(dict(path=str(current), original_size=int(size), original_mtime=int(mtime),
                current_size=stat.st_size if stat else None, current_mtime=int(stat.st_mtime) if stat else None))
    return dict(counts=counts, differences=differences)


def parity_position():
    """Electron diagnostic from saved moments; shower CoG is a position proxy."""
    rows = []
    for beam in (20, 30, 40):
        path = OUT.parent / 'data/diagnostics' / f'{beam}GeV/events.npz'
        with np.load(path) as saved:
            a = {key: saved[key] for key in saved.files}
        base = (~a['tail_rich'].astype(bool)) & (a['cog_extent_mm'] < 200)
        nref = float(np.median(a['hits'][base & (a['cog_extent_mm'] < 20)]))
        ordinary = base & (a['hits'] > .7 * nref) & (a['hits'] < 1.35 * nref)
        x, y = a['cog_x_mm'][:, 9], a['cog_y_mm'][:, 10]
        for axis, coordinate, other, long_parity in [('X9', x, y, 0), ('Y10', y, x, 1)]:
            select = ordinary & (np.abs(other) < 5)
            for parity in (0, 1):
                energy = a['layer_energy_MeV'][:, parity::2].sum(axis=1)
                hits = a['layer_hits'][:, parity::2].sum(axis=1)
                response = np.divide(energy, hits, out=np.full(len(hits), np.nan), where=hits > 0)
                central = select & (coordinate >= -5) & (coordinate < 5) & np.isfinite(response)
                ref = float(np.median(response[central]))
                for lo, hi in [(-30, -20), (-5, 5), (20, 30)]:
                    use = select & (coordinate >= lo) & (coordinate < hi) & np.isfinite(response)
                    count = int(use.sum())
                    median = float(np.median(response[use])) if count else np.nan
                    rows.append(dict(beam_GeV=beam, axis=axis, other_abs_lt_mm=5,
                        lo_mm=lo, hi_mm=hi, parity='odd' if parity else 'even',
                        direction='long' if parity == long_parity else 'fine',
                        n=count, reference_n=int(central.sum()), normal_nhit_ref=nref,
                        median_MeV_per_hit=median, central_median_MeV_per_hit=ref,
                        ratio_to_central=median/ref,
                        median_parity_energy_MeV=float(np.median(energy[use])) if count else np.nan,
                        median_parity_hits=float(np.median(hits[use])) if count else np.nan))
    write_csv(OUT / 'comparison/parity_position.csv', rows)


def main():
    parity_position()
    mask_source = PROJECT / 'analysis/resolution/SingleEnergyAnalysis.cc'
    text = mask_source.read_text().split('const unordered_set<int> bad_channels = {', 1)[1].split('};', 1)[0]
    bad = set(map(int, re.findall(r'\d+', re.sub(r'//[^\n]*', '', text))))
    source_stat = fingerprint(SOURCE / 'bins.tsv')
    groups = defaultdict(dict)
    with (SOURCE / 'bins.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['axis'] != 'u' or row['condition'] not in ('unadjusted', 'position_matched'):
                continue
            cell = int(row['cellid'])
            if cell in bad:
                continue
            key = row['sample'], row['condition'], row['stage'], int(row['split']), cell
            groups[key][int(row['bin'])] = row
    ratios, summaries = [], []
    for sample in ('data', 'mc'):
        for condition in ('unadjusted', 'position_matched'):
            for stage in ('pedestal_adc', 'path_adc', 'temperature_adc'):
                # A fixed common channel cohort across both identity-hash splits.
                eligible = []
                cells = {key[4] for key in groups if key[:3] == (sample, condition, stage)}
                for cell in cells:
                    if all(len(groups[(sample, condition, stage, split, cell)]) == 3 and
                           all(int(row['usable']) >= 100 for row in groups[(sample, condition, stage, split, cell)].values())
                           for split in (0, 1)):
                        eligible.append(cell)
                for split in (0, 1):
                    for cell in sorted(eligible):
                        bins = groups[(sample, condition, stage, split, cell)]
                        ref = bins[1]
                        for side in (0, 2):
                            row = bins[side]
                            ratio = float(row['mean_ADC']) / float(ref['mean_ADC'])
                            error = ratio * np.hypot(float(row['mean_ADC_SE']) / float(row['mean_ADC']),
                                                     float(ref['mean_ADC_SE']) / float(ref['mean_ADC']))
                            ratios.append(dict(sample=sample, condition=condition, stage=stage,
                                split=split, cellid=cell, layer=cell // 100000, bin=side,
                                u_low_mm=-22.5 if side == 0 else 7.5,
                                u_high_mm=-7.5 if side == 0 else 22.5,
                                usable=int(row['usable']), central_usable=int(ref['usable']),
                                mean_ADC=float(row['mean_ADC']), central_mean_ADC=float(ref['mean_ADC']),
                                mean_ratio=ratio, mean_ratio_SE=float(error),
                                acquisition=float(row['acquisition']), central_acquisition=float(ref['acquisition'])))
                    for side in (0, 2):
                        for parity in ('all', 'even', 'odd'):
                            selected = [r for r in ratios if (r['sample'], r['condition'], r['stage'], r['split'], r['bin']) ==
                                        (sample, condition, stage, split, side) and
                                        (parity == 'all' or r['layer'] % 2 == (parity == 'odd'))]
                            if not selected:
                                continue
                            values = np.array([r['mean_ratio'] for r in selected])
                            q = np.quantile(values, [.1, .5, .9])
                            summaries.append(dict(sample=sample, condition=condition, stage=stage,
                                split=split, bin=side, parity=parity, channels=len(selected),
                                ratio_q10=float(q[0]), ratio_median=float(q[1]), ratio_q90=float(q[2]),
                                fraction_ratio_below_0p9=float(np.mean(values < .9)),
                                median_acquisition=float(np.median([r['acquisition'] for r in selected])),
                                median_central_acquisition=float(np.median([r['central_acquisition'] for r in selected]))))
    assert source_stat == fingerprint(SOURCE / 'bins.tsv'), 'Muon aggregates changed during summary'
    write_csv(OUT / 'comparison/summary.csv', summaries)
    for sample in ('data', 'mc'):
        for condition in ('unadjusted', 'position_matched'):
            write_csv(OUT / sample / condition / 'channel_ratios.csv',
                      [r for r in ratios if r['sample'] == sample and r['condition'] == condition])
    # Same physical channels for the qualitative data/MC comparison. Each sample
    # is matched internally; data and MC have no common temperature support.
    paired_summary = []
    for condition in ('unadjusted', 'position_matched'):
        for stage in ('pedestal_adc', 'path_adc', 'temperature_adc'):
            selected = [r for r in ratios if r['condition'] == condition and r['stage'] == stage]
            cohort = set.intersection(*[{r['cellid'] for r in selected if r['sample'] == s} for s in ('data', 'mc')])
            for sample in ('data', 'mc'):
                for split in (0, 1):
                    for side in (0, 2):
                        values = np.array([r['mean_ratio'] for r in selected if r['sample'] == sample and
                            r['split'] == split and r['bin'] == side and r['cellid'] in cohort])
                        if len(values):
                            q = np.quantile(values, [.1, .5, .9])
                            paired_summary.append(dict(sample=sample, condition=condition, stage=stage,
                                split=split, bin=side, common_channels=len(values), ratio_q10=float(q[0]),
                                ratio_median=float(q[1]), ratio_q90=float(q[2])))
    write_csv(OUT / 'comparison/common_channel_summary.csv', paired_summary)
    inp = OUT / 'inputs'
    inp.mkdir(parents=True, exist_ok=True)
    input_audit = audit_inputs()
    (inp / 'input_stat_audit.json').write_text(json.dumps(input_audit, indent=2) + '\n')
    provenance = dict(command='python investigations/nhit_left_population/light_collection.py',
        source=source_stat, source_status='COLLECTED exists; COMPLETE absent; original MPV fits unfinished',
        settings=dict(axis='u, strip long axis', bins_mm=[-22.5, -7.5, 7.5, 22.5],
            min_usable_per_bin_per_split=100, cohort='common channels satisfying all 3 bins in both splits',
            mask='same original CellID mask as current resolution source',
            weighting='equal channel; each ratio relative to its own central bin',
            position_matching='existing within-sample run, temperature, slope, other-position matching',
            electron_parity=dict(selection='tail excluded; valid CoG<200mm; 0.7<Nhit/Nref<1.35; Nref=CoG20 median',
                position='X9 scan requires abs(Y10)<5mm, Y10 scan requires abs(X9)<5mm',
                statistic='median of event-level parity signal/parity hits; central [-5,5) reference',
                directions='X long=even, Y long=odd; CoG is an energy-dependent proxy')),
        limitations=['Mean ADC of acquired usable HG hits, not plateau-fitted MIP MPV or optical attenuation.',
            'Threshold and acquisition effects can bias observed mean ADC ratios.',
            'Track position resolution and bin migration not independently measured.',
            'No verified mapping of positive u to the SiPM end.',
            'Electron parity scan can mix channel, strip-segment-boundary, shower-shape and CoG-response effects.',
            'Muon sample is independent of electron sample; no electron response correction or causality claim.'],
        dependencies=[fingerprint(mask_source), fingerprint(Path(__file__)),
            fingerprint(SOURCE / 'inputs/configuration.tsv'), fingerprint(SOURCE / 'inputs/manifest.tsv'),
            *[fingerprint(path) for path in sorted((SOURCE / 'inputs/constants').glob('*.root'))],
            *[fingerprint(OUT.parent / 'data/diagnostics' / f'{beam}GeV/events.npz') for beam in (20, 30, 40)]])
    (inp / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    shutil.copyfile(Path(__file__), inp / 'light_collection.py')
    (OUT / 'run.log').write_text(provenance['command'] + '\nExisting aggregate input unchanged; summary completed.\n')
    for row in summaries:
        if row['condition'] == 'position_matched' and row['stage'] == 'temperature_adc' and row['parity'] == 'all':
            print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
