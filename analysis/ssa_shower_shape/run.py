#!/usr/bin/env python3
"""Pair inputs, record provenance, aggregate and draw the SSA shape comparison."""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DATA = Path('/home/murata_t/data_beamtest')
MC = DATA / 'simulation/CEPCScECAL_SML_Portable_update_new/Result_MC'
ENERGIES = ['0.5', '1', '2', '3', '4', '5', '10', '20', '30', '40', '50', '60', '70', '80', '100', '120', '150', '200', '250']
SOURCES = ['Common.hh', 'Aggregate.cc', 'Plot.cc', 'Makefile', 'run.py', 'run.sh', 'submit.sh', 'README.md']

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def metadata(path):
    p = Path(path).resolve(strict=True)
    st = p.stat()
    return {'path': str(p), 'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def discover(pre, post, truth, max_files):
    files = sorted(pre.glob('*.root'))
    require(files, f'No calibration ROOT files: {pre}')
    missing = [p.name for p in files if not (post / p.name).is_file()]
    require(not missing, f'Missing matching SSA files in {post}: {missing[:5]}')
    if max_files and len(files) > max_files:
        # Equal file-index strata; select their midpoints (not just the first runs).
        files = [files[((2 * i + 1) * len(files)) // (2 * max_files)] for i in range(max_files)]
    rows = []
    for p in files:
        paths = [p.resolve(), (post / p.name).resolve(), (truth / p.name).resolve() if truth else None]
        require(all('\t' not in str(x) and '\n' not in str(x) for x in paths), 'Tabs/newlines in input paths are unsupported')
        if truth:
            require(paths[2].is_file(), f'Missing truth: {paths[2]}')
        rows.append(paths)
    return rows

def execute(command, log):
    print(shlex.join(map(str, command)), flush=True)
    with log.open('a') as stream:
        stream.write('\nCOMMAND: ' + shlex.join(map(str, command)) + '\n')
        stream.flush()
        subprocess.run(list(map(str, command)), check=True, stdout=stream, stderr=subprocess.STDOUT)

def snapshot(destination, files):
    destination.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for source in files:
        source = Path(source)
        shutil.copy2(source, destination / source.name)
        hashes[source.name] = sha(source)
    (destination / 'checksums.json').write_text(json.dumps(hashes, indent=2) + '\n')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energy', choices=ENERGIES, default='100')
    parser.add_argument('--label', default='baseline')
    parser.add_argument('--result-dir', type=Path, default=PROJECT / 'analysis/result/ssa_shower_shape')
    parser.add_argument('--max-events', type=int, default=6000, help='per sample, before axis selection; 0=all')
    parser.add_argument('--max-files', type=int, default=4, help='per sample, evenly selected; 0=all')
    parser.add_argument('--axis-layers', type=int, default=4)
    parser.add_argument('--seed', type=int, default=20261005)
    parser.add_argument('--truth', choices=['auto', 'required', 'off'], default='auto')
    parser.add_argument('--data-calib-dir', type=Path)
    parser.add_argument('--data-ssa-dir', type=Path)
    parser.add_argument('--mc-calib-dir', type=Path)
    parser.add_argument('--mc-ssa-dir', type=Path)
    parser.add_argument('--mc-truth-dir', type=Path)
    parser.add_argument('--beam-map', type=Path)
    parser.add_argument('--plots-only', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    require(re.fullmatch(r'[A-Za-z0-9_-]+', args.label), 'Use letters, digits, underscore or hyphen for label')
    require(args.max_events >= 0 and args.max_files >= 0 and 2 <= args.axis_layers <= 30 and 0 <= args.seed < 2**32 - 1, 'Invalid sampling/axis/seed option')
    energy_dir = args.energy + 'GeV'
    base = args.result_dir.resolve() / args.label
    record = base / 'inputs' / energy_dir
    logs = PROJECT / 'analysis/log/ssa_shower_shape' / args.label / energy_dir
    beam = 'ps' if float(args.energy) <= 5 else 'sps'
    db = DATA / 'ECAL_data/analysed/2023' / beam
    rows = {}
    truth_dir = None
    beam_map = None
    if not args.plots_only:
        pre_data = args.data_calib_dir or db / 'calib/e-' / energy_dir
        post_data = args.data_ssa_dir or db / 'ssaCalib/e-' / energy_dir
        pre_mc = args.mc_calib_dir or MC / 'calib/e-/e_ssa/threshold' / energy_dir
        post_mc = args.mc_ssa_dir or MC / 'ssaCalib/e-/e_ssa/threshold' / energy_dir
        candidate = args.mc_truth_dir or MC / 'generate/e-/e_ssa' / energy_dir
        if args.truth != 'off':
            if candidate.is_dir() and list(candidate.glob('*.root')):
                truth_dir = candidate
            else:
                require(args.truth == 'auto' and args.mc_truth_dir is None, f'Truth directory unavailable: {candidate}')
                print(f'Truth disabled: no files in {candidate}', flush=True)
        rows['data'] = discover(pre_data, post_data, None, args.max_files)
        rows['mc'] = discover(pre_mc, post_mc, truth_dir, args.max_files)
        candidate = args.beam_map or PROJECT / 'analysis/result/trackFit/beamsize_ssa' / beam / 'e-' / (energy_dir + '.root')
        if candidate.is_file():
            beam_map = candidate
        else:
            require(args.beam_map is None, f'Beam template missing: {candidate}')
        require(not (record / 'configuration.json').exists(), 'This energy is already initialized; use --plots-only or a different label')
        for sample in rows:
            require(not (base / sample / 'diagnostics' / energy_dir / 'diagnostics.root').exists(), 'Output exists; choose another label')
    else:
        require((record / 'complete.json').is_file(), 'Aggregation is not complete; cannot replot')
    plan = {'base': str(base), 'energy': args.energy, 'axis_layers': args.axis_layers,
            'max_events': args.max_events, 'max_files': args.max_files, 'seed': args.seed,
            'truth': str(truth_dir) if truth_dir else None, 'beam_template': str(beam_map) if beam_map else None,
            'files': {s: [[str(p) if p else '-' for p in row] for row in rr] for s, rr in rows.items()}}
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return
    record.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    with (record / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        require(args.plots_only or not (record / 'configuration.json').exists(), 'Output was initialized by another process')
        with (HERE / '.build.lock').open('a') as build_lock:
            fcntl.flock(build_lock, fcntl.LOCK_EX)
            execute(['make', '-C', HERE, 'all'], logs / 'build.log')
        # Client/user timezone, independent of execution host timezone.
        now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
        with (record / 'run.log').open('a') as log:
            log.write(now + '\n' + shlex.join([sys.executable, str(HERE / 'run.py'), *sys.argv[1:]]) + '\n')
        if not args.plots_only:
            snapshot(record / 'sources', [HERE / n for n in SOURCES])
            snapshot(record / 'geometry', [PROJECT / 'analysis/calibration_parameters/EBUdecode.cxx', PROJECT / 'analysis/calibration_parameters/EBUdecode.h', PROJECT / 'analysis/channel_response/root.mk'])
            provenance = [DATA / 'analysis/ECAL_Analysis_LCIO/src/StripSplit.cxx', MC.parent / 'Analysis_edit/src/MCDigi.cxx', MC.parent / 'CEPCScECAL/src/PrimaryGeneratorAction.cc', MC.parent / 'CEPCScECAL/generate/mac/e-' / energy_dir / (energy_dir + '_0.mac')]
            snapshot(record / 'current_upstream_sources', [p for p in provenance if p.is_file()])
            plan['input_metadata'] = {s: [[metadata(p) if p else None for p in row] for row in rr] for s, rr in rows.items()}
            plan['beam_template_metadata'] = metadata(beam_map) if beam_map else None
            plan['beam_template_sha256'] = sha(beam_map) if beam_map else None
            plan['input_integrity'] = 'size and mtime; no full ROOT checksum; selected event IDs and parent signed energies checked'
            plan['upstream_source_caveat'] = 'Current source snapshots do not prove the historical producer version.'
            plan['root_version'] = subprocess.check_output(['root-config', '--version'], text=True).strip()
            (record / 'configuration.json').write_text(json.dumps(plan, indent=2) + '\n')
            for sample, rr in rows.items():
                manifest = record / (sample + '.tsv')
                manifest.write_text(''.join('\t'.join(str(p) if p else '-' for p in row) + '\n' for row in rr))
                cmd = [HERE / 'Aggregate', '--manifest', manifest, '--out', base, '--sample', sample,
                       '--energy', args.energy, '--max-events', args.max_events, '--seed', args.seed,
                       '--axis-layers', args.axis_layers]
                if sample == 'mc' and beam_map:
                    cmd += ['--beam-map', beam_map]
                execute(cmd, logs / (sample + '.log'))
            for sample, rr in rows.items():
                require([[metadata(p) if p else None for p in row] for row in rr] == plan['input_metadata'][sample], 'An input changed during aggregation')
            if beam_map:
                require(metadata(beam_map) == plan['beam_template_metadata'] and sha(beam_map) == plan['beam_template_sha256'], 'Beam template changed during aggregation')
            (record / 'complete.json').write_text(json.dumps({'aggregation': 'complete', 'energy': args.energy}) + '\n')
        snapshot(record / 'plot_sources', [HERE / n for n in ['Common.hh', 'Plot.cc', 'run.py']])
        execute([HERE / 'Plot', base, args.energy], logs / 'plot.log')
        print(f'Results: {base}\nOverview: {base / "comparison/overview" / energy_dir}\nLogs: {logs}', flush=True)

if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, OSError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
