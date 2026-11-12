#!/usr/bin/env python3
"""Evaluate one saved MC beam condition against data; no generation or recalibration."""
import argparse
import csv
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
sys.dont_write_bytecode = True
import numpy as np
import uproot
from evaluate import evaluate, json_safe, TOLERANCES, refresh_ssa_profile_bins, SSA_PROFILE_BINNING
from threshold import add_diagnostics, parse_cuts, compare_conditions

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DATA = Path('/home/murata_t/data_beamtest')
MC = DATA/'simulation/CEPCScECAL_SML_Portable_update_new/Result_MC'
ENERGIES = ['0.5','1','2','3','4','5','10','20','30','40','50','60','70','80','100','120','150','200','250']

def require(ok, message):
    if not ok: raise RuntimeError(message)

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def metadata(p):
    p = Path(p).resolve(strict=True); s = p.stat()
    return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)

def inputs(pre, post, truth, max_files, max_events, seed, phase):
    files = sorted(pre.glob('*.root'))
    require(files, f'No calibration ROOT files: {pre}')
    require(all((post/f.name).is_file() for f in files),f'Missing paired SSA files: {post}')
    if max_files and len(files)>max_files:
        files = [files[(2*i+1)*len(files)//(2*max_files)] for i in range(max_files)]
    pairs, pool = [], []
    for i, f in enumerate(files):
        ps = [f,post/f.name,truth/f.name if truth else None]
        require(not truth or ps[2].is_file(),f'Missing matching MC truth: {ps[2]}')
        require(all('\t' not in str(p) and '\n' not in str(p) for p in ps if p),'Tabs/newlines in input paths are unsupported')
        with uproot.open(f) as root: n = root['Calib_Hit'].num_entries
        entries = np.arange(n,dtype=np.int64)
        if phase != 'all':
            half = np.array([hashlib.blake2b(f'{f.name}:{j}'.encode(),digest_size=8).digest()[0]&1 for j in entries],dtype=bool)
            entries = entries[half if phase=='validation' else ~half]
        pool.extend((i,int(e)) for e in entries)
        pairs.append(dict(paths=[str(p.resolve()) if p else '-' for p in ps],entries=[],metadata=[metadata(p) if p else None for p in ps]))
    require(pool,'No events in the requested phase')
    target = min(max_events,len(pool)) if max_events else len(pool)
    rng = np.random.default_rng(seed)
    # Strata over the complete eligible pool; never just read the first run/events.
    for i in range(target):
        j = int(rng.integers(i*len(pool)//target,(i+1)*len(pool)//target))
        fi, entry = pool[j]; pairs[fi]['entries'].append(entry)
    return pairs

def write_manifest(p, pairs):
    p.write_text(''.join('\t'.join(row['paths'])+'\t'+','.join(map(str,row['entries']))+'\n' for row in pairs))

def execute(cmd, log):
    print(shlex.join(map(str,cmd)),flush=True)
    with log.open('a') as f:
        f.write('\nCOMMAND '+shlex.join(map(str,cmd))+'\n'); f.flush()
        subprocess.run(list(map(str,cmd)),check=True,stdout=f,stderr=subprocess.STDOUT)

def unchanged(pairs):
    return all([metadata(p) if p!='-' else None for p in row['paths']]==row['metadata'] for row in pairs)

def reference_matches(stored, current):
    # Plot/evaluation changes do not invalidate the event extraction cache.
    stored = dict(stored)
    stored['source_hashes'] = {k: stored['source_hashes'].get(k) for k in current['source_hashes']}
    return stored == current

def save_json(p, value): p.write_text(json.dumps(json_safe(value),indent=2,allow_nan=False)+'\n')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energy',choices=ENERGIES,default='100')
    parser.add_argument('--label',default='beam_scan',help='comparison group; default: beam_scan')
    parser.add_argument('--condition',default='default',help='MC candidate name; current production: default')
    parser.add_argument('--phase',choices=['scan','validation','all'],default='scan')
    parser.add_argument('--result-dir',type=Path,default=PROJECT/'analysis/result/beam_quantitative')
    parser.add_argument('--max-files',type=int,default=4)
    parser.add_argument('--max-events',type=int,default=1500,help='per sample after phase partition, before axis selection; 0=all')
    parser.add_argument('--axis-layers',type=int,default=4)
    parser.add_argument('--ssa-selection',choices=['split','all'],default='split')
    parser.add_argument('--bootstrap',type=int,default=500)
    parser.add_argument('--seed',type=int,default=20261011)
    parser.add_argument('--min-events',type=int,default=200)
    parser.add_argument('--gate-layers',default='0,1,2,3')
    parser.add_argument('--tolerances-json',type=Path)
    parser.add_argument('--truth',choices=['auto','required','off'],default='auto')
    for name in ['data-calib-dir','data-ssa-dir','mc-calib-dir','mc-ssa-dir','mc-truth-dir']:
        parser.add_argument('--'+name,type=Path)
    parser.add_argument('--material-label',default='unspecified',help='e.g. no_upstream_scint; provenance annotation, not a geometry change')
    parser.add_argument('--production-note',default='',help='generation version/seed range; independent validation production if available')
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--plots-only',action='store_true')
    parser.add_argument('--threshold-diagnostics-only',action='store_true',help='add pre-SSA diagnostics to an existing completed evaluation')
    parser.add_argument('--skip-threshold-diagnostics',action='store_true',help='run only the original beam evaluation')
    parser.add_argument('--mip-cuts',default='0.5,1,2',help='additional pre-SSA signal cuts in MIP; uncut is always included')
    parser.add_argument('--mip-energy-mev',type=float,default=0.305)
    parser.add_argument('--mc-threshold-scale',type=float,help='provenance annotation for an already digitized MC; does not modify thresholds')
    parser.add_argument('--mc-digitization-log-dir',type=Path,help='matching digitization stdout directory; inferred for the current MC layout')
    parser.add_argument('--mc-threshold-file',type=Path,help='threshold parameter ROOT used for this MC, if known')
    parser.add_argument('--compare-thresholds',help='comma-separated completed diagnostic conditions; requires the same generation input')
    parser.add_argument('--threshold-comparison-name',default='threshold_scan')
    args = parser.parse_args()
    cuts = parse_cuts(args.mip_cuts)
    require(np.isfinite(args.mip_energy_mev) and args.mip_energy_mev>0,'Invalid MIP energy')
    require(args.mc_threshold_scale is None or (np.isfinite(args.mc_threshold_scale) and args.mc_threshold_scale>0),'Invalid threshold scale')
    for value in [args.label,args.condition]: require(re.fullmatch('[A-Za-z0-9_-]+',value),'Invalid label/condition')
    require(args.condition != 'reference','The condition name reference is reserved for the shared data cache')
    require(args.max_files>=0 and args.max_events>=0 and 2<=args.axis_layers<=30,'Invalid input/axis limits')
    require(args.bootstrap>=50 and args.min_events>=20 and 0<=args.seed<2**32,'Invalid bootstrap/min-events/seed')
    layers = list(map(int,args.gate_layers.split(',')))
    require(layers and len(layers)==len(set(layers)) and all(0<=v<30 for v in layers),'Invalid gate layers')
    tolerances = json.loads(args.tolerances_json.read_text()) if args.tolerances_json else {}
    require(isinstance(tolerances,dict) and not set(tolerances)-set(TOLERANCES),'Invalid tolerance dictionary')
    require(all(isinstance(v,(int,float)) and np.isfinite(v) and v>0 for v in tolerances.values()),'Invalid tolerances')
    base = args.result_dir.resolve()/args.label
    energy = args.energy+'GeV'
    record = base/'inputs'/args.phase/energy/args.condition
    shared = base/'inputs'/args.phase/energy/'reference'
    summary_dir = base/'summary'/args.phase/energy
    summary_file = summary_dir/(args.condition+'.json')
    logs = PROJECT/'analysis/log/beam_quantitative'/args.label/args.phase/args.condition/energy
    plots = base/'plots'/args.phase/args.condition/energy
    roots = {sample:base/'root'/sample/args.phase/('reference' if sample=='data' else args.condition)/energy for sample in ('data','mc')}
    comparison = base/'root/comparison'/args.phase/args.condition/energy
    if args.compare_thresholds:
        require(not args.threshold_diagnostics_only and not args.plots_only,'Use compare-thresholds by itself; it replots saved diagnostics')
        compare_conditions(HERE,base,args.label,args.phase,args.energy,args.compare_thresholds,
                           args.threshold_comparison_name,args.dry_run,execute)
        return
    if args.threshold_diagnostics_only:
        require(summary_file.is_file() and json.loads(summary_file.read_text()).get('complete'),'No completed beam evaluation to augment')
        plan = json.loads((record/'configuration.json').read_text())
        require((plan['label'],plan['condition'],plan['phase'],plan['energy']) == (args.label,args.condition,args.phase,args.energy),'Saved configuration identity mismatch')
        if args.dry_run:
            print(json.dumps(add_diagnostics(HERE,base,plan,cuts,args.mip_energy_mev,args.bootstrap,args.seed,args.min_events,
                             args.mc_threshold_scale,args.mc_threshold_file,args.mc_digitization_log_dir,True,args.plots_only,execute,unchanged),indent=2))
            return
        with (base/'inputs'/args.phase/energy/'.run.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            logs.mkdir(parents=True,exist_ok=True)
            with (HERE/'.build.lock').open('a') as buildlock:
                fcntl.flock(buildlock,fcntl.LOCK_EX);execute(['make','-C',HERE,'all'],logs/'build.log')
            add_diagnostics(HERE,base,plan,cuts,args.mip_energy_mev,args.bootstrap,args.seed,args.min_events,
                            args.mc_threshold_scale,args.mc_threshold_file,args.mc_digitization_log_dir,False,args.plots_only,execute,unchanged)
        print(f'Threshold diagnostics: {summary_dir/(args.condition+".threshold.json")}\nPlots: {plots}',flush=True)
        return
    if args.plots_only:
        require(summary_file.is_file(),'No completed evaluation to replot')
        saved = json.loads(summary_file.read_text())
        require(saved.get('complete'),'Evaluation is incomplete')
        if args.dry_run: print(json.dumps(dict(replot=str(summary_file)),indent=2)); return
        logs.mkdir(parents=True,exist_ok=True)
        with (HERE/'.build.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX); execute(['make','-C',HERE,'all'],logs/'build.log')
        with (base/'inputs'/args.phase/energy/'.run.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            rebinned={sample:refresh_ssa_profile_bins(saved[sample+'_histograms'],roots[sample]/'events.root') for sample in ('data','mc')}
            plot_sources = record/'plot_sources'; plot_sources.mkdir(exist_ok=True)
            plot_files=[HERE/'Plot.cc',HERE/'PlotCommon.hh',HERE/'run.py',HERE/'evaluate.py']
            for p in plot_files:shutil.copy2(p,plot_sources/p.name)
            save_json(plot_sources/'checksums.json',{p.name:digest(p) for p in plot_files})
            save_json(plot_sources/'profile_binning.json',SSA_PROFILE_BINNING)
            with (record/'run.log').open('a') as log:
                log.write('SSA_PROFILE_BINNING '+json.dumps(dict(binning=SSA_PROFILE_BINNING,rebinned=rebinned))+'\n')
            with (record/'run.log').open('a') as log: log.write('REPLOT '+shlex.join([sys.executable,str(HERE/'run.py'),*sys.argv[1:]])+'\n')
            execute([HERE/'Plot',saved['data_histograms'],saved['mc_histograms'],
                     summary_dir/(args.condition+'.metrics.tsv'),plots,comparison/'comparison.root',args.energy,args.condition],logs/'plot.log')
        return
    require(not record.exists() and not summary_file.exists(),'Output already initialized; use --plots-only or a new condition/label')
    beam = 'ps' if float(args.energy)<=5 else 'sps'
    db = DATA/'ECAL_data/analysed/2023'/beam
    truth = args.mc_truth_dir or MC/'generate/e-/e_ssa'/energy
    if args.truth=='off': truth=None
    elif not truth.is_dir() or not list(truth.glob('*.root')):
        require(args.truth=='auto' and args.mc_truth_dir is None,f'MC truth unavailable: {truth}')
        print('MC truth disabled: no default input',flush=True); truth=None
    pairs = dict(data=inputs(args.data_calib_dir or db/'calib/e-'/energy,args.data_ssa_dir or db/'ssaCalib/e-'/energy,None,
                            args.max_files,args.max_events,args.seed,args.phase),
                 mc=inputs(args.mc_calib_dir or MC/'calib/e-/e_ssa/threshold'/energy,
                           args.mc_ssa_dir or MC/'ssaCalib/e-/e_ssa/threshold'/energy,truth,
                           args.max_files,args.max_events,args.seed,args.phase))
    source_files = list(HERE.glob('*.py'))+list(HERE.glob('*.cc'))+list(HERE.glob('*.hh'))+[HERE/'Makefile',HERE/'run.sh',HERE/'submit.sh',HERE/'README.md',
                    HERE/'../ssa_shower_shape/Common.hh',HERE/'../calibration_parameters/EBUdecode.cxx',HERE/'../calibration_parameters/EBUdecode.h']
    source_hashes = {str(p.resolve()):digest(p) for p in source_files}
    dataplan = dict(pairs=pairs['data'],axis_layers=args.axis_layers,ssa_selection=args.ssa_selection,source_hashes={k:v for k,v in source_hashes.items() if Path(k).name in ['Extract.cc','Common.hh','EBUdecode.cxx','EBUdecode.h']})
    require(not (shared/'configuration.json').exists() or reference_matches(json.loads((shared/'configuration.json').read_text()),dataplan),
            'Shared data reference differs; choose a new label')
    plan = dict(energy=args.energy,label=args.label,condition=args.condition,phase=args.phase,
                axis_layers=args.axis_layers,ssa_selection=args.ssa_selection,max_events=args.max_events,max_files=args.max_files,
                bootstrap=args.bootstrap,seed=args.seed,min_events=args.min_events,gate_layers=layers,
                tolerances=dict(TOLERANCES,**tolerances),material_label=args.material_label,production_note=args.production_note,
                pairs=pairs,source_hashes=source_hashes,result=str(base),logs=str(logs),
                phase_rule='Stable filename:entry hash halves. Validation MC is not automatically a new seed production.',
                threshold_diagnostics=dict(enabled=not args.skip_threshold_diagnostics,cuts=cuts,mip_energy_mev=args.mip_energy_mev,mc_threshold_scale=args.mc_threshold_scale,mc_threshold_file=str(args.mc_threshold_file) if args.mc_threshold_file else None),
                input_integrity='Selected IDs/energy checked; input size/mtime recorded, no full ROOT hash',
                geometry_caveat='Current source snapshots are not proof of the producer version.')
    if args.dry_run:
        print(json.dumps(dict(plan,selected_events={s:sum(len(r['entries']) for r in rr) for s,rr in pairs.items()}),indent=2)); return
    lockdir = base/'inputs'/args.phase/energy
    lockdir.mkdir(parents=True,exist_ok=True)
    with (lockdir/'.run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        require(not record.exists(),'Output initialized by another process')
        record.mkdir(); shared.mkdir(exist_ok=True)
        logs.mkdir(parents=True,exist_ok=True); summary_dir.mkdir(parents=True,exist_ok=True)
        plots.mkdir(parents=True,exist_ok=True); comparison.mkdir(parents=True,exist_ok=True)
        for p in roots.values(): p.mkdir(parents=True,exist_ok=True)
        save_json(record/'configuration.json',plan)
        (record/'run.log').write_text(datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()+'\n'+shlex.join([sys.executable,str(HERE/'run.py'),*sys.argv[1:]])+'\n')
        snapshots = record/'sources'; snapshots.mkdir()
        for i,p in enumerate(source_files): shutil.copy2(p,snapshots/(str(i)+'_'+p.name))
        upstream = MC.parent
        for rel in ['CEPCScECAL/src/DetectorConstruction.cc','CEPCScECAL/src/PrimaryGeneratorAction.cc',
                    'Analysis_edit/src/MCDigi.cxx']:
            p=upstream/rel
            if p.is_file(): shutil.copy2(p,snapshots/p.name)
        with (HERE/'.build.lock').open('a') as buildlock:
            fcntl.flock(buildlock,fcntl.LOCK_EX)
            execute(['make','-C',HERE,'all'],logs/'build.log')
        if not (shared/'complete.json').exists():
            require(not (roots['data']/'events.root').exists(),'Incomplete shared data extraction exists; choose new label')
            save_json(shared/'configuration.json',dataplan); write_manifest(shared/'data.tsv',pairs['data'])
            execute([HERE/'Extract',shared/'data.tsv',roots['data']/'events.root',args.energy,args.axis_layers,args.ssa_selection,'data'],logs/'data.log')
            require(unchanged(pairs['data']),'Data input changed during extraction')
            save_json(shared/'complete.json',dict(extraction='complete'))
        else:
            require(reference_matches(json.loads((shared/'configuration.json').read_text()),dataplan),'Shared data changed while waiting')
        write_manifest(record/'mc.tsv',pairs['mc'])
        execute([HERE/'Extract',record/'mc.tsv',roots['mc']/'events.root',args.energy,args.axis_layers,args.ssa_selection,'mc'],logs/'mc.log')
        require(all(unchanged(rr) for rr in pairs.values()),'Input changed during extraction')
        # Candidate-specific data histograms avoid simultaneous writes to a shared reference.
        data_hist = comparison/'data_histograms.root'
        result = evaluate(roots['data']/'events.root',roots['mc']/'events.root',data_hist,roots['mc']/'histograms.root',
                          comparison/'metrics.root',summary_dir/(args.condition+'.metrics.tsv'),
                          args.bootstrap,args.seed,args.min_events,layers,tolerances)
        execute([HERE/'Plot',data_hist,roots['mc']/'histograms.root',summary_dir/(args.condition+'.metrics.tsv'),
                 plots,comparison/'comparison.root',args.energy,args.condition],logs/'plot.log')
        require(all(unchanged(rr) for rr in pairs.values()),'Input changed during evaluation')
        result.update(energy=args.energy,condition=args.condition,label=args.label,phase=args.phase,complete=True,
                      data_histograms=str(data_hist),mc_histograms=str(roots['mc']/'histograms.root'),
                      evaluation_signature=dict(tolerances=plan['tolerances'],gate_layers=layers,ssa_selection=args.ssa_selection,
                                                axis_layers=args.axis_layers,material_label=args.material_label))
        save_json(summary_file,result)
        candidates=[json.loads(p.read_text()) for p in summary_dir.glob('*.json')]
        candidates=[r for r in candidates if r.get('complete') and 'beam_status' in r]
        candidates=[r for r in candidates if r.get('evaluation_signature')==result['evaluation_signature']]
        order={'PASS':0,'UNDETERMINED':1,'FAIL':2}
        candidates.sort(key=lambda r:(order[r['beam_status']],r['score'] if r['score'] is not None else float('inf')))
        with (summary_dir/'ranking.tsv').open('w',newline='') as f:
            writer=csv.writer(f,delimiter='\t');writer.writerow(['rank','condition','beam_status','score','all_layer_lab_status','aligned_diagnostic_status','data_events','mc_events'])
            for i,r in enumerate(candidates,1):writer.writerow([i,r['condition'],r['beam_status'],r['score'],r['all_layer_lab_status'],r['aligned_diagnostic_status'],r['data']['accepted'],r['mc']['accepted']])
        if not args.skip_threshold_diagnostics:
            add_diagnostics(HERE,base,plan,cuts,args.mip_energy_mev,args.bootstrap,args.seed,args.min_events,
                            args.mc_threshold_scale,args.mc_threshold_file,args.mc_digitization_log_dir,False,False,execute,unchanged)
        print(f'Result: {summary_file}\nBeam assessment: {result["beam_status"]}; score={result["score"]}\nPlots: {plots}\nLogs: {logs}',flush=True)

if __name__=='__main__':
    try: main()
    except (RuntimeError,OSError,ValueError,KeyError,subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}',file=sys.stderr);sys.exit(1)
