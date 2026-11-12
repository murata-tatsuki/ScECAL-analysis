#!/usr/bin/env python3
"""Compare MC stages with a common threshold-quality channel mask and fixed beam cohort."""
import argparse
import contextlib
import fcntl
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
import numpy as np
import uproot

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
sys.path.insert(0,str(HERE.parent/'beam_quantitative'))
from metrics import evaluate, has_sipm_histograms, ensure_sipm_histograms
from threshold_quality import add_arguments, select, write_table


def require(ok,message):
    if not ok:raise RuntimeError(message)


def metadata(path):
    p=Path(path).resolve(strict=True);s=p.stat()
    return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def save_json(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


def execute(cmd,log):
    print(shlex.join(map(str,cmd)),flush=True)
    with log.open('a') as f:
        f.write('\nCOMMAND '+shlex.join(map(str,cmd))+'\n');f.flush()
        subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,check=True)


def same_inputs(plan):
    for row in plan['pairs']:
        require(list(map(metadata,row['paths']))==row['metadata'],'Stage input changed during analysis')
    for row in plan['beam_cache']:
        require(metadata(row['path'])=={k:row[k] for k in ('path','bytes','mtime_ns')},'Beam cache changed during analysis')
    require(metadata(plan['beam_configuration']['path'])==plan['beam_configuration'],'Beam configuration changed during analysis')


@contextlib.contextmanager
def matched_cache(plan,common_record,common_roots,beam_roots,log):
    """Shared readers run together; a changed extraction waits for all readers."""
    common_record.mkdir(parents=True,exist_ok=True);common_roots.mkdir(parents=True,exist_ok=True)
    events=common_roots/'events.root';stamp=common_record/'cache_configuration.json'
    extractor_sources=[HERE/'Extract.cc',HERE.parent/'ssa_shower_shape/Common.hh',
                       HERE.parent/'calibration_parameters/EBUdecode.cxx',HERE.parent/'calibration_parameters/EBUdecode.h']
    source_hashes={str(p):digest(p) for p in extractor_sources}
    signature=dict(pairs=plan['pairs'],axis_layers=plan['axis_layers'],beam_cache=plan['beam_cache'],source_hashes=source_hashes)
    def ready():
        if not events.is_file() or not stamp.is_file():return False
        old=json.loads(stamp.read_text())
        return old['signature']==signature and metadata(events)==old['events'] and digest(events)==old['events_sha256']
    with (common_record/'.extract.lock').open('a+') as guard:
        fcntl.flock(guard,fcntl.LOCK_SH)
        if not ready():
            fcntl.flock(guard,fcntl.LOCK_UN);fcntl.flock(guard,fcntl.LOCK_EX)
            if not ready():
                # Adopt the unchanged original allChannel extraction after the
                # directory migration. Plot/evaluation edits do not invalidate it.
                old_path=common_record/'configuration.json';adopt=False
                if events.is_file() and old_path.is_file() and not stamp.is_file():
                    old=json.loads(old_path.read_text())
                    adopt=(old.get('pairs')==plan['pairs'] and old.get('axis_layers')==plan['axis_layers'] and
                           old.get('beam_cache')==plan['beam_cache'] and
                           all(old.get('source_hashes',{}).get(k)==v for k,v in source_hashes.items()))
                if not adopt:
                    manifest=common_record/'mc.tsv'
                    manifest.write_text(''.join('\t'.join(r['paths'])+'\t'+','.join(map(str,r['entries']))+'\n' for r in plan['pairs']))
                    execute([HERE/'Extract',manifest,beam_roots['mc']/'events.root',beam_roots['mc']/'pre_ssa_events.root',
                             events,plan['axis_layers'],plan['energy']],log)
                same_inputs(plan)
                save_json(stamp,dict(signature=signature,events=metadata(events),events_sha256=digest(events),
                                     adopted_original_cache=adopt))
                with log.open('a') as f:f.write('MATCHED_CACHE '+('adopted original cache' if adopt else 'extracted')+' '+str(events)+'\n')
            fcntl.flock(guard,fcntl.LOCK_SH)
        require(ready(),'Matched cache changed before read')
        with uproot.open(events) as f:
            require(f['events'].num_entries==sum(len(r['entries']) for r in plan['pairs']),'Matched cache event count differs')
        yield events,json.loads(stamp.read_text())
        require(ready(),'Matched cache changed during analysis')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--energy',default='100')
    p.add_argument('--condition',default='point')
    p.add_argument('--label',default='beam_scan')
    p.add_argument('--phase',choices=['scan','validation','all'],default='scan')
    p.add_argument('--beam-result-dir',type=Path,default=PROJECT/'analysis/result/beam_quantitative')
    p.add_argument('--result-dir',type=Path,default=PROJECT/'analysis/result/strip_stage_response')
    p.add_argument('--mc-digitized-dir',type=Path,help='Raw_Hit directory including the energy directory; otherwise generate -> decode')
    p.add_argument('--bootstrap',type=int,default=500)
    p.add_argument('--seed',type=int,default=20261011)
    p.add_argument('--mip-energy-mev',type=float,default=.305)
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--plots-only',action='store_true')
    add_arguments(p);a=p.parse_args()
    for value in (a.label,a.condition):require(re.fullmatch('[A-Za-z0-9_-]+',value),'Invalid label/condition')
    require(a.condition!='reference','reference is reserved')
    require(re.fullmatch('[0-9]+(?:\.[0-9]+)?',a.energy) and float(a.energy)>0,'Invalid energy')
    require(a.bootstrap>=50 and 0<=a.seed<2**32,'Invalid bootstrap/seed')
    require(np.isfinite(a.mip_energy_mev) and a.mip_energy_mev>0,'Invalid nominal MIP energy')
    require(re.fullmatch('[A-Za-z0-9_-]+',a.threshold_cut),'Invalid threshold-cut name')
    energy=a.energy+'GeV';base=a.result_dir.resolve()/a.label
    record=base/'inputs'/a.phase/a.condition/a.threshold_cut/energy
    roots=base/'root/mc'/a.phase/a.condition/a.threshold_cut/energy
    plots=base/'plots/mc'/a.phase/a.condition/a.threshold_cut/energy
    summary=base/'summary/mc'/a.phase/a.condition/a.threshold_cut/energy
    common_record=base/'inputs'/a.phase/a.condition/'allChannel'/energy
    common_roots=base/'root/mc'/a.phase/a.condition/'allChannel'/energy
    summary_file=summary/'summary.json'
    beam_base=a.beam_result_dir.resolve()/a.label
    beam_record=beam_base/'inputs'/a.phase/energy/a.condition
    beam_roots={s:beam_base/'root'/s/a.phase/('reference' if s=='data' else a.condition)/energy for s in ('data','mc')}
    if a.plots_only:
        require(summary_file.is_file() and json.loads(summary_file.read_text()).get('complete'),'No completed result to replot')
        require((roots/'histograms.root').is_file(),'Missing saved histograms')
        plan=dict(replot=str(roots/'histograms.root'),plots=str(plots),threshold_cut=a.threshold_cut)
    else:
        selection,channel_rows=select(a)
        config_path=beam_record/'configuration.json'
        require(config_path.is_file(),f'Run beam_quantitative for this condition/energy first: {config_path}')
        config=json.loads(config_path.read_text())
        require((config['label'],config['condition'],config['energy'],config['phase'])==(a.label,a.condition,a.energy,a.phase),'Beam identity mismatch')
        complete=beam_base/'summary'/a.phase/energy/(a.condition+'.json')
        require(complete.is_file() and json.loads(complete.read_text()).get('complete'),'Beam evaluation is incomplete')
        for sample in ('data','mc'):
            for name in ('events.root','pre_ssa_events.root'):
                require((beam_roots[sample]/name).is_file(),f'Missing beam cache: {beam_roots[sample]/name}; run beam diagnostics first')
            for pair in config['pairs'][sample]:
                require([metadata(v) if v!='-' else None for v in pair['paths']]==pair['metadata'],
                        'Beam inputs changed; rerun beam_quantitative: '+pair['paths'][0])
        pairs=[]
        for row in config['pairs']['mc']:
            calib,ssa,truth=row['paths'];require(truth!='-','MC truth is required')
            t=Path(truth);parts=list(t.parts)
            if a.mc_digitized_dir:raw=a.mc_digitized_dir.resolve()/t.name
            else:
                require('generate' in parts,'Specify --mc-digitized-dir');parts[parts.index('generate')]='decode';raw=Path(*parts)
            paths=[t.resolve(),raw.resolve(strict=True),Path(calib).resolve()]
            for path,tree in zip(paths,('MC_Truth','Raw_Hit','Calib_Hit')):
                with uproot.open(path) as f:
                    require(tree in f,f'Missing {tree}: {path}')
                    require(all(0<=i<f[tree].num_entries for i in row['entries']),f'Invalid sampled entry: {path}')
            require(all('\t' not in str(v) and '\n' not in str(v) for v in paths),'Unsupported path characters')
            pairs.append(dict(paths=list(map(str,paths)),entries=row['entries'],metadata=list(map(metadata,paths))))
        cached=[beam_roots[sample]/n for sample in ('data','mc') for n in ('events.root','pre_ssa_events.root')]
        plan=dict(schema=2,label=a.label,condition=a.condition,energy=a.energy,phase=a.phase,
                  axis_layers=config['axis_layers'],bootstrap=a.bootstrap,seed=a.seed,mip_energy_mev=a.mip_energy_mev,
                  pairs=pairs,beam_configuration=metadata(config_path),beam_configuration_sha256=digest(config_path),
                  beam_cache=[dict(**metadata(v),sha256=digest(v)) for v in cached],
                  selection='Original beam event cohort and original axes remain fixed for every channel cut',
                  channel_selection=selection,roots=str(roots),plots=str(plots),summary=str(summary))
    if a.dry_run:print(json.dumps(plan,indent=2));return
    record.mkdir(parents=True,exist_ok=True)
    with (record/'.run.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('The same condition/cut/energy is already running')
        for path in (roots,plots,summary):path.mkdir(parents=True,exist_ok=True)
        log=record/'run.log'
        with log.open('a') as f:f.write('\nINVOCATION '+shlex.join([sys.executable,*sys.argv])+'\n')
        with (HERE/'.build.lock').open('a') as build:
            fcntl.flock(build,fcntl.LOCK_EX);execute(['make','-C',HERE,'all'],log)
        source_dir=record/('plot_sources' if a.plots_only else 'sources');source_dir.mkdir(exist_ok=True)
        sources=[HERE/n for n in ('run.py','metrics.py','Extract.cc','Plot.cc','Makefile','run.sh','submit.sh','README.md')]
        sources += [HERE.parent/'ssa_shower_shape/Common.hh',HERE.parent/'calibration_parameters/EBUdecode.cxx',
                    HERE.parent/'calibration_parameters/EBUdecode.h',HERE.parent/'beam_quantitative/PlotCommon.hh',HERE.parent/'beam_quantitative/evaluate.py',HERE/'threshold_quality.py']
        hashes={}
        for i,src in enumerate(sources):
            shutil.copy2(src,source_dir/(str(i)+'_'+src.name));hashes[str(src)]=digest(src)
        save_json(source_dir/'checksums.json',hashes)
        if a.plots_only:
            report=json.loads(summary_file.read_text())
            hist_path=roots/'histograms.root'
            if not has_sipm_histograms(hist_path):
                # A saved pooled histogram cannot be split retrospectively.
                # Read the original sampled events and the saved channel mask.
                saved=json.loads((record/'configuration.json').read_text())
                require(a.threshold_cut=='allChannel' or report.get('channel_selection') is not None,
                        'Missing saved channel selection for SiPM replot')
                events=saved.get('matched_events')
                if events is None:
                    stamp=json.loads((common_record/'cache_configuration.json').read_text())
                    events=dict(**stamp['events'],sha256=stamp['events_sha256'])
                data_beam,data_strip=saved['beam_cache'][:2]
                require(Path(data_beam['path']).name=='events.root' and
                        Path(data_strip['path']).name=='pre_ssa_events.root' and
                        '/root/data/' in data_beam['path'],'Invalid saved data cache order')
                def check_caches():
                    for row in (events,data_beam,data_strip):
                        require(metadata(row['path'])=={k:row[k] for k in ('path','bytes','mtime_ns')} and
                                digest(row['path'])==row['sha256'],
                                'Saved event cache changed; rerun analysis before replotting: '+row['path'])
                with (common_record/'.extract.lock').open('a+') as guard:
                    fcntl.flock(guard,fcntl.LOCK_SH);check_caches()
                    sipm=ensure_sipm_histograms(hist_path,events['path'],data_strip['path'],data_beam['path'],report,validate=check_caches)
                    check_caches()
                report['sipm']=sipm;save_json(summary_file,report)
                with log.open('a') as f:f.write('SIPM_DISPLAY_ADDED from original event caches and saved channel mask\n')
            execute([HERE/'Plot',hist_path,plots,roots/'plots.root',a.energy,a.condition,a.threshold_cut],log)
        else:
            save_json(record/'complete.json',dict(complete=False));save_json(summary_file,dict(complete=False))
            plan['source_hashes']=hashes
            if record!=common_record:save_json(record/'configuration.json',plan)
            shutil.copy2(config_path,record/'beam_configuration.json')
            write_table(record/'threshold_channels.tsv',channel_rows)
            with matched_cache(plan,common_record,common_roots,beam_roots,log) as (events,cache):
                plan['matched_events']=dict(**cache['events'],sha256=cache['events_sha256'])
                save_json(record/'configuration.json',plan)
                report=evaluate(events,beam_roots['data']/'pre_ssa_events.root',beam_roots['data']/'events.root',
                                roots/'histograms.root',summary_file,summary/'metrics.tsv',a.bootstrap,a.seed,a.mip_energy_mev,selection)
                execute([HERE/'Plot',roots/'histograms.root',plots,roots/'plots.root',a.energy,a.condition,a.threshold_cut],log)
                same_inputs(plan)
                for row in plan['beam_cache']:require(digest(row['path'])==row['sha256'],'Beam cache changed during evaluation')
                if selection['source'] is not None:
                    src=selection['source'];require(metadata(src['path'])=={k:src[k] for k in ('path','bytes','mtime_ns')},'Threshold fit input changed during analysis')
            report['complete']=True;save_json(summary_file,report);save_json(record/'complete.json',dict(complete=True))
            print(f"Channels: {selection['selected_channels']}/6300 | cut: {a.threshold_cut}",flush=True)
    print(f'Summary: {summary_file}\nPlots: {plots}',flush=True)


if __name__=='__main__':
    try:main()
    except Exception as e:print(f'ERROR: {e}',file=sys.stderr);sys.exit(1)
