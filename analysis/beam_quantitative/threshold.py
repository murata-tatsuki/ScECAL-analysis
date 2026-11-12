"""Additional pre-SSA / MIP-cut diagnostics; never alter the original beam verdict."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import re

import awkward as ak
import numpy as np
import uproot
from evaluate import (buffered_root, json_safe, matrix, relative,
                      resamples, statistics)

SCHEMA = 1
SIGNAL_EDGES = np.linspace(-2., 50., 521)
# EBUdecode: 42 strips across the fine axis, 5 along the strip-length axis.
# Contiguous cells include the 0.3/0.4 mm inter-strip gaps in their pitch.
STRIP_FINE_EDGES = (np.arange(43, dtype=float) - 21.) * 5.3
STRIP_LONG_EDGES = (np.arange(6, dtype=float) - 2.5) * 45.4
STRIP_PLOT_BINNING = dict(version=1, fine_bins=42, long_bins=5,
                          fine_pitch_mm=5.3, long_pitch_mm=45.4,
                          nominal_strip_size_mm=[5., 45.],
                          even_layers='X: long axis; Y: fine axis',
                          odd_layers='X: fine axis; Y: long axis',
                          rule='One contiguous bin per physical strip; pitch includes gaps')


def strip_xy_edges(layer):
    return (STRIP_FINE_EDGES, STRIP_LONG_EDGES) if layer % 2 else (STRIP_LONG_EDGES, STRIP_FINE_EDGES)


def strip_plot_histograms(sample, cuts):
    """Only display bins: numerical assessments use the exact strip centers."""
    histograms={}
    for layer in range(30):
        axis='x' if layer % 2 else 'y'
        selected=sample['layer']==layer
        for cut in cuts:
            keep=selected if cut is None else selected & (sample['mip']>=cut)
            coord=sample[axis][keep]
            histograms['pre/'+cut_key(cut)+'/layer_'+str(layer)+'/u']=np.histogram(
                coord,bins=STRIP_FINE_EDGES,weights=np.full(len(coord),1./sample['n']))
        edges=strip_xy_edges(layer)
        histograms['occupancy/layer_'+str(layer)]=(np.histogram2d(
            sample['x'][selected],sample['y'][selected],bins=edges)[0]/sample['n'],*edges)
    return histograms


def refresh_strip_plot_bins(histogram_path, events_path, cuts, mip_energy):
    """Rebin cached events for plotting; preserve signal spectra and bootstrap curves."""
    with uproot.open(histogram_path) as source:
        if 'binning/strip_cells' in source and json.loads(str(source['binning/strip_cells']))==STRIP_PLOT_BINNING:
            return False
        retained={}
        for key,classname in source.classnames(cycle=False).items():
            if classname=='TDirectory':continue
            retained[key]=str(source[key]) if classname=='TObjString' else source[key]
    replacements=strip_plot_histograms(read_sample(events_path,mip_energy),cuts)
    with buffered_root(histogram_path) as output:
        for key,value in retained.items():
            if key not in replacements and key!='binning/strip_cells':output[key]=value
        for key,value in replacements.items():output[key]=value
        output['binning/strip_cells']=json.dumps(STRIP_PLOT_BINNING,sort_keys=True)
    return True


def remove_legacy_plot_names(plots):
    # Remove only the old name of a figure that has just been redrawn as strip_*.
    for old in plots.glob('pre_ssa_*.png'):
        current=old.with_name('strip_'+old.name[len('pre_ssa_'):])
        if current.is_file():old.unlink()


def cut_key(cut):
    return 'uncut' if cut is None else 'mip_' + format(cut, '.8g').replace('.', 'p')


def parse_cuts(value):
    cuts = [float(x) for x in value.split(',')]
    if not cuts or len(cuts) > 6 or len(set(cuts)) != len(cuts) or any(not np.isfinite(x) or x < 0 for x in cuts):
        raise RuntimeError('MIP cuts must be 1-6 distinct finite nonnegative values')
    return [None] + sorted(cuts)


def read_sample(path, mip_energy, valid_override=None):
    with uproot.open(path) as root:
        a = root['events'].arrays(library='ak')
    original_valid = np.asarray(a.axisValid, dtype=bool)
    valid = original_valid.copy()
    if valid_override is not None:
        selected=np.asarray(valid_override,dtype=bool)
        if selected.shape!=valid.shape or np.any(selected & ~original_valid):raise RuntimeError('Invalid common-event selection')
        valid &= selected
    audit = dict(sampled=len(valid), accepted=int(valid.sum()), axis_invalid=int((~original_valid).sum()), shared_selection_excluded=int((original_valid & ~valid).sum()),
                 invalid_pre_hits=int(ak.sum(a.invalidPre)), merged_entries=int(ak.sum(a.mergedEntries)))
    a = a[valid]
    n = len(a)
    if not n:
        raise RuntimeError('No axis-valid events for pre-SSA diagnostics')
    sizes = np.asarray(ak.num(a.layer), dtype=np.int32)
    return dict(n=n, audit=audit, event=np.repeat(np.arange(n, dtype=np.int32), sizes),
                layer=np.asarray(ak.flatten(a.layer), dtype=np.int32),
                cell=np.asarray(ak.flatten(a.cellID), dtype=np.int64),
                x=np.asarray(ak.flatten(a.x)), y=np.asarray(ak.flatten(a.y)),
                mip=np.asarray(ak.flatten(a.energy))/mip_energy,
                post=np.asarray(ak.to_numpy(a.postCounts), dtype=float))


def wilson(successes, n):
    z=1.959963984540054
    p=successes/n
    center=(p+z*z/(2*n))/(1+z*z/n)
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return p, np.maximum(0.,center-half), np.minimum(1.,center+half)


def write_rows(path, rows):
    with path.open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def root_table(root, name, rows):
    schema={k: 'string' if isinstance(v,str) else 'bool' if isinstance(v,bool) else 'int64' if isinstance(v,(int,np.integer)) else 'float64'
            for k,v in rows[0].items()}
    tree=root.mktree(name,schema)
    tree.extend({k: [r[k] for r in rows] if typ=='string' else np.asarray([r[k] for r in rows]) for k,typ in schema.items()})


def assess_pre(data_path, mc_path, data_hist, mc_hist, metrics_root, table_path, occupancy_path,
               cuts, mip_energy, bootstraps, seed, min_events, tolerances, audit_status="unavailable", mc_mask=None):
    d,m=read_sample(data_path,mip_energy),read_sample(mc_path,mip_energy,mc_mask)
    rng=np.random.default_rng(seed)
    bd,bm=resamples(d['n'],bootstraps,rng),resamples(m['n'],bootstraps,rng)
    rows=[]; occupancy=[]; totals=[]
    sufficient=min(d['n'],m['n'])>=min_events

    def add(layer,cut,metric,da,ma,delta,unit):
        finite=np.isfinite(delta[1:]); complete=np.isfinite(delta[0]) and finite.all()
        lo,hi=np.quantile(delta[1:],[.025,.975]) if complete else (np.nan,np.nan)
        limit=tolerances[metric]; status='UNDETERMINED'
        if complete and sufficient:
            if max(abs(lo),abs(hi))<=limit: status='PASS'
            elif lo>limit or hi < -limit: status='FAIL'
        rows.append(dict(group='pre_ssa',layer=layer,cut_key=cut_key(cut),cut_mip=np.nan if cut is None else cut,
                         metric=metric,value_unit=unit,difference_unit='mm' if metric=='mean' else 'relative',
                         data_value=float(da[0]),mc_value=float(ma[0]),difference=float(delta[0]),
                         ci_low=float(lo),ci_high=float(hi),tolerance=float(limit),status=status,
                         valid_bootstraps=int(finite.sum()),data_events=d['n'],mc_events=m['n']))

    with buffered_root(data_hist) as dh, buffered_root(mc_hist) as mh:
        stores=((dh,d,bd),(mh,m,bm))
        mh['audit/threshold_mode']=audit_status
        for store,sample,boot in stores:
            store['binning/strip_cells']=json.dumps(STRIP_PLOT_BINNING,sort_keys=True)
            store['signal/all']=np.histogram(sample['mip'],bins=SIGNAL_EDGES,weights=np.full(len(sample['mip']),1./sample['n']))
            store['signal/flow_per_event']=(np.array([np.sum(sample['mip']<SIGNAL_EDGES[0]),np.sum(sample['mip']>SIGNAL_EDGES[-1])])/sample['n'],np.arange(3.))
            post=boot@sample['post']/sample['n']
            store['post/count_mean']=(post[0],np.arange(-.5,30.5))
            store['post/count_low']=(np.quantile(post[1:],.025,axis=0),np.arange(-.5,30.5))
            store['post/count_high']=(np.quantile(post[1:],.975,axis=0),np.arange(-.5,30.5))
            # Per-channel occupancy is Bernoulli after the physical-ID merge.
            chan,counts=np.unique(sample['cell'],return_counts=True)
            for cell,count in zip(chan,counts):
                if count>sample['n']: raise RuntimeError('Physical channel is repeated within an event')
            sample['occupancy']=dict(zip(chan,counts))
        cells=[layer*100000+chip*10000+channel for layer in range(30) for chip in range(6) for channel in range(30 if chip==5 else 36)]
        for cell in cells:
            dp,dl,du=wilson(d['occupancy'].get(cell,0),d['n'])
            mp,ml,mu=wilson(m['occupancy'].get(cell,0),m['n'])
            occupancy.append(dict(cell_id=int(cell),layer=int(cell//100000),
                                  data_fraction=float(dp),data_ci_low=float(dl),data_ci_high=float(du),
                                  mc_fraction=float(mp),mc_ci_low=float(ml),mc_ci_high=float(mu),
                                  difference=float(mp-dp),data_events=d['n'],mc_events=m['n']))
        for cut in cuts:
            key=cut_key(cut); totals_for_cut=[]
            for sample,boot in ((d,bd),(m,bm)):
                selected=np.ones(len(sample['mip']),dtype=bool) if cut is None else sample['mip']>=cut
                per_event=np.bincount(sample['event'][selected],minlength=sample['n'])
                totals_for_cut.append(boot@per_event/sample['n'])
            dn,mn=totals_for_cut
            add(-1,cut,'hits_relative',dn,mn,relative(dn,mn),'strips/event')
            totals.append(dict(cut_key=key,cut_mip=cut,data_hits=float(dn[0]),mc_hits=float(mn[0]),
                               difference=float(relative(dn,mn)[0])))
            for layer in range(30):
                axis='x' if layer%2 else 'y'
                ds=[];ms=[]
                for store,sample,boot in stores:
                    selected=sample['layer']==layer
                    if cut is not None: selected &= sample['mip']>=cut
                    ev=sample['event'][selected]; coord=sample[axis][selected]
                    values=np.unique(np.round(coord,6))
                    if len(values):
                        p=matrix(sample['n'],ev,np.searchsorted(values,np.round(coord,6)),len(values),np.ones(len(ev)))
                        stats=statistics(values,boot@p)
                    else:
                        stats={k:np.full(boot.shape[0],np.nan) for k in ('mean','rms','q90')}
                    counts=np.bincount(ev,minlength=sample['n']); means=boot@counts/sample['n']
                    stats['count']=means
                    (ds if sample is d else ms).append(stats)
                    edges=np.arange(-.5,30.5)
                    # Save one layer-at-a-time then assemble layer curves below.
                    prefix='pre/'+key+'/layer_'+str(layer)
                    store[prefix+'/u']=np.histogram(coord,bins=STRIP_FINE_EDGES,weights=np.full(len(coord),1./sample['n']))
                    if cut is None:
                        store['signal/layer_'+str(layer)]=np.histogram(sample['mip'][selected],bins=SIGNAL_EDGES,
                                                                      weights=np.full(selected.sum(),1./sample['n']))
                        xy_edges=strip_xy_edges(layer)
                        store['occupancy/layer_'+str(layer)]=(np.histogram2d(sample['x'][selected],sample['y'][selected],bins=xy_edges)[0]/sample['n'],*xy_edges)
                    sample.setdefault('curves',{})[key,layer]=(means,stats['rms'])
                da,ma=ds[0],ms[0]
                add(layer,cut,'hits_relative',da['count'],ma['count'],relative(da['count'],ma['count']),'strips/event')
                add(layer,cut,'mean',da['mean'],ma['mean'],ma['mean']-da['mean'],'mm')
                add(layer,cut,'rms_relative',da['rms'],ma['rms'],relative(da['rms'],ma['rms']),'mm')
            for store,sample,boot in stores:
                for index,name in ((0,'count'),(1,'rms')):
                    curve=np.column_stack([sample['curves'][key,layer][index] for layer in range(30)])
                    store['pre/'+key+'/'+name+'_mean']=(curve[0],np.arange(-.5,30.5))
                    for q,part in ((.025,'low'),(.975,'high')):
                        low=np.array([np.quantile(v[1:],q) if np.isfinite(v[1:]).all() else np.nan for v in curve.T])
                        store['pre/'+key+'/'+name+'_'+part]=(low,np.arange(-.5,30.5))
                if cut is None:
                    pre=np.column_stack([sample['curves'][key,l][0] for l in range(30)])
                    with np.errstate(divide='ignore',invalid='ignore'): multiple=(boot@sample['post']/sample['n'])/pre
                    store['ssa_multiplicity/mean']=(multiple[0],np.arange(-.5,30.5))
        mh['comparison/count_tolerance']=(np.array([tolerances['hits_relative']]),np.array([0.,1.]))
        for cut in cuts:
            rr=[r for r in rows if r['cut_key']==cut_key(cut) and r['layer']>=0 and r['metric']=='hits_relative']
            mh['comparison/'+cut_key(cut)+'/count_difference']=(np.array([r['difference'] for r in rr]),np.arange(-.5,30.5))
            mh['comparison/'+cut_key(cut)+'/count_low']=(np.array([r['ci_low'] for r in rr]),np.arange(-.5,30.5))
            mh['comparison/'+cut_key(cut)+'/count_high']=(np.array([r['ci_high'] for r in rr]),np.arange(-.5,30.5))
    write_rows(table_path,rows);write_rows(occupancy_path,occupancy)
    with buffered_root(metrics_root) as root:
        root_table(root,'metrics',rows);root_table(root,'channel_occupancy',occupancy)
    return json_safe(dict(schema=SCHEMA,complete=True,data=d['audit'],mc=m['audit'],mip_energy_mev=mip_energy,
                          cuts=cuts,bootstrap_replicates=bootstraps,min_events=min_events,
                          confidence_method='95% percentile event bootstrap; channel occupancy uses Wilson intervals',
                          metric_count=len(rows),totals=totals,
                          interpretation='Offline cuts on saved pre-SSA strips diagnose low-signal sensitivity; missing hardware hits cannot be restored. SSA is not rerun by an offline cut.'))


def audit_digitization(plan, log_dir=None):
    # Filename-matched logs are evidence, not a proof of the exact producer.
    if log_dir is None:
        first=Path(plan['pairs']['mc'][0]['paths'][0])
        try:
            parts=first.parts;index=parts.index('Result_MC')
            relative=Path(*parts[index+2:-1])
            log_dir=Path(*parts[:index])/'Analysis_edit/job/scan'/relative/'out'
        except ValueError:
            return dict(status='unavailable',logs=[],interpretation='No production log directory inferred')
    log_dir=Path(log_dir)
    reports=[]
    for row in plan['pairs']['mc']:
        pre=Path(row['paths'][0])
        matches=list(log_dir.glob(pre.name+'-*.out'))
        matches=[p for p in matches if p.stat().st_mtime_ns <= row['metadata'][0]['mtime_ns']]
        if not matches:continue
        p=max(matches,key=lambda f:f.stat().st_mtime_ns)
        content=p.read_text(errors='replace')
        fallback='GET Threshold FILE FAILED' in content
        value=re.search(r'^threshold ([-+0-9.eE]+)\s*$',content,re.M)
        parameter=re.search(r'^threshold file (.+)$',content,re.M)
        reports.append(dict(path=str(p.resolve()),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,
                            sha256=hashlib.sha256(p.read_bytes()).hexdigest(),mc_file=pre.name,
                            mode='fallback' if fallback else 'map_loaded' if 'mean of Threshold' in content else 'unknown',
                            fallback_adc=float(value.group(1)) if fallback and value else None,
                            parameter_path=parameter.group(1) if parameter else None))
    status='fallback_detected' if any(r['mode']=='fallback' for r in reports) else 'map_loaded' if len(reports)==len(plan['pairs']['mc']) and all(r['mode']=='map_loaded' for r in reports) else 'unavailable_or_partial'
    return dict(status=status,logs=reports,selected_mc_files=len(plan['pairs']['mc']),
                interpretation='Latest filename-matched log ending before calibration mtime; producer provenance is not cryptographically linked. Fallback keeps pedestal/calibration cuts.')


def threshold_map(path):
    if path is None: return None
    path=Path(path).resolve(strict=True)
    with uproot.open(path) as root:
        arrays=root['expErfThre'].arrays(['CellID','threshold','sigma'],library='np')
    ids=arrays['CellID'];means=arrays['threshold'];sigmas=arrays['sigma']
    if not len(ids) or len(set(ids.tolist()))!=len(ids) or not np.isfinite(means).all() or not np.isfinite(sigmas).all() or np.any(sigmas<0):
        raise RuntimeError('Invalid or duplicate threshold-map parameters')
    order=np.argsort(ids)
    digest=hashlib.sha256()
    for value in (ids[order].astype('<i8'),means[order].astype('<f8'),sigmas[order].astype('<f8')):digest.update(value.tobytes())
    return dict(path=str(path),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns,
                parameter_sha256=digest.hexdigest(),channels=len(ids),
                mean_adc=float(np.mean(means)),mean_sigma_adc=float(np.mean(sigmas)))


def add_diagnostics(here, base, plan, cuts, mip_energy, bootstrap, seed, min_events,
                    threshold_scale, threshold_file, log_dir, dry_run, plots_only, execute, unchanged):
    condition,phase,energy=plan['condition'],plan['phase'],plan['energy']+'GeV'
    record=base/'inputs'/phase/energy/condition
    shared=base/'inputs'/phase/energy/'reference'
    summary=base/'summary'/phase/energy/(condition+'.threshold.json')
    logs=here.parent/'log/beam_quantitative'/plan['label']/phase/condition/energy
    plots=base/'plots'/phase/condition/energy
    roots={s:base/'root'/s/phase/('reference' if s=='data' else condition)/energy for s in ('data','mc')}
    comparison=base/'root/comparison'/phase/condition/energy
    signature=dict(schema=SCHEMA,cuts=cuts,mip_energy_mev=mip_energy,bootstrap=bootstrap,seed=seed,min_events=min_events,
                   tolerances=plan['tolerances'],axis_layers=plan['axis_layers'],
                   source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [here/'ExtractThreshold.cc',here/'threshold.py']})
    if plots_only:
        if not summary.is_file(): raise RuntimeError('No completed threshold diagnostics to replot')
        saved=json.loads(summary.read_text())
        if not saved.get('complete'): raise RuntimeError('Threshold diagnostics are incomplete')
        if dry_run:return dict(replot=str(summary),strip_plot_binning=STRIP_PLOT_BINNING)
        rebinned={}
        for sample in ('data','mc'):
            events=roots[sample]/'pre_ssa_events.root'
            rebinned[sample]=refresh_strip_plot_bins(saved[sample+'_histograms'],events,saved['cuts'],saved['mip_energy_mev'])
        snapshot=record/'threshold_plot_sources';snapshot.mkdir(exist_ok=True)
        plot_files=[here/'PlotThreshold.cc',here/'ThresholdPlotCommon.hh',here/'PlotCommon.hh',here/'threshold.py']
        for p in plot_files:shutil.copy2(p,snapshot/p.name)
        hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in plot_files}
        (snapshot/'checksums.json').write_text(json.dumps(hashes,indent=2)+'\n')
        with (record/'run.log').open('a') as f:f.write('REPLOT_THRESHOLD_DIAGNOSTICS '+json.dumps(hashes)+'\n')
        execute([here/'PlotThreshold',saved['data_histograms'],saved['mc_histograms'],plots,
                 comparison/'threshold_comparison.root',plan['energy'],condition,
                 ','.join(cut_key(c) for c in saved['cuts'])],logs/'threshold_plot.log')
        remove_legacy_plot_names(plots)
        saved['strip_plot_binning']=STRIP_PLOT_BINNING
        summary.write_text(json.dumps(json_safe(saved),indent=2,allow_nan=False)+'\n')
        with (record/'run.log').open('a') as f:
            f.write('STRIP_PLOT_BINNING '+json.dumps(dict(binning=STRIP_PLOT_BINNING,rebinned=rebinned))+'\n')
        return saved
    if summary.exists():raise RuntimeError('Threshold diagnostics already exist; use --plots-only with --threshold-diagnostics-only')
    for pairs in plan['pairs'].values():
        if not unchanged(pairs):raise RuntimeError('Original calib/SSA/truth inputs changed; cannot augment this saved evaluation')
    metadata=threshold_map(threshold_file)
    audit=audit_digitization(plan,log_dir)
    if dry_run:return dict(condition=condition,cuts=cuts,summary=str(summary),threshold_scale=threshold_scale,threshold_map=metadata,digitization_threshold_audit=audit)
    cache=shared/'threshold_configuration.json'
    data_signature=dict(extraction_source=signature['source_hashes']['ExtractThreshold.cc'],pairs=plan['pairs']['data'],axis_layers=plan['axis_layers'])
    pre_data=roots['data']/'pre_ssa_events.root'
    if cache.exists():
        if json.loads(cache.read_text())!=data_signature or not pre_data.exists():raise RuntimeError('Shared pre-SSA cache differs or is incomplete')
    elif pre_data.exists():raise RuntimeError('Incomplete pre-SSA data extraction exists')
    else:
        execute([here/'ExtractThreshold',shared/'data.tsv',roots['data']/'events.root',pre_data,plan['axis_layers']],logs/'threshold_data.log')
        if not unchanged(plan['pairs']['data']):raise RuntimeError('Data changed during pre-SSA extraction')
        cache.write_text(json.dumps(data_signature,indent=2)+'\n')
    pre_mc=roots['mc']/'pre_ssa_events.root'
    if pre_mc.exists():raise RuntimeError('Incomplete pre-SSA MC extraction exists')
    execute([here/'ExtractThreshold',record/'mc.tsv',roots['mc']/'events.root',pre_mc,plan['axis_layers']],logs/'threshold_mc.log')
    data_hist=comparison/'pre_ssa_data_histograms.root';mc_hist=roots['mc']/'pre_ssa_histograms.root'
    result=assess_pre(pre_data,pre_mc,data_hist,mc_hist,comparison/'threshold_metrics.root',
                      summary.with_name(condition+'.threshold.metrics.tsv'),summary.with_name(condition+'.channel_occupancy.tsv'),
                      cuts,mip_energy,bootstrap,seed,min_events,plan['tolerances'],audit['status'])
    execute([here/'PlotThreshold',data_hist,mc_hist,plots,comparison/'threshold_comparison.root',plan['energy'],condition,
             ','.join(cut_key(c) for c in cuts)],logs/'threshold_plot.log')
    for pairs in plan['pairs'].values():
        if not unchanged(pairs):raise RuntimeError('Inputs changed during threshold diagnostics')
    remove_legacy_plot_names(plots)
    snapshot=record/'threshold_sources';snapshot.mkdir()
    for p in [here/'ExtractThreshold.cc',here/'threshold.py',here/'PlotThreshold.cc',here/'PlotCommon.hh',here/'ThresholdPlotCommon.hh']:shutil.copy2(p,snapshot/p.name)
    for i,log in enumerate(audit['logs']):shutil.copy2(log['path'],snapshot/('digitization_'+str(i)+'.log'))
    result.update(condition=condition,label=plan['label'],phase=phase,energy=plan['energy'],signature=signature,
                  threshold_scale=threshold_scale,threshold_map=metadata,digitization_threshold_audit=audit,
                  data_histograms=str(data_hist),mc_histograms=str(mc_hist),strip_plot_binning=STRIP_PLOT_BINNING,
                  generation_inputs=[dict(path=r['paths'][2],metadata=r['metadata'][2],entries=r['entries']) for r in plan['pairs']['mc']],
                  data_inputs=plan['pairs']['data'],material_label=plan['material_label'],ssa_selection=plan['ssa_selection'])
    summary.write_text(json.dumps(json_safe(result),indent=2,allow_nan=False)+'\n')
    with (record/'run.log').open('a') as f:f.write('THRESHOLD_DIAGNOSTICS '+json.dumps(json_safe(dict(signature=signature,threshold_scale=threshold_scale,threshold_map=metadata)))+'\n')
    return result


def unchanged(pairs):
    for row in pairs:
        for path,stored in zip(row['paths'],row['metadata']):
            if path=='-':
                if stored is not None:return False
                continue
            p=Path(path).resolve(strict=True);stat=p.stat()
            if dict(path=str(p),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns)!=stored:return False
    return True


def compare_conditions(here,base,label,phase,energy,conditions,name,dry_run,execute):
    """Compare already evaluated threshold candidates from the SAME generation."""
    import fcntl
    import re
    names=conditions.split(',')
    if not 2<=len(names)<=8 or len(set(names))!=len(names) or any(not re.fullmatch('[A-Za-z0-9_-]+',n) for n in names+[name]):
        raise RuntimeError('Use 2-8 distinct valid condition names and a valid comparison name')
    summaries=base/'summary'/phase/(energy+'GeV')
    results=[];all_rows=[]
    for condition in names:
        p=summaries/(condition+'.threshold.json')
        if not p.exists():raise RuntimeError('Run --threshold-diagnostics-only first for '+condition)
        result=json.loads(p.read_text())
        if result.get('digitization_threshold_audit',{}).get('status')=='fallback_detected':raise RuntimeError('Threshold loading failed for '+condition+'; correct digitization before a parameter-map scan')
        if not result.get('complete'):raise RuntimeError('Incomplete diagnostics: '+condition)
        if result.get('threshold_scale') is None:raise RuntimeError('Missing --mc-threshold-scale annotation: '+condition)
        if not all(r['path']!='-' and r['metadata'] for r in result['generation_inputs']):
            raise RuntimeError('Threshold comparison requires original generation ROOT/truth for every input')
        if not unchanged(result['data_inputs']):raise RuntimeError('Data used in the saved evaluation changed: '+condition)
        for r in result['generation_inputs']:
            f=Path(r['path']);meta=r['metadata'];stat=f.stat()
            if stat.st_size!=meta['bytes'] or stat.st_mtime_ns!=meta['mtime_ns']:
                raise RuntimeError('Generation input changed: '+str(f))
        source=base/'inputs'/phase/(energy+'GeV')/condition/'configuration.json'
        plan=json.loads(source.read_text())
        if not all(unchanged(rr) for rr in plan['pairs'].values()):raise RuntimeError('MC/data inputs changed: '+condition)
        if results:
            reference=results[0]
            for key in ['signature','data_inputs','generation_inputs','material_label','ssa_selection']:
                if result[key]!=reference[key]:raise RuntimeError('Threshold candidates differ in '+key+'; keep generation, data and diagnostic definitions fixed')
        results.append(result)
        with (summaries/(condition+'.threshold.metrics.tsv')).open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                all_rows.append(dict(condition=condition,threshold_scale=result['threshold_scale'],**row))
    map_status='not verified; scales are user-supplied provenance annotations'
    if all(r.get('threshold_map') for r in results):
        baseline=None
        for r in results:
            info=r['threshold_map'];p=Path(info['path']);current=threshold_map(p)
            if current!=info:raise RuntimeError('Threshold parameter file changed: '+str(p))
            with uproot.open(p) as root:a=root['expErfThre'].arrays(['CellID','threshold','sigma'],library='np')
            order=np.argsort(a['CellID']);values=(a['CellID'][order],a['threshold'][order]/r['threshold_scale'],a['sigma'][order])
            if baseline is None:baseline=values
            elif not np.array_equal(values[0],baseline[0]) or not np.allclose(values[1],baseline[1],rtol=1e-10,atol=1e-10) or not np.allclose(values[2],baseline[2],rtol=1e-10,atol=1e-10):
                raise RuntimeError('Threshold maps are not a mean-only scale scan with fixed sigma')
        map_status='parameter maps verified: common CellIDs, means scaled, sigmas fixed; producer use remains a provenance claim'
    common=None
    for result in results:
        pre=Path(result['mc_histograms']).parent/'pre_ssa_events.root'
        with uproot.open(pre) as root:valid=np.asarray(root['events']['axisValid'].array(library='np'),dtype=bool)
        if common is None:common=valid.copy()
        elif valid.shape!=common.shape:raise RuntimeError('MC sampled event counts differ')
        else:common &= valid
    if not common.any():raise RuntimeError('No common axis-valid MC events across threshold candidates')
    target=base/'summary/threshold_comparison'/phase/(energy+'GeV')/(name+'.json')
    plots=base/'plots/threshold_comparison'/phase/name/(energy+'GeV')
    roots=base/'root/threshold_comparison'/phase/name/(energy+'GeV')
    record=base/'inputs/threshold_comparison'/phase/(energy+'GeV')/name
    logs=here.parent/'log/beam_quantitative'/label/phase/'threshold_comparison'/name/(energy+'GeV')
    if record.exists() or target.exists():raise RuntimeError('Threshold comparison already exists; choose a different comparison name')
    config=dict(complete=True,label=label,phase=phase,energy=energy,comparison=name,conditions=names,
                threshold_parameter_verification=map_status,
                sources=[str(summaries/(n+'.threshold.json')) for n in names],
                interpretation='Same sampled generation events and data reference. MC comparisons use the intersection of axis-valid events across all candidates. Original acceptance is reported separately. Offline MIP-cut scans do not replace digitization.',
                accepted_events_original={r['condition']:r['mc']['accepted'] for r in results},common_mc_events=int(common.sum()))
    if dry_run:
        print(json.dumps(dict(config,plots=str(plots),summary=str(target)),indent=2));return
    lockdir=base/'inputs'/phase/(energy+'GeV')
    with (lockdir/'.run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if record.exists():raise RuntimeError('Threshold comparison initialized by another process')
        record.mkdir(parents=True);target.parent.mkdir(parents=True,exist_ok=True)
        plots.mkdir(parents=True,exist_ok=True);roots.mkdir(parents=True,exist_ok=True);logs.mkdir(parents=True,exist_ok=True)
        all_rows=[];manifest_rows=[];cohort_results={}
        with uproot.recreate(roots/'common_selection.root') as output:
            tree=output.mktree('selection',{'sample_index':'int64','accepted':'bool'})
            tree.extend({'sample_index':np.arange(len(common)),'accepted':common})
        for r in results:
            source_plan=json.loads((base/'inputs'/phase/(energy+'GeV')/r['condition']/'configuration.json').read_text())
            data_pre=base/'root/data'/phase/'reference'/(energy+'GeV')/'pre_ssa_events.root'
            mc_pre=Path(r['mc_histograms']).parent/'pre_ssa_events.root'
            table=target.with_name(name+'.'+r['condition']+'.metrics.tsv')
            occupancy=target.with_name(name+'.'+r['condition']+'.channel_occupancy.tsv')
            dh=roots/(r['condition']+'.data_histograms.root');mh=roots/(r['condition']+'.mc_histograms.root')
            spec=r['signature']
            assessed=assess_pre(data_pre,mc_pre,dh,mh,roots/(r['condition']+'.metrics.root'),table,occupancy,
                                r['cuts'],r['mip_energy_mev'],spec['bootstrap'],spec['seed'],spec['min_events'],
                                source_plan['tolerances'],r.get('digitization_threshold_audit',{}).get('status','unavailable'),common)
            cohort_results[r['condition']]=dict(data=assessed['data'],mc=assessed['mc'],totals=assessed['totals'])
            with table.open() as f:
                for row in csv.DictReader(f,delimiter='\t'):all_rows.append(dict(condition=r['condition'],threshold_scale=r['threshold_scale'],**row))
            manifest_rows.append('\t'.join([r['condition'],str(r['threshold_scale']),str(dh),str(mh),str(table)]))
        config['common_cohort_results']=cohort_results
        manifest=record/'conditions.tsv'
        manifest.write_text('\n'.join(manifest_rows)+'\n')
        with (here/'.build.lock').open('a') as build:
            fcntl.flock(build,fcntl.LOCK_EX);execute(['make','-C',here,'all'],logs/'build.log')
        execute([here/'PlotThresholdScan',manifest,plots,roots/'comparison.root',energy,name,','.join(cut_key(c) for c in results[0]['cuts'])],logs/'plot.log')
        write_rows(target.with_suffix('.metrics.tsv'),all_rows)
        with buffered_root(roots/'metrics.root') as root:
            converted=[]
            for row in all_rows:
                converted.append({k: int(v) if k in ('layer','valid_bootstraps','data_events','mc_events') else float(v) if k in ('threshold_scale','cut_mip','data_value','mc_value','difference','ci_low','ci_high','tolerance') else v for k,v in row.items()})
            root_table(root,'metrics',converted)
        snapshot=record/'sources';snapshot.mkdir()
        for p in [here/'threshold.py',here/'PlotThresholdScan.cc',here/'ThresholdPlotCommon.hh',here/'PlotCommon.hh']:shutil.copy2(p,snapshot/p.name)
        (record/'configuration.json').write_text(json.dumps(config,indent=2)+'\n')
        (record/'run.log').write_text('COMPARE_THRESHOLDS '+conditions+'\n'+map_status+'\n')
        target.write_text(json.dumps(config,indent=2)+'\n')
    print(f'Threshold comparison: {target}\nPlots: {plots}',flush=True)
