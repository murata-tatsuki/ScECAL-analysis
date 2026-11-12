"""Paired physical-strip diagnostics for truth -> digitization -> calibration."""
import csv
import json
import io
import tempfile
from pathlib import Path
import awkward as ak
import numpy as np
import uproot
from evaluate import buffered_root, json_safe as _json_safe, resamples

LAYERS = np.arange(31, dtype=float) - .5
FINE = (np.arange(43, dtype=float)-21.)*5.3
ALIGNED = (np.arange(86, dtype=float)-42.5)*5.3
MIP_BINS = np.r_[-2., -1., 0., np.geomspace(.001, 2000., 147)]
TRUTH_BINS = np.geomspace(.001, 2000., 81)
ADC_BINS = np.arange(-.5, 5019.5, 20.)
# SiPM cell size, supplied by the detector layer assignment. Layers 30/31
# are classified here, but the existing analysis selects layers 0..29 only.
SIPM_LAYERS = {'10um': tuple(range(4,28)), '15um': (0,1,2,3,28,29,30,31)}
SIPM_SURVIVAL_GROUPS = {
    '10um': {'all': tuple(range(4,28)), 'layers_4_9': tuple(range(4,10)),
             'layers_10_19': tuple(range(10,20)), 'layers_20_27': tuple(range(20,28))},
    '15um': {'all': SIPM_LAYERS['15um'], 'layers_0_3': (0,1,2,3),
             'layers_28_29': (28,29), 'layer_29': (29,)}}
SIPM_SCHEMA = 1
GROUPS = {'all': (0, 30), 'layers_0_3': (0, 4), 'layers_4_9': (4, 10),
          'layers_10_19': (10, 20), 'layers_20_28': (20, 29), 'layer_29': (29, 30)}


def json_safe(value):
    if isinstance(value,np.ndarray): return json_safe(value.tolist())
    if isinstance(value,dict): return {k:json_safe(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [json_safe(v) for v in value]
    return _json_safe(value)


def load_mc(path, channel_ids=None):
    with uproot.open(path) as f: a = f['events'].arrays(library='ak')
    valid = np.asarray(a.axisValid, dtype=bool)
    audit = dict(sampled=len(valid), accepted=int(valid.sum()), axis_invalid=int((~valid).sum()))
    for name in ('invalidTruth', 'invalidRaw', 'invalidCalib'):
        audit[name] = int(ak.sum(a[name]))
    a = a[valid]; n = len(a)
    if not n: raise RuntimeError('No accepted MC events')
    s = {name: np.asarray(ak.flatten(a[name])) for name in
         ('cellID', 'layer', 'x', 'y', 'truthMultiplicity', 'rawMultiplicity', 'taggedMultiplicity',
          'calibMultiplicity', 'truthEnergy', 'hgADC', 'lgADC', 'calibEnergy')}
    s.update(n=n, event=np.repeat(np.arange(n), np.asarray(ak.num(a.layer))),
             axisX=np.asarray(a.axisX), axisY=np.asarray(a.axisY), audit=audit)
    apply_channel_mask(s,channel_ids)
    s['truth'] = (s['truthMultiplicity']>0) & (s['truthEnergy']>0)
    s['digitized'] = s['rawMultiplicity']>0; s['calibrated'] = s['calibMultiplicity']>0
    for name in ('truth', 'raw', 'calib'):
        v=s[name+'Multiplicity']; audit[name+'_duplicate_rows'] = int(np.maximum(v-1, 0).sum())
    audit['truth_nonpositive_strips'] = int(((s['truthMultiplicity']>0) & (s['truthEnergy']<=0)).sum())
    audit['raw_untagged_strips'] = int((s['digitized'] & (s['taggedMultiplicity']==0)).sum())
    audit['calib_nonpositive_strips'] = int((s['calibrated'] & (s['calibEnergy']<=0)).sum())
    audit['digitized_without_truth'] = int((s['digitized'] & (s['truthMultiplicity']==0)).sum())
    audit['calibrated_without_digitized'] = int((s['calibrated'] & ~s['digitized']).sum())
    audit['calibrated_without_truth'] = int((s['calibrated'] & (s['truthMultiplicity']==0)).sum())
    audit['calibrated_without_tagged_raw'] = int((s['calibrated'] & (s['taggedMultiplicity']==0)).sum())
    return s


def load_data(path,beam_path,channel_ids=None):
    with uproot.open(path) as f: a = f['events'].arrays(library='ak')
    valid=np.asarray(a.axisValid,dtype=bool)
    audit=dict(sampled=len(valid),accepted=int(valid.sum()),axis_invalid=int((~valid).sum()),
               invalid_calib=int(ak.sum(a.invalidPre)),merged_rows=int(ak.sum(a.mergedEntries)))
    with uproot.open(beam_path) as f: axes=f['events'].arrays(['axisValid','axisX','axisY','fileIndex','entry'],library='np')
    if not np.array_equal(valid,axes['axisValid']) or not np.array_equal(np.asarray(a.fileIndex),axes['fileIndex']) or not np.array_equal(np.asarray(a.entry),axes['entry']):
        raise RuntimeError('Data axis/strip cache identity mismatch')
    a=a[valid]; n=len(a)
    if not n: raise RuntimeError('No accepted data events')
    s={name:np.asarray(ak.flatten(a[name])) for name in ('cellID','layer','x','y','energy')}
    s.update(n=n,event=np.repeat(np.arange(n),np.asarray(ak.num(a.layer))),audit=audit,
             axisX=axes['axisX'][valid],axisY=axes['axisY'][valid])
    apply_channel_mask(s,channel_ids)
    return s


def apply_channel_mask(sample,channel_ids):
    # Preserve the complete original event denominator and entrance axes.
    # Channel arrays have one element per physical strip; event indices are
    # never renumbered after the mask, so zero-hit events remain in the sample.
    if channel_ids is None:return
    keep=np.isin(sample['cellID'],np.asarray(channel_ids,dtype=np.int32))
    before=len(keep)
    for key in list(sample):
        if key not in ('n','audit','axisX','axisY'):
            sample[key]=sample[key][keep]
    sample['audit']['channel_mask']=dict(input_strips=before,retained_strips=int(keep.sum()),
        events_without_any_selected_strip=int((np.bincount(sample['event'],minlength=sample['n'])==0).sum()))


def matrix(s, mask, weights=None, bins=30, indices=None):
    index=s['layer'] if indices is None else indices
    v=np.ones(mask.sum()) if weights is None else np.asarray(weights)[mask]
    return np.bincount(s['event'][mask]*bins+index[mask],weights=v,minlength=s['n']*bins).reshape(s['n'],bins)


def interval(v):
    low=np.full(v.shape[1],np.nan); high=low.copy()
    for i in range(v.shape[1]):
        x=v[1:,i][np.isfinite(v[1:,i])]
        if len(x): low[i],high[i]=np.quantile(x,[.025,.975])
    return dict(mean=v[0],low=low,high=high)


def divide(a,b):
    out=np.full(np.broadcast_shapes(a.shape,b.shape),np.nan)
    return np.divide(a,b,out=out,where=b!=0)


def sipm_histograms(mc,data,wm,mip):
    """Partition pooled diagnostics by SiPM size without changing the event cohort."""
    hist={}; report=dict(schema=SIPM_SCHEMA,classification=SIPM_LAYERS,
        analyzed_layers={k:[l for l in v if l<30] for k,v in SIPM_LAYERS.items()},
        excluded_layers=[30,31],flows={},truth_energy_survival={})
    for sample in (mc,data):
        if np.any(~np.isin(sample['layer'],np.arange(30))):
            raise RuntimeError('SiPM plots expect the existing layer 0..29 analysis')
    index=np.searchsorted(TRUTH_BINS,mc['truthEnergy']/mip,side='right')-1
    in_range=(index>=0)&(index<len(TRUTH_BINS)-1)
    corr=np.linspace(0,50,101)
    for kind,layers in SIPM_LAYERS.items():
        mm=np.isin(mc['layer'],layers); dm=np.isin(data['layer'],layers)
        flows={}; survival={}
        def spectrum(key,values,edges):
            hist['sipm/'+kind+'/signal/'+key]=(np.histogram(values,bins=edges)[0].astype(float),edges)
            flows[key]=dict(entries=len(values),underflow=int((values<edges[0]).sum()),
                            overflow=int((values>=edges[-1]).sum()))
        spectrum('truth',mc['truthEnergy'][mm & mc['truth']]/mip,MIP_BINS)
        spectrum('truth_digitized',mc['truthEnergy'][mm & mc['truth'] & mc['digitized']]/mip,MIP_BINS)
        spectrum('truth_calibrated',mc['truthEnergy'][mm & mc['truth'] & mc['calibrated']]/mip,MIP_BINS)
        spectrum('calibrated',mc['calibEnergy'][mm & mc['calibrated']]/mip,MIP_BINS)
        spectrum('data',data['energy'][dm]/mip,MIP_BINS)
        single=mm & (mc['rawMultiplicity']==1)
        spectrum('hg_adc',mc['hgADC'][single],ADC_BINS)
        spectrum('lg_adc',mc['lgADC'][single],ADC_BINS)
        matched=mm & mc['truth'] & mc['calibrated']
        values=np.histogram2d(mc['truthEnergy'][matched]/mip,mc['calibEnergy'][matched]/mip,bins=(corr,corr))[0]
        hist['sipm/'+kind+'/response/truth_calib']=(values,corr,corr)
        flows['response']=dict(matched_strips=int(matched.sum()),in_display=int(values.sum()))
        flows['survival_truth']=dict(underflow=int((mm & mc['truth'] & (index<0)).sum()),
                                    overflow=int((mm & mc['truth'] & (index>=len(TRUTH_BINS)-1)).sum()))
        for group,group_layers in SIPM_SURVIVAL_GROUPS[kind].items():
            mask=mc['truth'] & in_range & np.isin(mc['layer'],group_layers)
            denominator=matrix(mc,mask,bins=len(TRUTH_BINS)-1,indices=index)
            den=wm@denominator; support=denominator.sum(axis=0)
            survival[group]=dict(layers=[l for l in group_layers if l<30],
                edges_mip=TRUTH_BINS.tolist(),truth_strip_entries=support.tolist())
            prefix='sipm/'+kind+'/survival/'+group+'/'
            for stage in ('digitized','calibrated'):
                ci=interval(divide(wm@matrix(mc,mask & mc[stage],bins=len(TRUTH_BINS)-1,indices=index),den))
                survival[group][stage]=ci
                for q,a in ci.items():hist[prefix+stage+'_'+q]=(a,TRUTH_BINS)
            hist[prefix+'support']=(support.astype(float),TRUTH_BINS)
        report['flows'][kind]=flows; report['truth_energy_survival'][kind]=survival
    return hist,json_safe(report)


def has_sipm_histograms(path):
    with uproot.open(path) as f:
        if 'metadata/sipm' not in f:return False
        report=json.loads(str(f['metadata/sipm']))
        return (report.get('schema')==SIPM_SCHEMA and
                report.get('classification')==json_safe(SIPM_LAYERS))


def ensure_sipm_histograms(hist_path,mc_path,data_path,data_beam_path,report,validate=None):
    """Add only missing display histograms to a completed result, atomically.

    Use the saved channel mask, bootstrap settings and matched event caches.
    All original ROOT objects, layer metrics and depth-group statistics survive.
    """
    if has_sipm_histograms(hist_path):
        with uproot.open(hist_path) as f:return json.loads(str(f['metadata/sipm']))
    selection=report.get('channel_selection')
    ids=None if selection is None or selection['criteria'] is None else selection['selected_ids']
    mc=load_mc(mc_path,ids); data=load_data(data_path,data_beam_path,ids)
    if (mc['n']!=report['audit']['mc']['accepted'] or data['n']!=report['audit']['data']['accepted']):
        raise RuntimeError('Saved event cohort differs; rerun the analysis before replotting')
    rng=np.random.default_rng(report['seed'])
    wm=resamples(mc['n'],report['bootstrap'],rng)/mc['n']
    hist,sipm=sipm_histograms(mc,data,wm,report['mip_energy_mev'])
    # ROOT performs many small metadata writes: update in RAM rather than
    # issuing each write on the shared filesystem, then replace atomically.
    class Buffer(io.BytesIO):
        def close(self):pass
    buffer=Buffer(hist_path.read_bytes())
    try:
        with uproot.update(buffer) as f:
            for key,v in hist.items():f[key]=v
            f['metadata/sipm']=json.dumps(sipm,allow_nan=False)
        payload=buffer.getvalue()
    finally:io.BytesIO.close(buffer)
    scratch=Path(__file__).resolve().parents[2]/'.scratch/strip_stage_sipm'
    scratch.mkdir(parents=True,exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(dir=scratch) as tmp:
            target=Path(tmp)/'histograms.root'
            target.write_bytes(payload)
            if validate is not None:validate()
            target.replace(hist_path)
    finally:
        try:scratch.rmdir()
        except OSError:pass
    return sipm


def evaluate(mc_path,data_path,data_beam_path,hist_path,summary_path,table_path,bootstrap,seed,mip,channel_selection=None):
    channel_ids=None if channel_selection is None or channel_selection['criteria'] is None else channel_selection['selected_ids']
    mc=load_mc(mc_path,channel_ids); data=load_data(data_path,data_beam_path,channel_ids)
    rng=np.random.default_rng(seed)
    wm=resamples(mc['n'],bootstrap,rng)/mc['n']; wd=resamples(data['n'],bootstrap,rng)/data['n']
    masks={k:mc[k] for k in ('truth','digitized','calibrated')}
    for cut in (.5,1.,2.): masks['calib_mip_'+str(cut).replace('.','_')]=mc['calibrated'] & (mc['calibEnergy']>=cut*mip)
    masks.update(truth_kept_digitized=mc['truth'] & mc['digitized'], truth_kept_calibrated=mc['truth'] & mc['calibrated'],
                 lost_digitization=mc['truth'] & ~mc['digitized'], lost_calibration=mc['digitized'] & ~mc['calibrated'],
                 added_digitization=mc['digitized'] & (mc['truthMultiplicity']==0),
                 added_calibration=mc['calibrated'] & ~mc['digitized'],
                 raw_untagged=mc['digitized'] & (mc['taggedMultiplicity']==0),
                 calib_nonpositive=mc['calibrated'] & (mc['calibEnergy']<=0))
    event_values={key:matrix(mc,mask) for key,mask in masks.items()}
    data_values={'data':matrix(data,np.ones(len(data['layer']),dtype=bool))}
    for cut in (.5,1.,2.): data_values['data_mip_'+str(cut).replace('.','_')]=matrix(data,data['energy']>=cut*mip)
    for key,mask in [('truth_energy',mc['truth']),('truth_energy_digitized',mc['truth'] & mc['digitized']),
                     ('truth_energy_calibrated',mc['truth'] & mc['calibrated'])]:
        event_values[key]=matrix(mc,mask,mc['truthEnergy'])
    event_values['calib_energy']=matrix(mc,mc['calibrated'],np.maximum(mc['calibEnergy'],0))
    event_values['matched_calib_energy']=matrix(mc,mc['truth'] & mc['calibrated'],mc['calibEnergy'])
    data_values['data_energy']=matrix(data,np.ones(len(data['layer']),dtype=bool),np.maximum(data['energy'],0))
    draws={k:wm@v for k,v in event_values.items()}; draws.update({k:wd@v for k,v in data_values.items()})
    draws['retention_digitization']=divide(draws['truth_kept_digitized'],draws['truth'])
    draws['retention_calibration']=divide(draws['digitized']-draws['lost_calibration'],draws['digitized'])
    draws['retention_total']=divide(draws['truth_kept_calibrated'],draws['truth'])
    draws['energy_retention_digitization']=divide(draws['truth_energy_digitized'],draws['truth_energy'])
    draws['energy_retention_total']=divide(draws['truth_energy_calibrated'],draws['truth_energy'])
    draws['calibration_response']=divide(draws['matched_calib_energy'],draws['truth_energy_calibrated'])
    draws['mc_data_difference']=divide(draws['calibrated'],draws['data'])-1.
    for cut in ('0_5','1_0','2_0'):
        draws['mc_data_difference_mip_'+cut]=divide(draws['calib_mip_'+cut],draws['data_mip_'+cut])-1.
    for key in ('truth','digitized','calibrated','data'):
        s=data if key=='data' else mc; mask=np.ones(len(s['layer']),dtype=bool) if key=='data' else masks[key]
        w=wd if key=='data' else wm
        u=np.where(s['layer']%2,s['x'],s['y'])-np.where(s['layer']%2,s['axisX'][s['event']],s['axisY'][s['event']])
        n=draws[key]; mu=divide(w@matrix(s,mask,u),n)
        draws['width_'+key]=np.sqrt(np.maximum(0,divide(w@matrix(s,mask,u*u),n)-mu*mu))
    curves={k:interval(v) for k,v in draws.items()}
    groups={}
    for group,(lo,hi) in GROUPS.items():
        g={k:interval(v[:,lo:hi].sum(axis=1)[:,None]) for k,v in draws.items() if k in event_values or k in data_values}
        for name,num,den in [('retention_digitization','truth_kept_digitized','truth'),
                             ('retention_calibration',None,'digitized'),('retention_total','truth_kept_calibrated','truth'),
                             ('calibration_response','matched_calib_energy','truth_energy_calibrated'),
                             ('mc_data_difference','calibrated','data')]:
            top=(draws['digitized']-draws['lost_calibration']) if num is None else draws[num]
            v=divide(top[:,lo:hi].sum(axis=1)[:,None],draws[den][:,lo:hi].sum(axis=1)[:,None])
            if name=='mc_data_difference':v-=1
            g[name]=interval(v)
        groups[group]={k:{q:float(a[0]) for q,a in v.items()} for k,v in g.items()}
    audit={'mc':mc['audit'],'data':data['audit']}; hist={}; flows={}
    def spectrum(key,values,edges):
        counts,edge=np.histogram(values,bins=edges); hist[key]=(counts.astype(float),edge)
        return dict(entries=len(values),underflow=int((values<edges[0]).sum()),overflow=int((values>=edges[-1]).sum()))
    flows['truth']=spectrum('signal/truth',mc['truthEnergy'][mc['truth']]/mip,MIP_BINS)
    flows['truth_digitized']=spectrum('signal/truth_digitized',mc['truthEnergy'][mc['truth'] & mc['digitized']]/mip,MIP_BINS)
    flows['truth_calibrated']=spectrum('signal/truth_calibrated',mc['truthEnergy'][mc['truth'] & mc['calibrated']]/mip,MIP_BINS)
    flows['calibrated']=spectrum('signal/calibrated',mc['calibEnergy'][mc['calibrated']]/mip,MIP_BINS)
    flows['data']=spectrum('signal/data',data['energy']/mip,MIP_BINS)
    single=mc['rawMultiplicity']==1
    flows['adc_hg']=spectrum('signal/hg_adc',mc['hgADC'][single],ADC_BINS)
    flows['adc_lg']=spectrum('signal/lg_adc',mc['lgADC'][single],ADC_BINS)
    truth_index=np.searchsorted(TRUTH_BINS,mc['truthEnergy']/mip,side='right')-1
    in_range=(truth_index>=0)&(truth_index<len(TRUTH_BINS)-1); survival_bins={}
    for group,(lo,hi) in GROUPS.items():
        mask=mc['truth'] & in_range & (mc['layer']>=lo)&(mc['layer']<hi)
        denominator=matrix(mc,mask,bins=len(TRUTH_BINS)-1,indices=truth_index)
        den=wm@denominator; support=denominator.sum(axis=0)
        survival_bins[group]={'edges_mip':TRUTH_BINS.tolist(),'truth_strip_entries':support.tolist()}
        for stage in ('digitized','calibrated'):
            v=divide(wm@matrix(mc,mask & mc[stage],bins=len(TRUTH_BINS)-1,indices=truth_index),den)
            ci=interval(v); survival_bins[group][stage]=ci
            for q,a in ci.items():hist['survival/'+group+'/'+stage+'_'+q]=(a,TRUTH_BINS)
        hist['survival/'+group+'/support']=(support.astype(float),TRUTH_BINS)
    flows['survival_truth']=dict(underflow=int((mc['truth'] & (truth_index<0)).sum()),
                                 overflow=int((mc['truth'] & (truth_index>=len(TRUTH_BINS)-1)).sum()))
    profile_sets={k:(mc,mc[k]) for k in ('truth','digitized','calibrated')}
    profile_sets['data']=(data,np.ones(len(data['layer']),dtype=bool))
    for name,(s,selected) in profile_sets.items():
        u=np.where(s['layer']%2,s['x'],s['y'])
        offset=np.where(s['layer']%2,s['axisX'][s['event']],s['axisY'][s['event']])
        for l in range(30):
            mask=selected & (s['layer']==l)
            for frame,coord,edge in [('lab',u,FINE),('aligned',u-offset,ALIGNED)]:
                if np.any((coord[mask]<edge[0])|(coord[mask]>=edge[-1])):raise RuntimeError('Profile outside plotting range')
                hist['profiles/'+frame+'/'+name+'/layer_'+str(l)]=(np.histogram(coord[mask],bins=edge)[0]/s['n'],edge)
    matched=mc['truth'] & mc['calibrated']; corr=np.linspace(0,50,101)
    hist['response/truth_calib']=(np.histogram2d(mc['truthEnergy'][matched]/mip,mc['calibEnergy'][matched]/mip,bins=(corr,corr))[0],corr,corr)
    flows['response']=dict(matched_strips=int(matched.sum()),in_display=int(hist['response/truth_calib'][0].sum()))
    for key,v in curves.items():
        for q,a in v.items():hist['curves/'+key+'_'+q]=(a,LAYERS)
    sipm_hist,sipm=sipm_histograms(mc,data,wm,mip)
    hist.update(sipm_hist)
    with buffered_root(hist_path) as f:
        for key,v in hist.items():f[key]=v
        f['metadata/audit']=json.dumps(audit);f['metadata/flows']=json.dumps(flows)
        f['metadata/mip_energy_mev']=str(mip)
        f['metadata/sipm']=json.dumps(sipm,allow_nan=False)
        if channel_selection is not None:
            f['metadata/threshold_selection']=json.dumps(channel_selection,allow_nan=False)
            f['quality/selected_counts']=(np.asarray(channel_selection['selected_by_layer'],dtype=float),LAYERS)
            f['quality/available_counts']=(np.full(30,210.),LAYERS)
            tree=f.mktree('quality/selected_channels',{'CellID':'int32'})
            tree.extend({'CellID':np.asarray(channel_selection['selected_ids'],dtype=np.int32)})
    report=dict(schema=1,complete=False,bootstrap=bootstrap,seed=seed,mip_energy_mev=mip,
                audit=audit,layer_metrics=curves,groups=groups,truth_energy_survival=survival_bins,flows=flows,
                definitions=dict(truth='Unique connected physical strips with summed Edep > 0',
                  digitized='Unique connected strips present in Raw_Hit after producer threshold',
                  calibrated='Unique connected strips present in Calib_Hit; signed energies merged',
                  retention='Ratio of total matched strip counts; paired whole-event bootstrap',
                  response='Signed reconstructed energy sum / truth Edep sum on matched truth-positive calibrated strips',
                  width='Fine-axis RMS about its pooled mean after the fixed entrance-centroid position subtraction'),
                limitations=['No pre-threshold ADC is stored; rejected ADC fluctuations and threshold draws cannot be recovered.',
                  'Data has no truth stage. A larger MC truth count alone does not identify the cause of a data/MC difference.',
                  'Truth MIP is deposited MeV / nominal MIP energy, not an observed ADC or calibrated MIP.',
                  'Raw ADC contains pedestal. ADC plots use single-row strips; no implicit MIP conversion.',
                  'Intervals cover event sampling, not parameter or file-to-file systematic uncertainty.'])
    report['sipm']=sipm
    if channel_selection is not None:report['channel_selection']=channel_selection
    summary_path.write_text(json.dumps(json_safe(report),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    with table_path.open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['scope','metric','mean','low_95','high_95'])
        for key,v in curves.items():
            for l in range(30):w.writerow([f'layer_{l}',key,*[v[q][l] for q in ('mean','low','high')]])
        for group,g in groups.items():
            for key,v in g.items():w.writerow([group,key,*[v[q] for q in ('mean','low','high')]])
    return json_safe(report)
