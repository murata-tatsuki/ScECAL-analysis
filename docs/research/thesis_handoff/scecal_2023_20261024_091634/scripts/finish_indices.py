"""Extend assembled indices with saved ROOT fits and exact artifact provenance."""
import pathlib,json,csv,hashlib,re,math,datetime,collections
B=pathlib.Path(__file__).absolute().parents[1]
def uid(prefix,s):return prefix+'_'+hashlib.sha256(str(s).encode()).hexdigest()[:12]
def rows(n):
 with (B/n).open() as f:return list(csv.DictReader(f))
def write(n,a,fields):
 with (B/n).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(a)
m=json.loads((B/'manifest.yaml').read_text());rr=rows('results.csv');run=rows('runs.csv');fig=rows('figures.csv');cols=list(rr[0]) if rr else 'result_id,campaign_id,sample_id,condition_id,comparison_group,data_or_mc,metric,value,uncertainty,uncertainty_kind,unit,source,source_locator,value_origin,evidence_status,status,valid,invalid_reason,events,selection,notes'.split(',')
src={x['source_id']:x for x in m['sources']};samples={x['sample_id']:x for x in m['samples']};conds={x['condition_id']:x for x in m['conditions']}
def source(path,meta=None):
 sid=uid('src',path)
 src.setdefault(sid,{'source_id':sid,'path':path,'availability':'metadata_only','sha256':None,'hash_reason':'Large ROOT not fully hashed; stored objects exported read-only'})
 if meta:
  for k in ['size','mtime','read_started_at','read_finished_at','changed_during_read','uuid','root_file_version','objects','export_path','hash_omission_reason']: 
   if k in meta:src[sid][k]=meta[k]
 return sid

def condition(path,label,kind,selection,meta=None):
 cid=uid('condition',label);sid=uid('sample',path)
 if sid not in samples:samples[sid]={'sample_id':sid,'campaign_id':'ps_sps2023','data_or_mc':kind,'particle':'e-','input_path':path,'role':'保存electron評価結果;run組成未確定','status':'暫定'}
 if cid not in conds:conds[cid]={'condition_id':cid,'label':label,'sample_ids':[sid],'source_ids':[source(path,meta)],'calibration_set_id':'unknown','generation_source_version':'unknown','generation_config':'unknown','reconstruction_id':'unknown','digitization_id':'unknown','selection_id':'unknown','selection_description':selection,'evaluation_id':'unknown','evidence_status':'確認済み','status':'暫定','limitation':'Saved-object identity and metadata only; generation provenance not inferred from current code.'}
 return cid,sid

def add(path,loc,metric,value,cid,sid,kind,unit='',err='',origin='saved_fit',valid='unknown',reason='',events='',notes=''):
 if value is None:return
 if isinstance(value,float) and not math.isfinite(value):return
 rr.append(dict(result_id=uid('result',path+'|'+loc+'|'+metric+'|'+cid),campaign_id=samples[sid]['campaign_id'],sample_id=sid,condition_id=cid,comparison_group=conds[cid]['label'],data_or_mc=kind,metric=metric,value=value,uncertainty=err if err is not None else '',uncertainty_kind='saved_fit_parameter_error; covariance/status unknown' if err not in ['',None] else 'unknown',unit=unit,source=path,source_locator=loc,value_origin=origin,evidence_status='確認済み',status='暫定',valid=valid,invalid_reason=reason,events=events,selection=conds[cid]['selection_description'],notes=notes))

for p in sorted((B/'evidence/root_exports').glob('*.json')):
 o=json.loads(p.read_text())
 if not isinstance(o,dict) or 'path' not in o:continue
 path=o['path'];source(path,o)
 if o.get('kind')!='resolution':continue
 kind='data' if '/resolution/data/' in path else 'mc'
 match=re.search(r'/([0-9.]+)GeV_([0-9.]+)mm.root$',path);label=match.group(1)+'GeV' if match else 'unknown';cog=match.group(2) if match else 'unknown'
 for stage in ['beforeEventCut','afterEventCut']:
  cid,sid=condition(path,'saved_single:'+path+':'+stage,kind,f'{stage}; CoG {cog}mm; masks/tail/CoG may already applied before beforeEventCut',o)
  conds[cid].update(saved_metadata=o.get('metadata',{}),beam_label=label,beam_label_definition='nominal filename label; momentum versus total energy not established',cog_mm=cog,root_readout_id=p.name)
  for obj,fit in o.get('saved_fits',{}).items():
   if not obj.startswith(stage+'/'):continue
   params=fit.get('parameters',[]);errs=fit.get('errors',[]);mu=params[1] if len(params)>1 else None;sig=params[2] if len(params)>2 else None
   bad=mu is not None and (mu<=0 or sig is None or sig<=0)
   excluded=kind=='data' and label=='200GeV'
   valid='false' if bad or excluded else 'unknown';reason='nonpositive saved Gaussian mean or sigma' if bad else 'data200GeV explicitly excluded from resolution value/parameter fits; retain saved distribution fit as provenance' if excluded else ''
   for i,v in enumerate(params):
    metric=('gaussian_'+['amplitude','mean','sigma'][i]) if i<3 else 'fit_parameter_'+str(i)
    unit=('hits' if 'nhit' in obj else 'MeV') if i in [1,2] else ''
    add(path,obj+f'/parameter[{i}]',metric,v,cid,sid,kind,unit,errs[i] if i<len(errs) else '',valid=valid,reason=reason,notes='Saved TF1, no refit. TF1 xmin/xmax may differ from final explicit fit range. This row is not sigma/mean calculated anew.')
   for k in ['chi2','ndf','xmin','xmax']:add(path,obj+'/'+k,'fit_'+k,fit.get(k),cid,sid,kind,'MeV' if k in ['xmin','xmax'] and 'nhit' not in obj else '',origin='saved_fit',valid=valid,reason=reason)
  for obj,h in o.get('histograms',{}).items():
   if not obj.startswith(stage+'/'):continue
   for k in ['entries','nbins','xmin','xmax']:
    add(path,obj+'/'+k,'histogram_'+k,h.get(k),cid,sid,kind,'entries' if k=='entries' else 'bins' if k=='nbins' else 'hits' if 'nhits' in obj else 'MeV',origin='saved_table',events=h.get('entries',''),notes='Histogram saved counters/axis metadata, not raw acquisition or complete cut flow.')
 if o.get('kind')=='constants':src[source(path,o)]['tables']=o.get('tables',[])


# Canonical saved resolution graphs: first full-series panel; other views remain in JSON.
cols += ['x_value','x_unit','x_definition']
p=B/'evidence/comparison_canvas.json'
if p.is_file():
 o=json.loads(p.read_text());path=o['source'];source(path,o)
 source_record=src[uid('src',path)];source_record['captured_at']=o.get('captured_at');source_record['mtime']=datetime.datetime.fromtimestamp(o['mtime_ns']/1e9,datetime.timezone(datetime.timedelta(hours=9))).isoformat()
 source_record['root_version']=o.get('root_version');source_record['archive_export']='evidence/comparison_canvas.json'
 for obj in o.get('objects',[]):
  loc=obj['path'];is_plain=loc.startswith('resolution/resolution_plots_1[0]/');is_noise=loc.startswith('resolution_with_noise/resolution_plots_noise_1[0]/')
  if not obj.get('points') or not (is_plain or is_noise):continue
  funs=obj.get('functions',[]);match=next((re.search(r'res_func_low_(?:noise_)?(\d+)$',f.get('name','')) for f in funs if re.search(r'res_func_low_(?:noise_)?(\d+)$',f.get('name',''))),None)
  if not match:continue
  ds=int(match.group(1));kind='data' if ds%2==0 else 'mc';cog=[200,20,10,5][ds//2]
  cid,sid=condition(path,f'comparison_oct12_series_{ds}',kind,f'beforeEventCut; CoG {cog}mm; saved comparison series {ds}; see log and metadata',o)
  conds[cid].update(cog_mm=cog,series_index=ds,series_mapping_source='Saved attached function index, first-panel legend order, and S_RESLOG line9',evaluation_id='resolution_current',evaluation_relation='saved estimator compatible with captured implementation; generation source version unknown')
  if not is_noise:
   for point in obj['points']:
    val=point['y'];bad=val<=0
    add(path,loc+f"/point[{point['index']}]",'gaussian_sigma_over_mu',val,cid,sid,kind,'dimensionless',point.get('ey',''),origin='saved_table',valid='false' if bad else 'unknown',reason='nonpositive saved resolution' if bad else '',notes='Saved TGraphErrors point. Linear/log and duplicate panels are one measurement, not independent points. Saved graph error propagation; full covariance/systematics and fit stability not validated.')
    rr[-1].update(x_value=point['x'],x_unit='GeV',x_definition='nominal filename energy-axis label; momentum/total energy interpretation unverified')
  for f in funs:
   model='with_noise' if is_noise else 'without_noise';conds[cid].setdefault('saved_curve_models',{})[model]={'formula_or_title':f.get('title'),'function_name':f['name'],'chi2':f.get('chi2'),'ndf':f.get('ndf'),'saved_range':[f.get('xmin'),f.get('xmax')],'fit_status':f.get('fit_status')}
   for par in f['parameters']:
    names=['stochastic_a','noise_b','constant_c'] if is_noise else ['stochastic_a','constant_c'];units=['sqrt(GeV)','GeV','dimensionless'] if is_noise else ['sqrt(GeV)','dimensionless']
    i=par['index'];add(path,loc+'/'+f['name']+f'/parameter[{i}]',model+'_'+names[i],par['value'],cid,sid,kind,units[i],par.get('error',''),notes='Stored TF1 parameter; raw sign retained because squared model can have sign degeneracy. Legend uses absolute values. Fit status/covariance and physical adoption unknown; function saved range need not be actual fit interval.')

# Fill table condition evidence links without treating proximity as proof of use.
for c in conds.values():
 for id in c.get('source_ids',[]):
  s=src.get(id,{})
  path=s.get('path','');parent=str(pathlib.Path(path).parent)
  nearby=[x['source_id'] for x in src.values() if x.get('path','').startswith(parent+'/inputs/')]
  if nearby:c['nearby_generation_evidence_ids']=nearby;c['nearby_evidence_relation']='adjacent input bundle; consult each manifest before assigning generation source'
# Verify known invalid saved table fits without deleting the transcript.
bycond=collections.defaultdict(list)
for r in rr:bycond[r['condition_id']].append(r)
for cid,a in bycond.items():
 for r in a:
  if r['metric']=='e_fit':
   try:bad=float(r['value'])<=0
   except Exception:bad=False
   if bad:
    for z in a:z['valid']='false';z['invalid_reason']='saved Gaussian energy mean is nonpositive'
# Every condition includes the sample rows indexed under it.
for r in run:
 if r['condition_id'] in conds and r['sample_id'] not in conds[r['condition_id']]['sample_ids']:conds[r['condition_id']]['sample_ids'].append(r['sample_id'])
# Normalize capture hash/time fields and root table metadata.
for s in src.values():
 if s.get('copy_sha256'):s['archive_sha256']=s['copy_sha256']
 if s.get('capture_extent')=='whole':s['extent']='whole_file'
 if s.get('archive_path') and not s.get('archive_sha256'):
  pp=B/s['archive_path']
  if pp.is_file():s['archive_sha256']=hashlib.sha256(pp.read_bytes()).hexdigest()
 if s.get('size') is None:s['size']=None
 if s.get('mtime') is None:s['mtime']=None
 if not s.get('captured_at'):s['captured_at']=s.get('read_started_at',s.get('retrieved_at','unknown'))
 if not s.get('sha256'):s.setdefault('hash_reason',s.get('hash_omission_reason','whole source not captured; directory index or metadata only'))
for p in sorted((B/'evidence/root_exports').glob('*.json')):
 o=json.loads(p.read_text())
 if isinstance(o,dict) and o.get('kind')=='constants':src[source(o['path'],o)]['tables']=o.get('tables',[])
# Manifest parameter and claim records explicitly refer to captured current implementation.
m['claims']=[{'claim_id':'claim_'+key.lower(),'source_id':m['source_aliases'][key],'source_locator':loc,'scope':'captured current implementation, not automatically past generation source','evidence_status':'確認済み','status':'暫定'} for key,loc in [('S_CAL','Calibration.cxx: MIPFactorExtract_simulation L401–463; InterFactorExtract_simulation L465–530; Calibration L659–709'),('S_MAP','EBUdecode.cxx L17–50; tempReconstruction L74–105'),('S_DIGI','MCDigi.cxx PixelAva L141–170, PixelCorrection L175–195, StepDigi L199–302'),('S_COMB','CombFit.cc common-spacing model and fit extraction'),('S_GEOM','DetectorConstruction.cc L138–221,346–358'),('S_G4MAIN','CEPCScECAL.cc L66–86'),('S_SINGLE','SingleEnergyAnalysis.cc saved histogram/fit definitions'),('S_MULTI','MultiEnergyAnalysis.cc resolution estimator and fit point exclusions')]]
vals=[('MIP_energy',.305,'MeV','all'),('PDE',.1,'probability','10um'),('PDE',.25,'probability','15um'),('pCT',.1,'neighbor_cascade_probability','10um'),('pCT',.15,'neighbor_cascade_probability','15um'),('pixel_count',10000,'pixels','10um'),('pixel_count',4489,'pixels','15um'),('uniform_light_half_width',.042,'relative_half_width','all'),('ADC_sigma_coefficient',3,'ADC/sqrt(fired_pixels)','all'),('fallback_gain',8,'ADC/pe','10um'),('fallback_gain',22.5,'ADC/pe','15um')]
m['digitization_parameters']=[dict(parameter=k,value=v,unit=u,group=g,source_id=m['source_aliases']['S_DIGI'],source_locator='MCDigi.cxx StepDigi L199–295',uncertainty=None,origin='current code; calibration/manufacturer/tuning history not fully established',status='暫定') for k,v,u,g in vals]
m['capture_resume']={'interrupted_before_completion':True,'resumed_on':'2026-10-26 Asia/Tokyo','reason':'Agent restarted; prior NFS waits prevented completion. Partial files inspected; current source recaptured with file-level times. No atomic snapshot claim.'}
m['sources']=list(src.values());m['samples']=list(samples.values());m['conditions']=list(conds.values());m['figures']=fig
m['completed_at']=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
(B/'manifest.yaml').write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n')
# Each row has a unique identity. Repeated object primitives remain distinct by path.
unique={r['result_id']:r for r in rr};rr=list(unique.values())
write('results.csv',rr,cols)
(B/'audit/final_index_counts.json').write_text(json.dumps({'results':len(rr),'runs':len(run),'figures':len(fig),'sources':len(src),'conditions':len(conds),'completed_at':m['completed_at']},indent=2))
print('extended',len(rr),'saved metrics',len(run),'run usages',len(fig),'figures',flush=True)
