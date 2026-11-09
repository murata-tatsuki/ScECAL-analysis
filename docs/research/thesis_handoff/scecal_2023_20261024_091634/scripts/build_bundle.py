#!/usr/bin/env python3
"""Assemble portable indices from copied evidence and read-only exports; no event calculations."""
import pathlib,json,csv,hashlib,re,os,math,datetime,collections,concurrent.futures
B=pathlib.Path(__file__).absolute().parents[1]
JST=datetime.timezone(datetime.timedelta(hours=9)); now=lambda:datetime.datetime.now(JST).isoformat()
R='/home/murata_t/ScECAL_BeamTest'; A='/home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO'; S='/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new'
def uid(prefix,s):return prefix+'_'+hashlib.sha256(str(s).encode()).hexdigest()[:12]
def readj(p,default=None):
 try:return json.loads(p.read_text())
 except Exception:return default

OUTPUT_POOL=concurrent.futures.ThreadPoolExecutor(max_workers=6)
OUTPUT_JOBS=[]
def _dumpcsv(name,rows,cols):
 with (B/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)
def dumpcsv(name,rows,cols):
 OUTPUT_JOBS.append((name,OUTPUT_POOL.submit(_dumpcsv,name,rows,cols)))
def save_text(name,contents):
 OUTPUT_JOBS.append((name,OUTPUT_POOL.submit((B/name).write_text,contents)))

def clean(x):
 if isinstance(x,float) and not math.isfinite(x):return None
 if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
 if isinstance(x,list):return [clean(v) for v in x]
 return x

# Normalize both capture records and each actual bundled entity. Files are always real copies.
sources={};captured=[]
for p in sorted((B/'audit').glob('*capture.json')):
 a=readj(p,[])
 if isinstance(a,dict):a=a.get('files',a.get('sources',a.get('records',[])))
 for x in a:
  if isinstance(x,dict):captured.append(x)
for x in captured:
 orig=x.get('source_path',x.get('path',x.get('source','')))
 if not orig:continue
 z=dict(x);z['source_id']=uid('src',orig);z['path']=orig
 if z.get('copy_sha256'):z['archive_sha256']=z['copy_sha256']
 if z.get('capture_extent')=='whole':z['extent']='whole_file'
 sources[orig]=z
for p in sorted((B/'evidence/files').rglob('*')) if (B/'evidence/files').exists() else []:
 if not p.is_file():continue
 orig='/'+str(p.relative_to(B/'evidence/files'))
 if orig in sources and sources[orig].get('availability')=='copied':continue
 z=sources.setdefault(orig,{'source_id':uid('src',orig),'path':orig})
 z.update(archive_path=str(p.relative_to(B)),archive_size=p.stat().st_size,archive_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),availability='copied',extent='whole_file')
 if not z.get('sha256'):z['sha256']=z['archive_sha256']
# Captures with alternate destination keep their metadata; map by copy_to when needed.
for pp in (B/'audit').glob('*plan.json'):
 plan=readj(pp,[])
 if not isinstance(plan,list):continue
 for x in plan:
  if not isinstance(x,dict) or 'path' not in x:continue
  orig=x['path'];dest=x.get('copy_to')
  if dest and (B/dest).is_file():
   p=B/dest;z=sources.setdefault(orig,{'source_id':uid('src',orig),'path':orig})
   z.update(archive_path=dest,archive_size=p.stat().st_size,archive_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),availability='copied',extent='whole_file')
for p in sorted((B/'audit').glob('external_sources*.json')):
 for x in readj(p,[]):
  z=dict(x);z.update(source_id=uid('src',x['url']),path=x['url'],extent='whole_file');sources[x['url']]=z

aliases={
'S_CAL':A+'/src/Calibration.cxx','S_MAP':A+'/src/EBUdecode.cxx','S_DECODE':A+'/src/Decode.cxx','S_PED':A+'/src/Pedestal.cxx','S_CONFIG':A+'/run/global_config','S_INTER':A+'/src/InterCalib.cxx',
'S_COMB':R+'/analysis/sipm_gain/common/CombFit.cc','S_GAIN_README':R+'/analysis/sipm_gain/README.md','S_GAIN_SELECT':R+'/analysis/sipm_gain_selection/README.md','S_MIP_NEW':R+'/analysis/mip_calibration/README.md',
'S_SINGLE':R+'/analysis/resolution/SingleEnergyAnalysis.cc','S_MULTI':R+'/analysis/resolution/MultiEnergyAnalysis.cc','S_RESLOG':R+'/analysis/log/resolution/compare/test-1326323.out',
'S_G4MAIN':S+'/CEPCScECAL/CEPCScECAL.cc','S_GEOM':S+'/CEPCScECAL/src/DetectorConstruction.cc','S_STEP':S+'/CEPCScECAL/src/SteppingAction.cc','S_EVENT':S+'/CEPCScECAL/src/EventAction.cc','S_PRIMARY':S+'/CEPCScECAL/src/PrimaryGeneratorAction.cc',
'S_DIGI':S+'/Analysis_edit/src/MCDigi.cxx','S_EXTRACT':S+'/Analysis_edit/src/Extract.cxx','S_GEN_SPS':S+'/CEPCScECAL/generate/generator_sps_e.sh','S_GEN_POINT':S+'/CEPCScECAL/generate/generator_e.sh',
'S_OLD_TABLE':'/home/murata_t/scecal_reassessment_20260929/saved_fit_values.csv'}
for alias,path in aliases.items():
 z=sources.setdefault(path,{'source_id':uid('src',path),'path':path,'availability':'unknown','hash_reason':'not captured or capture incomplete'})
 z.setdefault('aliases',[]).append(alias)

campaigns=[{'campaign_id':'ps2023','facility':'CERN PS','beamline':'T9','site_label':'157','period_presentation':'2023-05-17/2023-05-31','period_schedule':'2023-05-16/2023-05-31','period_difference':'setup versus acquisition definition unconfirmed','status':'暫定','evidence_status':'確認済み','source_locator':'L6 slide 6; L8 slide 7'}, {'campaign_id':'sps2023','facility':'CERN SPS','beamline':'H2','site_label':'887','period_presentation':'2023-04-26/2023-05-10','period_schedule':'2023-04-24/2023-05-10','period_difference':'setup versus acquisition definition unconfirmed','status':'暫定','evidence_status':'確認済み','source_locator':'L6 slide 6; L8 slide 7'}, {'campaign_id':'ps_sps2023','facility':'CERN PS/SPS','beamline':'mixed or not individually established','status':'暫定','evidence_status':'不明・未検証'}]
samples={};conditions={};results=[];runs=[];figures=[]
def sample(path,kind='unknown',role='保存結果参照'):
 sid=uid('sample',path)
 if sid in samples:return sid
 cp='sps2023' if '/sps/' in path.lower() else 'ps2023' if '/ps/' in path.lower() else 'ps_sps2023'
 particle=next((p for p in ['e-','mu-','pi-','p+'] if '/'+p+'/' in path or ('_'+p+'_') in path),'unknown')
 samples[sid]={'sample_id':sid,'campaign_id':cp,'data_or_mc':kind,'particle':particle,'input_path':path,'role':role,'beam_quantity_definition':'unknown unless explicitly given by production record','status':'暫定'}
 return sid

def condition(label,source,kind='unknown',selection='unknown'):
 cid=uid('condition',label)
 sources.setdefault(source,{'source_id':uid('src',source),'path':source,'availability':'metadata_only','sha256':None,'hash_reason':'Directory listing or read-only export; original whole-file hash not obtained'})
 if cid not in conditions:
  sid=sample(source,kind)
  conditions[cid]={'condition_id':cid,'label':label,'sample_ids':[sid],'source_ids':[uid('src',source)],'calibration_set_id':'unknown','generation_source_version':'unknown','generation_config':'unknown','reconstruction_id':'unknown','digitization_id':'unknown','selection_id':'unknown','selection_description':selection,'evaluation_id':'unknown','evidence_status':'確認済み','status':'暫定','limitation':'Saved result readout; current implementation is not substituted for its generation version.'}
 return cid

resultcols='result_id,campaign_id,sample_id,condition_id,comparison_group,data_or_mc,metric,value,uncertainty,uncertainty_kind,unit,source,source_locator,value_origin,evidence_status,status,valid,invalid_reason,events,selection,notes'.split(',')
def result(source,locator,metric,value,cid,kind='unknown',unit='',error='',origin='saved_table',valid='unknown',reason='',events='',notes=''):
 if value is None:return
 try:
  if isinstance(value,str) and value.lower() in ['nan','inf','-inf','none','null','']:return
  if isinstance(value,float) and not math.isfinite(value):return
 except Exception:return
 sid=conditions[cid]['sample_ids'][0]
 results.append(dict(result_id=uid('result',source+'|'+locator+'|'+metric+'|'+cid),campaign_id=samples[sid]['campaign_id'],sample_id=sid,condition_id=cid,comparison_group=uid('group',conditions[cid]['label']),data_or_mc=kind,metric=metric,value=value,uncertainty=error if error is not None else '',uncertainty_kind='fit_parameter_error; covariance not exported' if error not in ['',None] else 'unknown',unit=unit,source=source,source_locator=locator,value_origin=origin,evidence_status='確認済み',status='暫定',valid=valid,invalid_reason=reason,events=events,selection=conditions[cid]['selection_description'],notes=notes))

# Existing summary tables. Read a bounded number of rows, never the millions of channel/method rows.
summary_names={'response_comparison.tsv','common_group_summary.tsv','reader_status.tsv','summary.tsv','impact.tsv','model_summary.tsv','truth_reference.tsv','saved_fit_values.csv','fit_summary.tsv','metrics.tsv','closure_summary.tsv'}
for orig,z in sorted(sources.items()):
 rel=z.get('archive_path');name=pathlib.Path(orig).name
 if not rel or name not in summary_names:continue
 if name=='saved_fit_values.csv':continue
 p=B/rel
 if p.stat().st_size>3000000:continue
 try:
  with p.open() as f:
   reader=csv.DictReader(f,delimiter='\t' if name.endswith('.tsv') else ','); rows=[]
   for i,row in enumerate(reader):
    if i>=1500:break
    rows.append(row)
  if len(rows)>=1500:continue
 except Exception:continue
 # Do not turn channel constants into results.
 if rows and any(c in rows[0] for c in ['CellID','cell_id','cellID','cellid']) and len(rows)>100:continue
 for i,row in enumerate(rows,2):
  labels={k:v for k,v in row.items() if k and v is not None and (k.lower() in ['energy','energy_gev','sipm','sipm_type','category','group','model','sample','cog','channel','layer','mode','population','match','matching','kind','sipm_um','chip','old_fallback_class','variant','threshold_shift'] or not re.fullmatch(r'[-+]?\d+(\.\d*)?([eE][-+]?\d+)?',str(v)))}
  group=orig+':'+json.dumps(labels,ensure_ascii=False,sort_keys=True)
  kind=row.get('sample',row.get('data_or_mc','unknown'))
  kind='mc' if kind in ['sim','simulation','mc'] else 'data' if kind=='data' else 'unknown'
  cid=condition(group,orig,kind,selection='table-specific selection; see copied inputs/config/source; no cross-condition merging')
  conditions[cid]['table_row_context']=labels
  for key,val in row.items():
   if key is None or val is None:continue
   if key.lower() in ['energy','energy_gev','cog','channel','layer','cellid','cell_id','run','run_id','index','sipm_um','chip','old_fallback_class','variant']:continue
   try:num=float(val)
   except (ValueError,TypeError):continue
   if not math.isfinite(num):continue
   unit=''
   if name=='response_comparison.tsv':unit='GeV' if 'mean' in key and 'ratio' not in key else 'dimensionless'
   elif name=='reader_status.tsv':unit='ADC' if 'mean' in key or 'mpv' in key else 'count' if any(k in key for k in ['accepted','fallback','channels']) else ''
   elif name=='saved_fit_values.csv':unit='MeV' if key in ['e_fit','e_sigma','e_fit_error','e_mean','e_fit_min','e_fit_max'] else 'hits' if key in ['n_fit','n_mean','n_fit_error'] else ''
   error=row.get(key+'_error','')
   result(orig,f'line {i}; column {key}',key,val,cid,kind,unit,error,events=row.get('entries',''),notes='Saved table verbatim. Units only assigned where definition is confirmed; blank units require source interpretation.')

# Selected numeric summaries retain counts in their declared scope (not unique channels/events).
for orig,z in sorted(sources.items()):
 if pathlib.Path(orig).name not in ['summary.json','reader_status.json']:continue
 rel=z.get('archive_path')
 if not rel:continue
 obj=readj(B/rel,{})
 if not isinstance(obj,dict):continue
 if not any(x in orig for x in ['sipm_gain','mip_calibration','muon_mip']):continue
 cid=condition('saved_summary:'+orig,orig)
 def scalars(o,p=''):
  if not isinstance(o,dict):return
  for k,v in o.items():
   key=p+'.'+k if p else k
   if isinstance(v,(int,float)) and not isinstance(v,bool):result(orig,'JSON '+key,key,v,cid,origin='saved_summary',notes='Scope follows summary key; no conversion from channel-conditions to unique channels.')
   elif isinstance(v,dict) and key.count('.')<2:scalars(v,key)
 scalars(obj)

# Raw directory inventory: names only; no event decoding or wholesale raw hashing.
inv=readj(B/'audit/raw_run_inventory.json',[])
for item in inv:
 d=item.get('directory','')
 if '/2023/' not in d or not any('/'+c+'/' in d for c in ['ps','sps']):continue
 for n in item.get('files',[]):
  path=d+'/'+n
  if not re.search(r'\.(dat|raw|bin|root)$',n,re.I):continue
  sid=sample(path,'data','取得済みraw候補;用途未確認')
  m=re.search(r'(?:Run|run)[_-]?(\d+)',n);energy=re.search(r'/([0-9.]+)GeV(?:/|$)',path)
  part=samples[sid]['particle'];charge='-1' if part in ['e-','mu-','pi-'] else '+1' if part=='p+' else ''
  cid=condition('raw_inventory:'+samples[sid]['campaign_id'],d,'data','none; directory inventory only')
  row=dict(campaign_id=samples[sid]['campaign_id'],sample_id=sid,run_id=m.group(1) if m else '',data_or_mc='data',particle=part,charge=charge,beam_momentum='',beam_momentum_unit='',beam_energy='',beam_energy_unit='',beam_energy_definition='unknown; filename/directory label is not a verified physical quantity',beam_label=energy.group(1)+'GeV' if energy else '',use='取得済みraw候補;較正/調整/性能/除外の区別は未確認',input_path=path,recorded_events='',reconstructed_events='',selected_events='',condition_id=cid,selection_id='unknown',source=d,source_locator='directory entry '+n,evidence_status='確認済み',associated_data_run='',notes='Bounded directory listing, not acquisition log. Do not sum repeated paths or aliases.')
  if '/pedestal/' in path:row['use']='取得済み;専用pedestal候補'
  runs.append(row)

# All small copied input manifests are kept intact; add each input path as a usage index.
seen=set()
for orig,z in sorted(sources.items()):
 rel=z.get('archive_path');name=pathlib.Path(orig).name.lower()
 if not rel or not any(x in name for x in ['manifest','input','run_list','runlist']):continue
 p=B/rel
 if p.suffix.lower() not in ['.json','.tsv','.csv','.txt'] or p.stat().st_size>3000000:continue
 try:lines=p.read_text().splitlines()
 except Exception:continue
 for i,line in enumerate(lines,1):
  for path in re.findall(r'/(?:home|megraid01|net)/[^\s\"\t,]+?\.(?:root|dat)',line):
   if (orig,path) in seen:continue
   seen.add((orig,path))
   if '/2022/' in path or '/2025/' in path:continue
   kind='mc' if any(t in path.lower() for t in ['/simulation/','/generate/','/digi/','/threshold_mpvwoct/','/threshold/']) else 'data'
   role='較正入力' if any(t in orig for t in ['sipm_gain','mip_calibration']) else '診断/比較入力;parameter調整・検証の独立性未確定'
   sid=sample(path,kind,role);cid=condition('input_manifest:'+orig,orig,kind,'manifest-specific; see source row')
   part=samples[sid]['particle'];m=re.search(r'(?:Run|run)[_-]?(\d+)',path);energy=re.search(r'/([0-9.]+)GeV(?:/|$)',path)
   runs.append(dict(campaign_id=samples[sid]['campaign_id'],sample_id=sid,run_id=m.group(1) if m and kind=='data' else '',data_or_mc=kind,particle=part,charge='-1' if part in ['e-','mu-','pi-'] else '+1' if part=='p+' else '',beam_momentum='',beam_momentum_unit='',beam_energy='',beam_energy_unit='',beam_energy_definition='unknown; nominal label preserved separately',beam_label=energy.group(1)+'GeV' if energy else '',use=role,input_path=path,recorded_events='',reconstructed_events='',selected_events='',condition_id=cid,selection_id='unknown',source=orig,source_locator=f'line {i}',evidence_status='確認済み',associated_data_run=m.group(1) if m and kind=='mc' else '',notes='Manifest row is a usage, not a new acquisition. Raw event count not inferred. MC Run label denotes association only.'))

# Root export records remain fully available. Fit-specific normalization is handled by a separate adapter.
rootexports=[]
for p in sorted((B/'evidence/root_exports').glob('*.json')) if (B/'evidence/root_exports').exists() else []:
 rootexports.append({'export_id':uid('export',str(p.relative_to(B))),'archive_path':str(p.relative_to(B)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'size':p.stat().st_size,'kind':'read_only_saved_ROOT_object_export','not_new_fit':True})

# Actual copied images/PDFs and explicit requests.
for orig,z in sorted(sources.items()):
 rel=z.get('archive_path')
 if not rel or pathlib.Path(rel).suffix.lower() not in ['.png','.pdf','.jpg','.jpeg','.svg']:continue
 cid=condition('figure:'+orig,orig,selection='unknown; filename correspondence is not generation proof')
 fid=uid('fig',orig)
 figures.append(dict(figure_id=fid,row_type='file',original_path=orig,archive_path=rel,sha256=z.get('archive_sha256',z.get('sha256','')),description='既存図・配布資料: '+pathlib.Path(orig).name,manuscript_use='装置/較正/評価の候補。本文G節と原図captionを照合',campaign_id='ps_sps2023',sample_id=conditions[cid]['sample_ids'][0],condition_id=cid,result_ids='',generation_source='unknown',generation_command='unknown',input_paths='unknown',source_locator='whole file; PDF page identification in HANDOFF/manifest',evidence_status='推定',availability='copied',limitations='Unmodified existing asset. Thesis adoption and exact input/result correspondence remain provisional.'))
for name,desc in [('setup2023','2023 PS/SPSそれぞれのas-built setup、trigger/DAQ同期、位置・角度・外部PIDを示す図と設定log'),('pedestal2023','実際に使用したpedestal2023_SPS.rootに対応するHG/LG fit例、run・channel・memory・fit範囲付き'),('gain2023','production gain tableの生成run/DAC/温度と対応するp.e. peak fit。新候補図との区別'),('hlratio_v3','all_hl_electron2023_hlratio_v3.rootの生成producerと横/縦軸・pedestal・fit範囲が分かるHG/LG図'),('mip_v4','production MIP v4 trackfitに対応するmuon/path length/threshold-aware fit図'),('linearity','同一productionとselectionによるresponse/linearity/residual、scale合わせと誤差の定義付き'),('mapping_asbuilt','32層の実機構成、追加2層の位置と30層解析への対応、SiPM型番とlayer map')]:
 figures.append(dict(figure_id='request_'+name,row_type='request',description=desc,manuscript_use='第3〜5章',campaign_id='ps_sps2023',sample_id='unknown',condition_id='unknown',result_ids='',generation_source='unknown',generation_command='unknown',input_paths='候補: 本文記載のsource/table/発表資料',source_locator='探索済み範囲はHANDOFF A/G節',evidence_status='不明・未検証',availability='not_identified',limitations='未発見を不存在と断定しない。新規図は今回作成していない。'))

repositories=[]
for name,path in [('research',R),('analysis','/home/murata_t/data_beamtest/analysis'),('simulation','/home/murata_t/data_beamtest/simulation')]:
 o=readj(B/'audit'/('git_'+name+'.json'),{})
 repositories.append({'repository_id':name,'path':path,'resolved_path':o.get('resolved_path','unknown'),'branch':o.get('branch','unknown').strip(),'commit':o.get('commit','unknown').strip(),'dirty_state':'unknown' if o.get('status_error') or 'status' not in o else 'dirty' if o.get('status','').strip() else 'clean_tracked_only','untracked_scope':'excluded from status query','status_error':o.get('status_error'),'source_hashes':'sources entries; actual copies hashed independently'})
repositories.append({'repository_id':'data_root','path':'/home/murata_t/data_beamtest','resolved_path':'/net/megsrv01/export/megraid01/users/data_beamtest','git_managed':False,'branch':None,'commit':None,'dirty_state':'not_applicable'})

calibrations=[]
for key,name,tree,role in [('ped_sps','pedestal2023_SPS.root','ChnLevel','pedestal'),('hl_v3','all_hl_electron2023_hlratio_v3.root','intercalibration tree; inspect exported schema','HG/LG'),('mip_v4','all_auto_muon_v4_trackfit.root','MIP_Fit','MIP'),('gain_old','20230504_1425_Calib_LED_all.root','T_Event','SiPM gain')]:
 path=S+'/Analysis_edit/share/'+name
 calibrations.append({'calibration_set_id':key,'role':role,'file':path,'tree':tree,'current_reader_aliases':['S_CAL','S_EXTRACT'],'generation_source':'unknown','generation_run':'unknown','reference_temperature':'unknown','actual_application_by_result':'unknown except explicit input manifests','status':'暫定','evidence_status':'確認済み'})
calibrations.append({'calibration_set_id':'threshold','role':'threshold turn-on','file':R+'/analysis/result/threshold/threshold.root','tree':'expErfThre','generation_source':'unknown','limitation':'electron-spectrum model, not independent hardware turn-on','status':'暫定'})
calibrations.append({'calibration_set_id':'mip_refit_candidate','role':'candidate MIP','file':R+'/analysis/result/mip_calibration/all_muons_refit/mip.root','generation_record':'copied inputs/config/source/COMPLETE','adopted_for_production':'unknown','status':'暫定'})

references=[
 {'reference_id':'L1','title':'Design of Sc-ECAL Prototype for CEPC and Performance of First Two Layers','author':'Yazhou Niu','journal':'JINST 15 (2020) C05036','doi':'10.1088/1748-0221/15/05/C05036','url':'https://arxiv.org/pdf/2002.01809v2','locator':'PDF pp.3–4, fig.1; calibration figs.3–5 are early two-layer prototype','scope':'2020 design; not 2023 measurements'},
 {'reference_id':'L2','title':'CALICE scintillator-SiPM calorimeter prototypes: R&D highlights and beamtests','authors':['Jianbei Liu','Yong Liu','Wataru Ootani','Tohru Takeshita','Haijun Yang','Yunlong Zhang'],'journal':'NIM A 1072 (2025) 170191','doi':'10.1016/j.nima.2024.170191','url':'https://www.sciencedirect.com/science/article/abs/pii/S0168900224011173','locator':'abstract/prototype overview','scope':'overview including 2022 and 2023; not imported as 2023 performance'},
 {'reference_id':'L3','author':'Yong Liu','title':'CERN Beamtests of CALICE Scintillator Calorimeter Prototypes','date':'2023-09-06','venue':'TIPP2023','url':'https://indico.tlabs.ac.za/event/112/contributions/2961/attachments/1161/1577/2023_0906_CERN_Beamtests_CALICE_ScCalo_Prototypes.pdf','locator':'slides 12,24','scope':'2023 campaign and LED figure credit'},
 {'reference_id':'L4','author':'Wataru Ootani','title':'Status of ScW-ECAL','date':'2025-04-01','venue':'DRD6','url':'https://indico.cern.ch/event/1531732/contributions/6448349/attachments/3042846/5375711/DRD_ScWECAL20250401WO.pdf','locator':'slide 4','scope':'later apparatus overview; exclude 2025 measurement results'},
 {'reference_id':'L6','author':'Tatsuki Murata','title':'Status of high granular Sc-ECAL for future electron positron colliders','venue':'CALOR2024','url':'https://indico.cern.ch/event/1339557/contributions/5898550/attachments/2859747/5003008/calor2024TatsukiMurata.pdf','locator':'slide 6','scope':'2023 campaign overview; exact settings need run logs'},
 {'reference_id':'L7','author':'Tatsuki Murata','title':'ILC-Japan presentation; verify exact title on cover','url':'https://agenda.linearcollider.org/event/10229/contributions/53761/attachments/39257/61905/ILC-Japan_murata.pdf','locator':'cover and campaign slides; page verification pending','scope':'photo/campaign reference candidate'},
 {'reference_id':'L8','title':'BeamTest-20230928.pdf; verify exact cover title and authors','url':'https://agenda.linearcollider.org/event/10098/contributions/53294/attachments/39113/61631/BeamTest-20230928.pdf','locator':'slide 7','scope':'campaign schedule; setup/acquisition distinction unconfirmed'}]

manifest={'schema_version':'scecal-thesis-handoff-1','snapshot_id':B.name,'created_at':'2026-10-24T09:16:34+09:00','assembled_at':now(),'timezone':'Asia/Tokyo','scope':{'included':'ScECAL 2023 PS/SPS apparatus, calibration, reconstruction, simulation/digitization, saved evaluation','excluded':['2022 SPS measurements','2025 KEK measurements','GNN/PFA/SiECAL other study'],'cross_year_dependencies':'older constants/comments and later apparatus references are separately identified'},'methods':{'read_only_originals':True,'new_event_processing':False,'new_fits':False,'root_readout':'saved objects/metadata/parameters only','atomic_snapshot':False,'coverage':'targeted source/result/2023 raw directories; bounded depth; no claim of exhaustive server search','network_filesystem':'NFS; read delays/timeouts recorded in audit','csv_missing':'empty is unknown, never zero'},'status_terms':{'evidence_status':['確認済み','推定','不明・未検証'],'status_default':'暫定','valid':['true','false','unknown']},'repositories':repositories,'sources':list(sources.values()),'source_aliases':{k:uid('src',v) for k,v in aliases.items()},'campaigns':campaigns,'samples':list(samples.values()),'calibration_sets':calibrations,'reconstructions':[{'reconstruction_id':'reco_current','source_alias':'S_CAL','status':'暫定','generation_assignment':'Not assigned to past results without execution evidence','schema':'Raw_Hit to Calib_Hit','formula':'qH=rawHG-pedHG; qL=rawLG-pedLG; C=1-(T-20)*alpha; reject C*qH<10; E=0.305*(rawHG<2600?C*qH:C*qL/s-b/s)/M [MeV]','source_locator':'Calibration.cxx L560–582,659–709'}],'simulations':[{'simulation_id':'g4_current','source_aliases':['S_G4MAIN','S_GEOM','S_STEP','S_EVENT','S_PRIMARY'],'physics_list':'QGSP_BERT','layers':30,'production_version':'unknown for historical samples','status':'暫定'}],'digitizations':[{'digitization_id':'digi_current','source_aliases':['S_DIGI','S_EXTRACT'],'source_locator':'MCDigi.cxx L141–302,304 onward','production_assignment':'unknown except frozen diagnostic input bundles','light_yield':'(1-pCT)*M/G in captured current source','status':'暫定'}],'selections':[{'selection_id':'electron_current','source_alias':'S_SINGLE','source_locator':'SingleEnergyAnalysis.cc; CoG and eventCutConditions','generation_assignment':'unknown without matching log/config','status':'暫定'}],'evaluations':[{'evaluation_id':'resolution_current','source_aliases':['S_SINGLE','S_MULTI'],'estimators':['Gaussian sigma/mu','RMS/mean; separately labeled'],'beam_spread_correction':'not confirmed','full_systematics':'not found','status':'暫定'}],'conditions':list(conditions.values()),'figures':figures,'results_file':'results.csv','runs_file':'runs.csv','figures_file':'figures.csv','root_exports':rootexports,'references':references,'unresolved':[{'issue_id':'u_layers','question':'32-layer as-built versus 30-layer analysis/MC; locate inventory and exclusion rationale','impact':'chapters 3,4,5,9,10'}, {'issue_id':'u_constants','question':'generation source/run/temperature and PS applicability of production tables','impact':'chapters 4,5,9,10'}, {'issue_id':'u_temperature','question':'coefficient sign/reference and pre/post-fix production mapping','impact':'chapters 4,5'}, {'issue_id':'u_ct','question':'M/G definition, CT normalization, tuning and independent validation samples','impact':'chapters 5,9,10'}, {'issue_id':'u_adoption','question':'author choice of final conditions and contribution boundary','impact':'all thesis claims'}],'handoff_destination':'docs/research/snapshots/'+B.name+'/; preserve originals and create docs/research/scecal_context.md separately'}
# JSON is valid YAML 1.2 and avoids adding a runtime dependency to assemble this bundle.
save_text('manifest.yaml',json.dumps(clean(manifest),ensure_ascii=False,indent=2)+'\n')
dumpcsv('results.csv',results,resultcols)
runcols='campaign_id,sample_id,run_id,data_or_mc,particle,charge,beam_momentum,beam_momentum_unit,beam_energy,beam_energy_unit,beam_energy_definition,beam_label,use,input_path,recorded_events,reconstructed_events,selected_events,condition_id,selection_id,source,source_locator,evidence_status,associated_data_run,notes'.split(',')
dumpcsv('runs.csv',runs,runcols)
figcols='figure_id,row_type,original_path,archive_path,sha256,description,manuscript_use,campaign_id,sample_id,condition_id,result_ids,generation_source,generation_command,input_paths,source_locator,evidence_status,availability,limitations'.split(',')
dumpcsv('figures.csv',figures,figcols)
lines=['# Source aliases and captured evidence','', '本文のsource aliasを絶対path・同梱実体・hashに結ぶ索引。全採取記録はmanifest.yamlを参照。未取得のfileを同梱済みと解釈しない。','']
for alias,path in aliases.items():
 z=sources[path];lines+=['## '+alias,'','Original: `'+path+'`','', 'Source ID: `'+z['source_id']+'`','',('Bundled: ['+z.get('archive_path','')+']('+z.get('archive_path','')+')') if z.get('archive_path') else 'Bundled: unknown / not captured','', 'SHA-256: `'+str(z.get('archive_sha256',z.get('sha256','unknown')))+'`','']
save_text('SOURCES.md','\n'.join(lines))
save_text('audit/assembly_counts.json',json.dumps({'sources':len(sources),'results':len(results),'runs':len(runs),'figures':len(figures),'conditions':len(conditions),'samples':len(samples),'assembled_at':now()},indent=2))
for name,future in OUTPUT_JOBS:
 future.result();print('output complete',name,flush=True)
OUTPUT_POOL.shutdown(wait=True)
print('assembled',len(sources),'sources',len(results),'results',len(runs),'runs',len(figures),'figures',flush=True)
