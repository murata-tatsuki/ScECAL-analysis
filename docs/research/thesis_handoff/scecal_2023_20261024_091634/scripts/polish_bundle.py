"""Final semantic checks and references using copied evidence only; no event analysis."""
import pathlib,json,csv,hashlib,re,datetime
B=pathlib.Path(__file__).absolute().parents[1]
def uid(k,s):return k+'_'+hashlib.sha256(s.encode()).hexdigest()[:12]
def readcsv(n):
 with (B/n).open() as f:
  r=csv.DictReader(f);return list(r),r.fieldnames
def writecsv(n,rows,cols):
 with (B/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)
m=json.loads((B/'manifest.yaml').read_text());rr,rc=readcsv('results.csv');runs,runc=readcsv('runs.csv');figs,fc=readcsv('figures.csv')
src={s['source_id']:s for s in m['sources']};ss={s['sample_id']:s for s in m['samples']};cc={c['condition_id']:c for c in m['conditions']}
legacy=[r for r in rr if r['source'].startswith('/home/murata_t/scecal_reassessment_20260929/')]
rr=[r for r in rr if not r['source'].startswith('/home/murata_t/scecal_reassessment_20260929/')]
m['legacy_result_tables']=[{'source_id':uid('src','/home/murata_t/scecal_reassessment_20260929/saved_fit_values.csv'),'result_index_inclusion':False,'reason':'Original ROOT moved/missing; campaign/generation chain incomplete. Copied transcript retained separately, not merged into the 2023 scan.'}]
for r in rr:
 if r['metric']=='gaussian_sigma_over_mu':r['uncertainty_kind']='saved graph error; Gaussian-parameter propagation, covariance/systematics unverified'
 c=cc[r['condition_id']];c['comparison_group']=r['comparison_group']
 if c['label'].startswith('comparison_oct12_series_'):
  sid=uid('sample',r['source']+'#'+r['condition_id'])
  ss.setdefault(sid,{'sample_id':sid,'campaign_id':r['campaign_id'],'data_or_mc':r['data_or_mc'],'particle':'e-','input_path':r['source'],'role':'saved comparison cohort; series/CoG distinguished by condition; raw-run composition not inferred','status':'暫定'})
  r['sample_id']=sid;c['sample_ids']=[sid]
 metric=r['metric'].lower()
 if r['data_or_mc']=='unknown':
  kind='comparison' if metric.startswith(('data_mc','data_over_mc')) else 'data' if metric.startswith('data_') else 'mc' if metric.startswith(('mc_','sim_')) else None
  if kind:
   sid=uid('sample',r['source']+'#'+r['condition_id']+'#'+kind)
   ss.setdefault(sid,{'sample_id':sid,'campaign_id':r['campaign_id'],'data_or_mc':kind,'particle':'mu-' if 'muon' in r['source'] else 'e-' if '/integrated/' in r['source'] else 'unknown','input_path':r['source'],'role':'summary cohort; raw-run mapping in copied input manifests','status':'暫定'})
   r['sample_id']=sid;r['data_or_mc']=kind
   if sid not in cc[r['condition_id']]['sample_ids']:cc[r['condition_id']]['sample_ids'].append(sid)
for s in src.values():
 if not s.get('path','').endswith('/event_audit.tsv') or not s.get('archive_path'):continue
 with (B/s['archive_path']).open() as f:
  for line,row in enumerate(csv.DictReader(f,delimiter='\t'),2):
   for r in runs:
    if pathlib.Path(r['input_path']).name==pathlib.Path(row['file']).name and '/decode/' in r['input_path'] and r['data_or_mc']=='data':
     r['reconstructed_events']=row['entries'];r['notes']+='; saved decode entries from '+s['path']+':line '+str(line)+'; not acquisition count or final selection'
for r in runs:
 c=cc[r['condition_id']]
 if c['label'].startswith('raw_inventory:'):
  sid=uid('src',r['source'])
  if sid in src and sid not in c['source_ids']:c['source_ids'].append(sid)
  c['source_index_archive_path']='audit/raw_run_inventory.json'
roots=json.loads((B/'audit/root_readout.json').read_text())
for c in m['calibration_sets']:
 c['source_id']=uid('src',c['file']);o=next((o for o in roots if o['path']==c['file']),None)
 if o:
  c['saved_tree_exports']=o.get('tables',[]);c['source_metadata']={k:o.get(k) for k in ['size','mtime','uuid','read_started_at','read_finished_at']}
  if c['calibration_set_id']=='hl_v3':c['tree']='InterCalib'
  if c['calibration_set_id']=='ped_sps':
   c['stored_value_branch_type']='std::vector<double>';c['captured_reader_value_type']='std::vector<float>';c['reader_type_note']='Extract.cxx L57–58 documents matching processCali rounding; stored branch type differs from reader type.'
comparison=[r for r in rr if '/resolution/comparsion/comparison_5mm.root' in r['source']]
for f in figs:
 p=f.get('original_path','')
 if p.endswith('/resolution.png') or p.endswith('/resolution_with_noise.png'):
  noise=p.endswith('/resolution_with_noise.png');prefix='with_noise_' if noise else 'without_noise_'
  f.update(result_ids=';'.join(r['result_id'] for r in comparison if r['metric']=='gaussian_sigma_over_mu' or r['metric'].startswith(prefix)),generation_source='/home/murata_t/ScECAL_BeamTest/analysis/resolution/MultiEnergyAnalysis.cc (compatible captured producer; generation version unknown)',input_paths='/home/murata_t/ScECAL_BeamTest/analysis/result/resolution/comparsion/comparison_5mm.root',source_locator='Saved ROOT resolution canvas; S_RESLOG lines3–9 and timestamps',evidence_status='推定',limitations='Saved canvas/name/timestamp/log correspondence; frozen generation source not established. Physics adoption remains provisional.')
for c in cc.values():
 for sid in c.get('source_ids',[]):
  p=src.get(sid,{}).get('path','')
  if '/20261013/' in p:
   family=p.split('/20261013/')[0]+'/20261013/'
   c['diagnostic_generation_evidence_ids']=[x['source_id'] for x in src.values() if x.get('path','').startswith(family) and any(k in pathlib.Path(x['path']).name for k in ['manifest','generation_conditions','configuration','source_sha'])]
   c['diagnostic_evidence_note']='Saved adjacent records; no assignment to older production or independent-validation claim.'
for ref in m['references']:
 if ref['reference_id']=='L8':ref.update(title='2023 Beam Test of Sci-W ECAL and AHCAL Prototypes',author='Yunlong Zhang',locator='slides7–11,18–19,27–29; footer date ambiguity retained')
 if ref['reference_id']=='L3':ref['offline_availability']='PDF retrieval failed; URL and inspected presentation locators retained'
m['samples']=list(ss.values());m['conditions']=list(cc.values());m['figures']=figs
m['capture_resume']['resumed_on']='2026-10-26 Asia/Tokyo';m['completed_at']=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
m['methods']['source_epoch']='Evidence captured 2026-10-24; assembly/verification resumed 2026-10-26. No atomic or latest-source claim.'
m['csv_columns']={'x_value':'saved graph horizontal coordinate','x_unit':'GeV nominal energy-axis label','x_definition':'label definition, no momentum/energy conversion','reconstructed_events':'saved decode entries from event_audit.tsv where filled; not acquisition or final-selection count'}
(B/'manifest.yaml').write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n');writecsv('results.csv',rr,rc);writecsv('runs.csv',runs,runc);writecsv('figures.csv',figs,fc)
p=B/'HANDOFF.md';s=p.read_text()
s=s.replace('afterEventCutはさらにenergy–Nhitの事象選択後。','afterEventCutはさらにenergy–Nhitの事象選択後。current eventCutConditions L115–118は ((E−μE)/(3σE))²+((N−μN)/(3σN))²>1 を除外し、境界=1を残す。')
s=s.replace('table内pathと現在のscript参照pathが異なり、元ROOTの移動/欠落と生成版を要確認','table内pathと現在のscript参照pathが異なり、元ROOTの移動/欠落・campaign・生成版を要確認。この旧tableは同梱するがresults.csvの2023性能点には統合していない')
addition='\n最終組立の補足：全channel定数はexport JSONとschemaへ、保存比較graph等はresults.csvへ収録する。linear/logおよび同じ点の複数panelは独立測定として重複させない。L3のPDFは取得失敗のためURLと参照箇所を保持する。\n'
if addition not in s:s+=addition
p.write_text(s)
counts={'results':len(rr),'runs':len(runs),'figures':len(figs),'sources':len(m['sources']),'conditions':len(cc),'excluded_legacy_metric_rows':len(legacy),'completed_at':m['completed_at']}
(B/'audit/final_index_counts.json').write_text(json.dumps(counts,indent=2));print(json.dumps(counts),flush=True)
