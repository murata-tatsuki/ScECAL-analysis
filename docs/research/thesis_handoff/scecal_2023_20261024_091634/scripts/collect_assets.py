import pathlib,json,os,requests,hashlib,datetime
B=pathlib.Path(__file__).resolve().parents[1];R='/home/murata_t/ScECAL_BeamTest/'
plans=[]
for p,n in [('analysis/result/resolution/figures/comparsion/dataThre/resolution.png','resolution.png'),('analysis/result/resolution/figures/comparsion/dataThre/resolution_with_noise.png','resolution_with_noise.png'),('analysis/result/resolution/figures/comparsion/dataThre/summary_100GeV.png','energy100_summary.png'),('analysis/result/resolution/figures/comparsion/dataThre/summary_5GeV.png','energy5_summary.png')]:plans.append({'path':R+p,'copy_to':'assets/results/'+n,'max_size':15000000})
roots=['analysis/result/mip_calibration/all_muons_refit/figures','analysis/result/sipm_gain/full_grouped/comparison/plots/layer_00','investigations/result/e_ssa_threshold_review_20260930/20261013/muon_ct','investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/model_scan']
for root in roots:
 try:
  imgs=[n for n in os.listdir(R+root) if n.endswith(('.png','.pdf'))]
  for n in imgs[:4]:plans.append({'path':R+root+'/'+n,'copy_to':'assets/results/'+root.split('/')[-2]+'_'+n,'max_size':15000000})
 except Exception:pass
for root in ['investigations/result/e_ssa_threshold_review_20260930/20261013/digitization','investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/inputs','investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/model_scan','investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/model_scan/inputs','investigations/result/e_ssa_threshold_review_20260930/20261013/integrated/inputs','investigations/result/e_ssa_threshold_review_20260930/20261013/muon_ct/inputs']:
 try:
  for n in os.listdir(R+root):
   if n.endswith(('.tsv','.json','.txt','.sh','.log')):plans.append({'path':R+root+'/'+n})
 except Exception:pass
(B/'audit/assets_plan.json').write_text(json.dumps(plans,ensure_ascii=False,indent=2));print('assets planned',len(plans),flush=True)
items=[('design_2020.pdf','https://arxiv.org/pdf/2002.01809v2'),('tipp_2023.pdf','https://indico.tlabs.ac.za/event/112/contributions/2961/attachments/1161/1577/2023_0906_CERN_Beamtests_CALICE_ScCalo_Prototypes.pdf'),('ootani_2025.pdf','https://indico.cern.ch/event/1531732/contributions/6448349/attachments/3042846/5375711/DRD_ScWECAL20250401WO.pdf'),('murata_calor2024.pdf','https://indico.cern.ch/event/1339557/contributions/5898550/attachments/2859747/5003008/calor2024TatsukiMurata.pdf'),('murata_ilcjapan.pdf','https://agenda.linearcollider.org/event/10229/contributions/53761/attachments/39257/61905/ILC-Japan_murata.pdf'),('calice_beamtest2023.pdf','https://agenda.linearcollider.org/event/10098/contributions/53294/attachments/39113/61631/BeamTest-20230928.pdf'),('liu_2024_overview.pdf','https://indico.ihep.ac.cn/event/22295/contributions/155135/attachments/77354/95987/2024_0416_CEPC_Calorimeters_Overview.pdf')]
out=[]
for n,u in items:
 z={'url':u,'retrieved_at':datetime.datetime.now(datetime.timezone.utc).isoformat()};p=B/'assets/literature'/n
 try:
  if p.is_file() and p.read_bytes().startswith(b'%PDF'):data=p.read_bytes();z['retrieval_note']='existing prior copy retained'
  else:
   r=requests.get(u,timeout=30);z.update(http_status=r.status_code,content_type=r.headers.get('content-type'));r.raise_for_status();data=r.content
   if not data.startswith(b'%PDF'):raise ValueError('not a PDF')
   p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
  z.update(archive_path=str(p.relative_to(B)),availability='copied',size=len(data),sha256=hashlib.sha256(data).hexdigest())
 except Exception as e:z.update(availability='unavailable',error=str(e))
 out.append(z);print(n,z['availability'],flush=True)
 (B/'audit/external_sources_resumed.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
