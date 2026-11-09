#!/usr/bin/env python3
"""Read saved objects only, via uproot. No macro, fit or event processing."""
import concurrent.futures, csv, hashlib, json, math, os, pathlib
from datetime import datetime
from zoneinfo import ZoneInfo
import awkward as ak
import uproot
BASE=pathlib.Path(__file__).resolve().parents[1]
OUT=BASE/'evidence/root_exports'; OUT.mkdir(parents=True,exist_ok=True)
def now():return datetime.now(ZoneInfo('Asia/Tokyo')).isoformat()
def clean(v):
    if isinstance(v,float) and not math.isfinite(v):return None
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if hasattr(v,'tolist'):return clean(v.tolist())
    if isinstance(v,(str,int,float,bool)) or v is None:return v
    return str(v)
def fit(obj):
    m=obj.all_members
    formula=m.get('fFormula')
    params=list(formula.member('fClingParameters')) if formula is not None else []
    return clean({'class':obj.classname,'name':m.get('fName'),'formula':str(formula.member('fFormula')) if formula is not None else None,'parameters':params,'errors':list(m.get('fParErrors',[])),'xmin':m.get('fXmin'),'xmax':m.get('fXmax'),'chi2':m.get('fChisquare'),'ndf':m.get('fNDF'),'fit_status':None,'fit_status_note':'TF1 parameters are saved; minimizer status not independently saved here'})
def read(item):
    p,kind=item; sid=hashlib.sha256(p.encode()).hexdigest()[:12]; rec={'path':p,'kind':kind,'source_id':'root_'+sid,'read_started_at':now(),'sha256':None,'hash_omission_reason':'Large ROOT not streamed for full hash; selected saved objects exported; original metadata retained'}
    try:
        st=os.stat(p);rec.update(size=st.st_size,mtime=datetime.fromtimestamp(st.st_mtime,ZoneInfo('Asia/Tokyo')).isoformat())
        with uproot.open(p) as f:
            rec['objects']=f.classnames(recursive=False);rec['root_file_version']=f.file.fVersion;rec['uuid']=str(f.file.uuid)
            if kind=='constants':
                rec['tables']=[]
                for key,cls in rec['objects'].items():
                    if cls!='TTree':continue
                    t=f[key];r={'object':key,'entries':t.num_entries,'branches':t.typenames()}
                    if t.num_entries>50000:r['export_error']='entry limit';rec['tables'].append(r);continue
                    rows=ak.to_list(t.arrays(library='ak'))
                    for i,row in enumerate(rows):row['_entry']=i
                    target=OUT/(sid+'_'+key.split(';')[0]+'.json')
                    target.write_text(json.dumps(clean(rows),ensure_ascii=False,separators=(',',':')))
                    r['archive_path']=str(target.relative_to(BASE));r['sha256']=hashlib.sha256(target.read_bytes()).hexdigest();r['extent']='all stored tree rows, no aggregation; vector element index retained by list order'
                    rec['tables'].append(r)
            elif kind=='resolution':
                rec['metadata']={};rec['saved_fits']={};rec['histograms']={}
                for key,cls in rec['objects'].items():
                    name=key.split(';')[0]
                    if cls=='TNamed':rec['metadata'][name]=f[name].member('fTitle')
                    elif cls.startswith('TParameter'):rec['metadata'][name]=f[name].member('fVal')
                for stage in ['beforeEventCut','afterEventCut']:
                    if stage not in f:continue
                    d=f[stage]
                    for name in ['fit_gaus','fit_gaus_nhit']:
                        if name in d:rec['saved_fits'][stage+'/'+name]=fit(d[name])
                    for key,cls in d.classnames(recursive=False).items():
                        name=key.split(';')[0]
                        if cls.startswith('TH1') and (name.startswith('edep_') or name.startswith('nhits_')) and 'layer' not in name and '1hit' not in name:
                            h=d[name];m=h.all_members
                            r=clean({'class':cls,'entries':m.get('fEntries'),'nbins':h.axis().member('fNbins'),'xmin':h.axis().member('fXmin'),'xmax':h.axis().member('fXmax'),'sumw':m.get('fTsumw'),'sumwx':m.get('fTsumwx'),'sumwx2':m.get('fTsumwx2')})
                            if pathlib.Path(p).name.startswith(('5GeV_','100GeV_')):
                                r['edges']=h.axis().edges().tolist();r['contents_with_flow']=h.values(flow=True).tolist();r['variances_with_flow']=h.variances(flow=True).tolist()
                            rec['histograms'][stage+'/'+name]=r
            else:
                rec['trees']={k:{'entries':f[k].num_entries,'branches':f[k].typenames()} for k,c in rec['objects'].items() if c=='TTree'}
        en=os.stat(p);rec['changed_during_read']=st.st_size!=en.st_size or st.st_mtime_ns!=en.st_mtime_ns
        rec['availability']='read'
    except Exception as e:rec.update(availability='read_error',error=str(e))
    rec['read_finished_at']=now();target=OUT/(sid+'.json');target.write_text(json.dumps(clean(rec),ensure_ascii=False,separators=(',',':')))
    rec['export_path']=str(target.relative_to(BASE));print(kind,pathlib.Path(p).name,rec['availability'],flush=True)
    return rec
if __name__=='__main__':
    import sys
    plan=json.loads((BASE/sys.argv[1]).read_text())
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(read,plan))
    (BASE/'audit'/sys.argv[2]).write_text(json.dumps(clean(records),ensure_ascii=False,indent=2))
