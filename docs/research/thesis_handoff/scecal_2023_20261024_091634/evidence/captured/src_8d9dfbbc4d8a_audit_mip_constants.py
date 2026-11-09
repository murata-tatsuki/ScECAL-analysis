"""Audit MIP fallback means and infer constants from saved MC raw/calib pairs."""
from pathlib import Path
import csv,json,statistics
import ROOT
ROOT.gROOT.SetBatch(True)
OUT=Path(__file__).parent
BASE=Path('/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new')
SHARE=BASE/'Analysis_edit/share'
f=ROOT.TFile.Open(str(SHARE/'all_auto_muon_v4_trackfit.root'));t=f.Get('MIP_Fit')
rows=[];mip={};quality={}
for e in t:
    c=int(e.CellID);l=c//100000;ndf=int(e.NDF);chi=float(e.ChiSquare)
    # Match C++ comparison: NaN > cut is false. No observed zero-NDF entry is
    # needed for the quoted finite means, but preserve its comparison semantics.
    q=chi/ndf if ndf else (float('nan') if chi==0 else float('inf') if chi>0 else -float('inf'))
    v=float(e.LandauMPV)
    rows.append(dict(cellid=c,layer=l,mpv=v,chi2_ndf=q,accepted_by_calib=(not q>(2 if 4<=l<=27 else 1.7))))
    if v>0:mip[c]=v;quality[c]=q
f.Close()
summary={}
for typ in [10,15]:
    a=[r for r in rows if r['layer']<30 and (4<=r['layer']<=27)==(typ==10) and r['accepted_by_calib']]
    summary[str(typ)]=dict(accepted=len(a),nonpositive=sum(r['mpv']<=0 for r in a),
                          fallback_mean=statistics.mean(r['mpv'] for r in a),
                          positive_only_mean=statistics.mean(r['mpv'] for r in a if r['mpv']>0))
valid={int(r['CellID']) for r in csv.DictReader((OUT/'HGLG_overlap_e-_100GeV_0.csv').open())}
summary['central_10um_observed']=dict(n=len(valid),rejected_by_quality=sum(quality[c]>2 for c in valid))
with (OUT/'mip_parameter_audit.csv').open('w',newline='') as out:
    w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
f=ROOT.TFile.Open(str(SHARE/'pedestal2023_SPS.root'));t=f.Get('ChnLevel');ped={}
for e in t:ped.update(zip(list(e.CellID),list(e.PedHighMean)))
f.Close()
fr=ROOT.TFile.Open(str(BASE/'Result_MC/decode/e-/sps/threshold/100GeV/e-_100GeV_0.root'));tr=fr.Get('Raw_Hit')
fc=ROOT.TFile.Open(str(BASE/'Result_MC/calib/e-/sps/threshold/100GeV/e-_100GeV_0.root'));tc=fc.Get('Calib_Hit')
ft=ROOT.TFile.Open(str(BASE/'Result_MC/generate/e-/sps/100GeV/e-_100GeV_0.root'));tt=ft.Get('MC_Truth')
for tree,names in [(tr,['CellID','HG_Charge','TriggerID']),(tc,['CellID','Hit_HG_Energy','Event_Num','NewTemperature']),(tt,['CellID','Hit_Energy'])]:
    tree.SetBranchStatus('*',0)
    for n in names:tree.SetBranchStatus(n,1)
digi=[];ratios=[];samples={}
for i in range(50):
    tr.GetEntry(i);tc.GetEntry(i);tt.GetEntry(i)
    assert tr.TriggerID==tc.Event_Num==i
    assert all(abs(x-20)<1e-6 for x in tc.NewTemperature)
    truth=dict(zip(tt.CellID,tt.Hit_Energy));cal=dict(zip(tc.CellID,tc.Hit_HG_Energy))
    for c,q in zip(tr.CellID,tr.HG_Charge):
        if c not in valid or c not in cal or c not in ped or c not in mip or c not in truth:continue
        if 1<truth[c]<3 and cal[c]>0 and q<2000:
            eff=(q-ped[c])*0.305/cal[c]
            digi.append((q-ped[c])/(truth[c]/.305*mip[c]));ratios.append(mip[c]/eff)
            samples.setdefault(c,[]).append(eff)
summary['raw_calib_probe']=dict(events=50,hits=len(digi),
                               median_raw_over_linear=statistics.median(digi),
                               median_file_mip_over_effective=statistics.median(ratios))
with (OUT/'mip_inferred_from_saved_mc.csv').open('w',newline='') as out:
    w=csv.writer(out);w.writerow(['cellid','mip_file_ADC','inferred_mip_median_ADC','chi2_ndf','n_hits'])
    for c,values in sorted(samples.items()):w.writerow([c,mip[c],statistics.median(values),quality[c],len(values)])
(OUT/'mip_audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
