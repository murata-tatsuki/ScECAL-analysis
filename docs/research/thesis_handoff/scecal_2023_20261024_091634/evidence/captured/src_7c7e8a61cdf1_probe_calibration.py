"""Bounded read-only probes of calibration trees and corresponding MC truth."""
from pathlib import Path
import csv, math, statistics, re, json
import ROOT
ROOT.gROOT.SetBatch(True)
OUT=Path(__file__).parent
SRC=Path('/home/murata_t/ScECAL_BeamTest/analysis/cog_channel_scan/SingleEnergyAnalysis.cc').read_text()
mask=set(map(int,re.findall(r'\d+',SRC.split('const unordered_set<int> bad_channels = {',1)[1].split('};',1)[0].split('// bad channel の CellID をここに列挙',1)[1])))
DATA=Path('/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/calib/e-/100GeV')
MC=Path('/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC')
rows=[]
for sample,p in [('data',p) for p in sorted(DATA.glob('*.root'))]+[('mc',MC/'calib/e-/sps/threshold/100GeV/e-_100GeV_0.root')]:
    f=ROOT.TFile.Open(str(p));t=f.Get('Calib_Hit')
    t.SetBranchStatus('*',0)
    for name in ['CellID','NewTemperature','Hit_HG_Energy','Hit_LG_Energy','Hit_Energy','Hit_X','Hit_Y']:
        t.SetBranchStatus(name,1)
    vals={};temps=[];passed=0;es=[];ns=[]
    for entry in range(min(2000,t.GetEntries())):
        t.GetEntry(entry)
        ids=list(t.CellID);ee=list(t.Hit_Energy);xx=list(t.Hit_X);yy=list(t.Hit_Y)
        hh=list(t.Hit_HG_Energy);ll=list(t.Hit_LG_Energy)
        tt=list(t.NewTemperature)
        sx={9:0.,10:0.};sy={9:0.,10:0.};sw={9:0.,10:0.};e_sum=0.;nh=0
        for c,e,x,y,h,l,temp in zip(ids,ee,xx,yy,hh,ll,tt):
            layer=c//100000
            if c in mask or layer>=30 or e<0:continue
            if layer in sw:
                sx[layer]+=e*x;sy[layer]+=e*y;sw[layer]+=e
            if abs(x)>=22.5 or abs(y)>=22.5:continue
            e_sum+=e;nh+=1
            temps.append(temp)
            if 4<=layer<=27 and 2<=h<5:
                vals.setdefault(c,[]).append(l-h)
        if all(sw[k]>0 and abs(sx[k]/sw[k])<5 and abs(sy[k]/sw[k])<5 for k in sw):
            passed+=1;es.append(e_sum);ns.append(nh)
    row=dict(sample=sample,path=str(p),entries_sampled=min(2000,t.GetEntries()),
             cog5_pass=passed,temp_median=statistics.median(temps),
             central_E_median_cog5=statistics.median(es),central_N_median_cog5=statistics.median(ns),
             central_HGLG_delta_median=statistics.median([v for a in vals.values() for v in a]))
    rows.append(row);print(row,flush=True)
    with (OUT/f'HGLG_overlap_{p.stem}.csv').open('w',newline='') as out:
        w=csv.writer(out);w.writerow(['CellID','n','LG_minus_HG_median_MeV'])
        w.writerows((c,len(a),statistics.median(a)) for c,a in sorted(vals.items()))
    f.Close()
with (OUT/'calibration_sample_100GeV.csv').open('w',newline='') as out:
    w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

fc=ROOT.TFile.Open(str(MC/'calib/e-/sps/threshold/100GeV/e-_100GeV_0.root'))
ft=ROOT.TFile.Open(str(MC/'generate/e-/sps/100GeV/e-_100GeV_0.root'))
tc=fc.Get('Calib_Hit');tt=ft.Get('MC_Truth')
print('MC truth/calib entries',tt.GetEntries(),tc.GetEntries(),flush=True)
assert tc.GetEntries()==tt.GetEntries()
tc.SetBranchStatus('*',0)
for name in ['CellID','Hit_Energy','Hit_HG_Energy','Hit_LG_Energy','Event_Num']:
    tc.SetBranchStatus(name,1)
tt.SetBranchStatus('*',0)
for name in ['CellID','Hit_Energy']:tt.SetBranchStatus(name,1)
edges=[.5,1,2,5,10,20,50,100,200,500]
acc={}
for entry in range(min(500,tc.GetEntries())):
    tc.GetEntry(entry);tt.GetEntry(entry)
    assert int(tc.Event_Num)==entry
    ids=list(tt.CellID);energies=list(tt.Hit_Energy)
    assert len(set(ids))==len(ids)
    truth=dict(zip(ids,energies))
    for c,e,h,l in zip(tc.CellID,tc.Hit_Energy,tc.Hit_HG_Energy,tc.Hit_LG_Energy):
        if c in mask or c//100000>=30:continue
        typ=10 if 4<=c//100000<=27 else 15
        q=truth[c]
        for a,b in zip(edges[:-1],edges[1:]):
            if a<=q<b:
                acc.setdefault((typ,a,b),[]).append((e/q,h/q,l/q));break
closure=[]
for (typ,a,b),v in sorted(acc.items()):
    row=dict(type_um=typ,truth_min_MeV=a,truth_max_MeV=b,n=len(v),
             selected_to_truth_mean=statistics.mean(x[0] for x in v),
             HG_to_truth_mean=statistics.mean(x[1] for x in v),
             LG_to_truth_mean=statistics.mean(x[2] for x in v))
    closure.append(row);print('MC closure',row,flush=True)
with (OUT/'mc_truth_calibration_closure_sample.csv').open('w',newline='') as out:
    w=csv.DictWriter(out,fieldnames=closure[0]);w.writeheader();w.writerows(closure)
