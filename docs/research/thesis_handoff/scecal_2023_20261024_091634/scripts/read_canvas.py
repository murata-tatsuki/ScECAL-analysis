#!/usr/bin/env python3
"""Read saved TCanvas/TGraph/TF1 primitives. Never Fit, Draw, or execute macros."""
import os,sys,json,pathlib,datetime,hashlib
import ROOT
ROOT.gROOT.SetBatch(True)
B=pathlib.Path(__file__).resolve().parents[1]
p='/home/murata_t/ScECAL_BeamTest/analysis/result/resolution/comparsion/comparison_5mm.root'
a={'source':p,'captured_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'root_version':ROOT.gROOT.GetVersion(),'objects':[],'operation':'read saved objects only; no Fit/Draw/new event processing'}
def fun(f):
 return {'class':f.ClassName(),'name':f.GetName(),'title':f.GetTitle(),'parameters':[{'index':i,'name':f.GetParName(i),'value':float(f.GetParameter(i)),'error':float(f.GetParError(i))} for i in range(f.GetNpar())],'chi2':float(f.GetChisquare()),'ndf':int(f.GetNDF()),'xmin':float(f.GetXmin()),'xmax':float(f.GetXmax()),'fit_status':'not available from TF1 alone'}
def walk(o,path,depth=0):
 if depth>8:return
 r={'path':path,'class':o.ClassName(),'name':o.GetName(),'title':o.GetTitle()}
 if o.InheritsFrom('TGraph'):
  r['points']=[]
  for i in range(o.GetN()):
   q={'index':i,'x':float(o.GetPointX(i)),'y':float(o.GetPointY(i))}
   if o.InheritsFrom('TGraphErrors'):q.update(ex=float(o.GetErrorX(i)),ey=float(o.GetErrorY(i)))
   if o.InheritsFrom('TGraphAsymmErrors'):q.update(ex_low=float(o.GetErrorXlow(i)),ex_high=float(o.GetErrorXhigh(i)),ey_low=float(o.GetErrorYlow(i)),ey_high=float(o.GetErrorYhigh(i)))
   r['points'].append(q)
  r['functions']=[fun(f) for f in o.GetListOfFunctions() if f.InheritsFrom('TF1')]
 elif o.InheritsFrom('TF1'):r.update(fun(o))
 elif o.InheritsFrom('TLegend'):
  r['labels']=[e.GetLabel() for e in o.GetListOfPrimitives()]
 elif o.InheritsFrom('TPaveText'):
  r['lines']=[e.GetTitle() for e in o.GetListOfLines()]
 a['objects'].append(r)
 if o.InheritsFrom('TPad'):
  for i,x in enumerate(o.GetListOfPrimitives()):
   if x.InheritsFrom('TExec'):continue
   walk(x,path+'/'+x.GetName()+'['+str(i)+']',depth+1)
st=os.stat(p);a.update(size=st.st_size,mtime_ns=st.st_mtime_ns)
f=ROOT.TFile.Open(p,'READ')
if not f or f.IsZombie():raise RuntimeError('Cannot open comparison ROOT')
for k in f.GetListOfKeys():
 o=k.ReadObj()
 if o.InheritsFrom('TExec'):continue
 walk(o,k.GetName())
f.Close();a['finished_at']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(B/'evidence/comparison_canvas.json').write_text(json.dumps(a,ensure_ascii=False,indent=2))
print('read',len(a['objects']),'objects',sum(len(x.get('points',[])) for x in a['objects']),'saved graph points',flush=True)
