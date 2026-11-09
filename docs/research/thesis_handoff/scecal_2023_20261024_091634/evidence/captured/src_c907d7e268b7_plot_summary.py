from pathlib import Path
from array import array
import csv, math
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
P=Path(__file__).parent
read=lambda n:list(csv.DictReader((P/n).open()))
fits=read('saved_fit_values.csv');layers=read('saved_layer_means_100GeV.csv');conditional=read('conditional_energy_100GeV.csv')
c=ROOT.TCanvas('reassessment','ScECAL data / MC reassessment',1500,500)
c.Divide(3,1)
keep=[]
def graph(x,y,err,color,marker):
    g=ROOT.TGraphErrors(len(x),array('d',x),array('d',y),array('d',[0.]*len(x)),array('d',err))
    g.SetMarkerStyle(marker);g.SetMarkerSize(1);g.SetMarkerColor(color);g.SetLineColor(color);g.SetLineWidth(2)
    keep.append(g);return g
def frame(xmin,ymin,xmax,ymax,title,xlabel,ylabel):
    ROOT.gPad.SetLeftMargin(.15);ROOT.gPad.SetBottomMargin(.14);ROOT.gPad.SetTopMargin(.12)
    h=ROOT.gPad.DrawFrame(xmin,ymin,xmax,ymax,f'{title};{xlabel};{ylabel}')
    h.GetXaxis().SetTitleSize(.043);h.GetYaxis().SetTitleSize(.043);h.GetYaxis().SetTitleOffset(1.6)
    keep.append(h)
def unity(x1,x2):
    l=ROOT.TLine(x1,1,x2,1);l.SetLineStyle(2);l.SetLineColor(ROOT.kGray+2);l.Draw();keep.append(l)
c.cd(1);frame(15,.99,265,1.18,'Central channel range: 22.5 mm','Beam energy [GeV]','Energy peak: data / MC');unity(15,265)
leg=ROOT.TLegend(.18,.72,.86,.86);leg.SetBorderSize(0);leg.SetFillStyle(0)
for cog,col,mar,label in [('0',ROOT.kBlue+1,20,'No CoG cut'),('5',ROOT.kRed+1,21,'CoG < 5 mm')]:
    x=[];y=[];err=[]
    for en in ['30','50','100','150','250']:
        pair={r['sample']:r for r in fits if r['energy']==en and r['cog']==cog and r['channel']=='22.5'}
        a,b=pair['data'],pair['mc'];v=float(a['e_fit'])/float(b['e_fit'])
        x.append(float(en));y.append(v);err.append(v*math.hypot(float(a['e_fit_error'])/float(a['e_fit']),float(b['e_fit_error'])/float(b['e_fit'])))
    g=graph(x,y,err,col,mar);g.Draw('PL SAME');leg.AddEntry(g,label,'lp')
leg.Draw();keep.append(leg)
c.cd(2);frame(-.5,.65,29.5,1.5,'100 GeV: CoG < 5 mm, channel < 22.5 mm','Layer (zero based)','Ratio of means: data / MC');unity(-.5,29.5)
leg=ROOT.TLegend(.19,.74,.79,.87);leg.SetBorderSize(0);leg.SetFillStyle(0)
for key,col,mar,label in [('e_mean',ROOT.kRed+1,20,'Layer energy'),('n_mean',ROOT.kBlue+1,24,'Layer hit count')]:
    vals=[]
    for l in range(30):
        a={r['sample']:r for r in layers if r['cog']=='5' and r['channel']=='22.5' and r['layer']==str(l)}
        vals.append(float(a['data'][key])/float(a['mc'][key]))
    g=graph(list(range(30)),vals,[0.]*30,col,mar);g.Draw('PL SAME');leg.AddEntry(g,label,'lp')
leg.Draw();keep.append(leg)
c.cd(3);frame(187,1.10,217,1.17,'100 GeV: CoG < 5 mm, channel < 22.5 mm','N_{hit} (groups of 5)','Mean energy at fixed N_{hit}: data / MC')
x=[];y=[];err=[]
for a in ['190','195','200','205','210']:
    p={r['sample']:r for r in conditional if r['n_low']==a};d,m=p['data'],p['mc']
    v=float(d['e_mean_MeV'])/float(m['e_mean_MeV']);x.append(int(a)+2);y.append(v)
    err.append(v*math.hypot(float(d['e_mean_error_MeV'])/float(d['e_mean_MeV']),float(m['e_mean_error_MeV'])/float(m['e_mean_MeV'])))
g=graph(x,y,err,ROOT.kBlack,20);g.Draw('P SAME')
c.SaveAs(str(P/'reassessment_summary.png'));c.SaveAs(str(P/'reassessment_summary.pdf'))
