"""Read existing comparison ROOT files; never modify analysis inputs."""
from pathlib import Path
import csv
import ROOT

ROOT.gROOT.SetBatch(True)
BASE = Path('/home/murata_t/ScECAL_BeamTest/analysis/result/cog_channel_scan')
OUT = Path(__file__).parent
rows = []
layers = []
for energy in [30, 50, 100, 150, 250]:
    for cog in [0, 5, 10, 20]:
        for radius in ['22.5', '68.2']:
            condition = f'cog{cog}mm_channel{radius}mm'
            pair = {}
            for sample, subdir in [('data', 'data'), ('mc', 'simulation/threshold')]:
                path = BASE / subdir / condition / f'{energy}GeV_{condition}.root'
                if not path.exists():
                    continue
                f = ROOT.TFile.Open(str(path), 'READ')
                d = f.Get('beforeEventCut')
                h = d.Get(f'edep_{energy}GeV')
                n = d.Get(f'nhits_{energy}GeV')
                ef = d.Get('fit_gaus')
                nf = d.Get('fit_gaus_nhit')
                # Reset axis ranges so full stored histogram statistics are used.
                h.GetXaxis().SetRange(0, 0)
                n.GetXaxis().SetRange(0, 0)
                row = dict(energy=energy, cog=cog, channel=radius, sample=sample,
                           entries=h.GetEntries(), e_fit=ef.GetParameter(1),
                           e_fit_error=ef.GetParError(1), e_sigma=ef.GetParameter(2),
                           e_fit_chi2=ef.GetChisquare(), e_fit_ndf=ef.GetNDF(),
                           e_fit_min=ef.GetXmin(), e_fit_max=ef.GetXmax(),
                           n_fit=nf.GetParameter(1), n_fit_error=nf.GetParError(1),
                           e_mean=h.GetMean(), n_mean=n.GetMean(),
                           e_overflow=h.GetBinContent(h.GetNbinsX()+1),
                           n_overflow=n.GetBinContent(n.GetNbinsX()+1), path=str(path))
                for typ in [10, 15]:
                    hh = d.Get(f'edep_1hit_{typ}um')
                    row[f'hit_overflow_{typ}um'] = hh.GetBinContent(hh.GetNbinsX()+1)/hh.GetEntries() if hh.GetEntries() else 0
                rows.append(row)
                pair[sample] = row
                if energy == 100:
                    for layer in range(30):
                        le = d.Get(f'edep_{energy}GeV_layer{layer}')
                        ln = d.Get(f'nhits_{energy}GeV_layer{layer}')
                        if not le or not ln:
                            continue
                        le.GetXaxis().SetRange(0,0)
                        ln.GetXaxis().SetRange(0,0)
                        layers.append(dict(cog=cog,channel=radius,sample=sample,layer=layer,
                                           e_mean=le.GetMean(),n_mean=ln.GetMean(),
                                           p_n8=ln.GetBinContent(ln.FindBin(8))/ln.GetEntries(),
                                           entries=ln.GetEntries()))
                f.Close()
            if len(pair) == 2:
                a,b=pair['data'],pair['mc']
                print(f'{energy:3} GeV cog={cog:2} channel={radius:4} '
                      f'R_E={a["e_fit"]/b["e_fit"]:.5f} R_N={a["n_fit"]/b["n_fit"]:.5f} '
                      f'R_Emean={a["e_mean"]/b["e_mean"]:.5f} '
                      f'events={a["entries"]:.0f}/{b["entries"]:.0f}')
for name, table in [('saved_fit_values.csv',rows),('saved_layer_means_100GeV.csv',layers)]:
    if table:
        with (OUT/name).open('w',newline='') as out:
            w=csv.DictWriter(out,fieldnames=table[0].keys());w.writeheader();w.writerows(table)

# Compare mean energy at the same central hit count, without a new fit.
conditional=[]
for sample,subdir in [('data','data'),('mc','simulation/threshold')]:
    path=BASE/subdir/'cog5mm_channel22.5mm/100GeV_cog5mm_channel22.5mm.root'
    f=ROOT.TFile.Open(str(path),'READ')
    h=f.Get('beforeEventCut/hit_vs_e_100GeV')
    assert h.GetYaxis().GetBinWidth(1)==1
    for low in [190,195,200,205,210]:
        high=low+4
        projection=h.ProjectionX(f'conditional_{sample}_{low}',h.GetYaxis().FindBin(low),h.GetYaxis().FindBin(high))
        conditional.append(dict(sample=sample,n_low=low,n_high=high,n_events=projection.Integral(),
                                e_mean_MeV=projection.GetMean(),e_mean_error_MeV=projection.GetMeanError()))
    f.Close()
with (OUT/'conditional_energy_100GeV.csv').open('w',newline='') as out:
    w=csv.DictWriter(out,fieldnames=conditional[0]);w.writeheader();w.writerows(conditional)
