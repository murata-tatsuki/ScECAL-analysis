#!/usr/bin/env python3
"""Eight data series: four CoG selections times tail cut on/off."""
from pathlib import Path
import math
import tempfile
import shlex
import ROOT
from test_sample_comparison import fixture, run, ENERGIES, legends, check_fit_curves
ROOT.gROOT.SetBatch(True)
SOURCE=Path(__file__).resolve().parents[1]

def metadata(path, flag):
    output=ROOT.TFile(str(path),'UPDATE')
    ROOT.TParameter('int')('hl_tail_excluded',flag).Write('hl_tail_excluded',ROOT.TObject.kOverwrite)
    ROOT.TNamed('hl_tail_status','complete').Write('hl_tail_status',ROOT.TObject.kOverwrite)
    output.Close()

def checks(work):
    labels=[]; groups=[]; flags=[]
    for cog in [200,20,10,5]:
        for mode,flag in [('exclude',1),('keep',0)]:
            index=len(groups)
            label=f'Data tail {mode}'
            files=[work/'data'/mode/f'{cog}mm'/f'{energy:g}GeV_{cog}mm.root' for energy in ENERGIES]
            for energy,path in zip(ENERGIES,files):
                fixture(path,energy,index);metadata(path,flag)
            groups.append(files);labels.append(label);flags.append(flag)
    output=work/'comparison.root';figures=work/'figures'
    config=work/'eight.sh'
    config.write_text('\n'.join([
        'multi_options=(--compare-tail-selections --skip-channel-plots)',
        *['multi_sample '+shlex.join([label,*map(str,reversed(files))]) for label,files in zip(labels,groups)],
        'multi_run '+shlex.join([str(output),str(figures)])]))
    dry=run(['bash',SOURCE/'run_samples.sh','--dry-run','multi',config])
    args=shlex.split(dry)
    assert args[2:4]==['8',str(len(ENERGIES))]
    assert '--compare-tail-selections' in args and '--exclude-tail-events' not in args
    run(['bash',SOURCE/'run_samples.sh','multi',config])
    file=ROOT.TFile(str(output))
    assert file.Get('hl_tail_excluded').GetVal()==-1
    for i,flag in enumerate(flags):assert file.Get(f'hl_tail_excluded_sample_{i}').GetVal()==flag
    colors=[ROOT.kAzure+1,ROOT.kViolet+1,ROOT.kGreen+2,ROOT.kRed+1]
    captions=[f'{label}, {"nocut" if cog==200 else str(cog)+" mm"}'
              for cog in [200,20,10,5] for label in ['Data tail exclude','Data tail keep']]
    for name in ['resolution','resolution_with_noise']:
      for log in [False,True]:
        key=name+('_loglog' if log else '')
        canvas=file.Get(key)
        check_fit_curves(canvas)
        assert canvas.GetLogx()==int(log) and canvas.GetLogy()==int(log)
        entries=legends(canvas)
        for legend in canvas.GetListOfPrimitives():
            if legend.InheritsFrom('TLegend'):
                assert legend.GetX1NDC() >= 1-canvas.GetRightMargin()
        assert all(any(entry.startswith(label) for entry in entries) for label in captions), entries
        assert not any('200 mm' in entry or 'no cut' in entry for entry in entries)
        graphs=[obj for obj in canvas.GetListOfPrimitives() if obj.InheritsFrom('TGraphErrors')]
        assert len(graphs)==8
        assert len({g.GetLineColor() for g in graphs})==4
        for index,g in enumerate(graphs):
            assert g.GetLineColor()==colors[index//2]
            assert g.GetLineStyle()==(ROOT.kSolid if index%2==0 else ROOT.kDotted)
            if log:
                assert 0<g.GetXaxis().GetXmin()<min(ENERGIES)
                assert g.GetXaxis().GetXmax()>max(ENERGIES)
                assert g.GetMinimum()>0
            else:
                assert g.GetXaxis().GetXmin()==0 and g.GetXaxis().GetXmax()==130
                assert g.GetMinimum()==0 and g.GetMaximum()==.3
            for point,energy in enumerate(ENERGIES):
                expected=math.sqrt((.18+index*.02)**2/energy+.01**2)
                assert abs(g.GetPointY(point)-expected)<1e-12
        assert (figures/(key+'.png')).exists()
    file.Close()
    print('PASS: eight series, CoG colors, solid/dotted conditions, nocut, linear/log axes, unchanged values and selection metadata',flush=True)
    args=[SOURCE/'MultiEnergyAnalysis',work/'invalid.root',8,len(ENERGIES),*[p for group in groups for p in group],work/'invalid_figs','--labels',*labels,'--skip-channel-plots']
    run(args,'tail selection mismatch')
    run(args+['--compare-tail-selections','--keep-tail-events'],'cannot be combined')
    metadata(groups[0][-1],0)
    run(args+['--compare-tail-selections'],'tail selection changes across energies')
    metadata(groups[0][-1],1)
    f=ROOT.TFile(str(groups[0][-1]),'UPDATE');f.Delete('hl_tail_excluded;*');f.Close()
    run(args+['--compare-tail-selections'],'missing hl_tail_excluded')
    run(['bash',SOURCE/'run_samples.sh','--keep-tail-events','multi',config],'cannot be combined')
    print('PASS: implicit mixing, conflicting options, mixed energy selections and unknown selection rejected',flush=True)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--work-dir',type=Path);options=parser.parse_args()
    if options.work_dir:
        options.work_dir.mkdir(parents=True,exist_ok=True);checks(options.work_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='tail-eight-') as directory:checks(Path(directory))
