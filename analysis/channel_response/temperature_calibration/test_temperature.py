"""Synthetic ROOT integration tests and numerical identifiability checks."""
from array import array
import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
from temperature_tools import factor, energies, slope

HERE=Path(__file__).resolve().parent


def table(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))


class Numerics(unittest.TestCase):
    def test_invalid_constants_are_audited_without_losing_valid_channels(self):
        from run import process, pedestal_row
        params={c:SimpleNamespace(pedPresent=True,mip=100.,ph=500.,pl=400.,gain=40.,offset=12.,mipFallback=False)
                for c in range(420008,420013)}
        params[420009].mip=-10
        params[420010].mip=0
        params[420011].pedPresent=False
        params[420012].mip=math.nan
        cells=np.repeat(list(params),3)
        h=dict(cellid=cells,run=np.full(len(cells),10),temperature=np.full(len(cells),29.),
               temperature_valid=np.ones(len(cells)),sensor_temperature=np.full(len(cells),29.),
               raw_hg=np.full(len(cells),1000.),raw_lg=np.full(len(cells),420.),
               saved_hg=np.full(len(cells),1.),saved_lg=np.full(len(cells),1.1),path_factor=np.ones(len(cells)))
        a=SimpleNamespace(reference_temperature=20.,coefficient_10=1.6/135,coefficient_15=3.5/230,
                          temp_width=.5,hg_max=2200.,compare_models=False,plateau_sigmas=3.,
                          fit_high_mip=3.,min_fit_hits=20,overlap_min=800.,overlap_max=2200.,adc_width=100.,min_hits=3)
        fit=dict(n=3,mpv=100.,mpv_err=2.,valid=True,reason='ok')
        for role in ['mip','response']:
            with patch('run.fit_spectrum',return_value=fit),np.errstate(divide='raise',invalid='raise'):
                mm,rr,cc,ee=process(None,a,0,role,h,params,{c:(20.,1.,'individual') for c in params},None)
            self.assertEqual([r['cellid'] for r in cc],[420008])
            self.assertEqual(sum(r['n'] for r in ee),12)
            reasons={r['cellid']:r['reason'] for r in ee}
            self.assertEqual(reasons,{420009:'nonpositive_mip',420010:'nonpositive_mip',
                                      420011:'missing_pedestal',420012:'nonfinite_mip'})
            self.assertEqual(len(mm),3 if role=='mip' else 0)
            self.assertEqual(len(rr),1 if role=='response' else 0)
        ped=dict(run='10',cellid='420009',memory='0',temp_bin='58',temperature='29',n='40',
                 hg_mean='502',lg_mean='401',hg_sem='.2',lg_sem='.1',missing_temperature='0')
        q=pedestal_row(ped,0,'mip',params[420009],3)
        self.assertEqual(q['status'],'ok');self.assertEqual(q['hg_residual'],2)
        q=pedestal_row(ped,0,'mip',params[420011],3)
        self.assertEqual(q['status'],'missing_pedestal');self.assertTrue(math.isnan(q['hg_residual']))
        self.assertEqual(q['hg_mean'],502)

    def test_sign_reference_and_order(self):
        self.assertAlmostEqual(float(factor(20,420008)),1)
        self.assertAlmostEqual(float(factor(29,420008)),1-9*1.6/135)
        self.assertAlmostEqual(float(factor(29,0)),1-9*3.5/230)
        e=energies(600,403,29,420008,500,400,100,40,12)
        f=float(factor(29,420008))
        self.assertAlmostEqual(e['lg_offset_before']-e['lg'],12*(f-1)*.305/100)

    def test_run_confounding(self):
        rows=[dict(run=r,temperature=t,y=5*r+.2*(t-20),err=.1) for r in [10,11] for t in [20,25,29]]
        fit=slope(rows,'y','err',fixed_run=True)
        self.assertEqual(fit['status'],'ok');self.assertAlmostEqual(fit['slope'],.2,places=10)
        rows=[dict(run=r,temperature=20+r,y=r,err=.1) for r in range(4)]
        self.assertEqual(slope(rows,'y','err',fixed_run=True)['status'],'run_temperature_confounded')
        self.assertEqual(slope(rows[:2],'y','err')['status'],'insufficient_points')


def constants(ROOT,base,mip=100):
    f=ROOT.TFile(str(base/'ped.root'),'RECREATE');t=ROOT.TTree('ChnLevel','')
    v={k:ROOT.std.vector(ty)() for k,ty in [('CellID','int'),('PedHighMean','double'),('PedLowMean','double')]}
    for k,a in v.items():t.Branch(k,a)
    for c in [0,420008]:
        for k,val in [('CellID',c),('PedHighMean',500),('PedLowMean',400)]:v[k].push_back(val)
    t.Fill();t.Write();f.Close()
    for file,tree,spec,rows in [
        ('mip.root','MIP_Fit',[('CellID','i'),('LandauMPV','d'),('ChiSquare','d'),('NDF','i')],[(0,100,1,10),(420008,mip,1,10)]),
        ('hl.root','InterCalib',[('CellID','i'),('Slope','d'),('Intercept','d')],[(0,.025,-.3),(420008,.025,-.3)]),
        ('threshold.root','expErfThre',[('CellID','i'),('threshold','d'),('sigma','d')],[(420008,10,1)])]:
        f=ROOT.TFile(str(base/file),'RECREATE');t=ROOT.TTree(tree,'')
        v={k:array(ty,[0]) for k,ty in spec}
        for k,ty in spec:t.Branch(k,v[k],f'{k}/{"I" if ty=="i" else "D"}')
        for row in rows:
            for (k,_),val in zip(spec,row):v[k][0]=val
            t.Fill()
        t.Write();f.Close()


def events(ROOT,base,name,rows,tracks=False,missing_sensor=False,duplicate=False):
    for mode in ['raw','cal']+(['track'] if tracks else []):
        f=ROOT.TFile(str(base/f'{name}_{mode}.root'),'RECREATE')
        t=ROOT.TTree(dict(raw='Raw_Hit',cal='Calib_Hit',track='T_Event')[mode],'')
        scalar=['Event_Time','Event_Num' if mode=='cal' else 'TriggerID']+(['Run_Num'] if mode!='track' else [])
        s={k:array('i',[0]) for k in scalar}
        for k,v in s.items():t.Branch(k,v,k+'/I')
        spec={'CellID':'int'}
        if mode=='raw':spec.update(HitTag='int',HG_Charge='double',LG_Charge='double')
        if mode=='cal':spec.update({k:'double' for k in ['Hit_HG_Energy','Hit_LG_Energy','NewTemperature','Hit_X','Hit_Y','Hit_Z']})
        if mode=='track':spec={'hitCellnew':'int','trackFitPars':'double'}
        v={k:ROOT.std.vector(ty)() for k,ty in spec.items()}
        for k,x in v.items():t.Branch(k,x)
        sensor=ROOT.std.vector('vector<double>')()
        if mode=='raw':t.Branch('Temperature',sensor)
        for i,run,temp,signal in (rows[::-1] if mode=='cal' else rows):
            for k in scalar:s[k][0]=run if k=='Run_Num' else (1000+i if k=='Event_Time' else i)
            for x in v.values():x.clear()
            fac=float(factor(temp,420008))
            h=500+signal/fac
            residual=.1+.02*(temp-20)+.03*(run-10)
            eh=signal*.305/100;el=eh+residual
            l=400+(el*100/.305-12)/(fac*40)
            if mode=='raw':
                values={'CellID':[420008,420108,100],'HitTag':[1,0,0],
                        'HG_Charge':[h,500+.5*(temp-20)+.2*(run-10)+(i%7-3),500+(i%7-3)],
                        'LG_Charge':[l,400+.2*(temp-20)+.1*(run-10)+(i%5-2),400+(i%5-2)]}
                if duplicate:
                    for x in values.values():x.append(x[0])
                sensor.clear()
                for layer in range(32):
                    ts=ROOT.std.vector('double')()
                    if not missing_sensor:
                        for _ in range(16):ts.push_back(temp)
                    sensor.push_back(ts)
            elif mode=='cal':
                values={'CellID':[420008],'Hit_HG_Energy':[eh],'Hit_LG_Energy':[el],
                        'NewTemperature':[temp],'Hit_X':[0.],'Hit_Y':[0.],'Hit_Z':[0.]}
            else:
                values={'hitCellnew':[-1]+[j*100000 for j in range(14)]+[420008],
                        'trackFitPars':[0.,0.,0.,0.,0.,0.,0.,0.,1.,10.]}
            for k,vals in values.items():
                for val in vals:v[k].push_back(val)
            t.Fill()
        t.Write();f.Close()


class EndToEnd(unittest.TestCase):
    def test_negative_legacy_mip_completes_with_explicit_exclusions(self):
        import ROOT,json
        ROOT.gROOT.SetBatch(True)
        with tempfile.TemporaryDirectory(prefix='scecal_negative_mip_') as folder:
            base=Path(folder);constants(ROOT,base,mip=-10)
            events(ROOT,base,'negative',[(i,10,29.,600.) for i in range(60)],tracks=True)
            (base/'manifest.tsv').write_text('kind\traw\tcalib\ttrack\n'
                'mip\tnegative_raw.root\tnegative_cal.root\tnegative_track.root\n')
            out=base/'out'
            cmd=[sys.executable,str(HERE/'run.py'),'--manifest',str(base/'manifest.tsv'),
                 '--output',str(out),'--no-plots','--channels','0,420008']
            for key in ['pedestal','mip','hl','threshold']:
                cmd.extend(['--'+key,str(base/('ped.root' if key=='pedestal' else key+'.root'))])
            result=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('WARNING: excluded',result.stdout)
            self.assertTrue((out/'COMPLETE').is_file())
            audit=next(r for r in table(out/'calibration_channels.tsv') if r['cellid']=='420008')
            self.assertEqual(float(audit['mip_adopted_ADC']),-10)
            excluded=table(out/'excluded_channels.tsv')
            self.assertEqual(sum(int(r['n']) for r in excluded),60)
            self.assertTrue(all(r['reason']=='nonpositive_mip' for r in excluded))
            self.assertEqual(table(out/'mip_fits.tsv'),[])
            self.assertEqual(table(out/'closure.tsv'),[])
            self.assertTrue(any(r['cellid']=='420008' and r['status']=='ok' for r in table(out/'pedestal.tsv')))
            summary=json.loads((out/'summary.json').read_text())
            self.assertEqual(summary['excluded_signal_channels'],1)
            self.assertEqual(summary['excluded_signal_hits'],60)
            self.assertEqual(summary['signal_status'],'no_valid_signal_constants')

    def test_pipeline(self):
        import ROOT
        ROOT.gROOT.SetBatch(True)
        with tempfile.TemporaryDirectory(prefix='scecal_temperature_test_') as folder:
            base=Path(folder);constants(ROOT,base)
            rng=ROOT.TRandom3(73091)
            rows=[]
            # Landau MPV=100, convolved Gaussian; raw response is divided by F(T).
            for run in [10,11]:
                for temp in [20.,25.,29.]:
                    for _ in range(2500):
                        signal=rng.Landau(100+.22278298*12,12)+rng.Gaus(0,10)
                        rows.append((len(rows),run,temp,signal))
            events(ROOT,base,'mip',rows,tracks=True)
            response=[(i,run,temp,float(500+5*(i%180))) for i,(run,temp) in enumerate([(r,t) for r in [10,11] for t in [20.,25.,29.] for _ in range(900)])]
            events(ROOT,base,'response',response)
            events(ROOT,base,'missing',[(i,12,29.,100.) for i in range(60)],missing_sensor=True)
            manifest=base/'manifest.tsv'
            manifest.write_text('kind\traw\tcalib\ttrack\n'+'\n'.join([
                f'mip\t{base}/mip_raw.root\t{base}/mip_cal.root\t{base}/mip_track.root',
                f'response\t{base}/response_raw.root\t{base}/response_cal.root\t-',
                f'response\t{base}/missing_raw.root\t{base}/missing_cal.root\t-'])+'\n')
            cmd=[sys.executable,str(HERE/'run.py'),'--manifest',str(manifest),'--output',str(base/'out'),
                 '--channels','0,420008','--min-hits','10','--min-fit-hits','200',
                 '--temp-width','1','--mip-measurement-temperature','29','--mip-constant-temperature','29']
            for key in ['pedestal','mip','hl','threshold']:
                cmd.extend(['--'+key,str(base/('ped.root' if key=='pedestal' else key+'.root'))])
            result=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            out=base/'out'
            self.assertTrue((out/'COMPLETE').exists())
            self.assertTrue((out/'figures/channel_420008.png').is_file())
            self.assertTrue((out/'figures/channel_0.png').is_file())
            png_timings=table(out/'png_timings.tsv')
            self.assertEqual({r['cellid'] for r in png_timings},{'0','420008'})
            for r in png_timings:
                png=out/f'figures/channel_{r["cellid"]}.png'
                self.assertTrue(png.read_bytes().startswith(b'\x89PNG\r\n\x1a\n'))
                self.assertEqual(png.stat().st_size,int(r['bytes']))
                self.assertGreater(float(r['raster_seconds']),0)
                self.assertGreater(float(r['encode_seconds']),0)
            self.assertEqual(table(out/'reference_audit.tsv')[0]['status'],'declared_mismatch')
            closure=table(out/'closure.tsv')
            self.assertTrue(all(float(r['production_hg_absmax'])<1e-10 and float(r['production_lg_absmax'])<1e-10 for r in closure))
            self.assertTrue(any(r['temperature_valid']=='False' for r in closure))
            self.assertTrue(all(r['matching']=='unique_event_key_or_raw_only' for r in table(out/'counters.tsv')))
            fits=table(out/'mip_fits.tsv')
            valid=[r for r in fits if r['valid']=='True']
            self.assertGreaterEqual(len(valid),15)
            corrected=[float(r['mpv_MIP']) for r in valid if r['stage']=='corrected']
            self.assertTrue(all(abs(x-1)<.06 for x in corrected),corrected)
            slopes=table(out/'temperature_slopes.tsv')
            raw=next(r for r in slopes if r['dataset']=='mip' and r['stage']=='raw' and r['scope']=='run_adjusted')
            corr=next(r for r in slopes if r['dataset']=='mip' and r['stage']=='corrected' and r['scope']=='run_adjusted')
            self.assertGreater(float(raw['slope']),.009)
            self.assertLess(abs(float(corr['slope'])),.004)
            residual=next(r for r in slopes if r['dataset']=='residual' and r['metric']=='saved_delta_mean' and r['scope']=='run_adjusted' and r['status']=='ok')
            self.assertAlmostEqual(float(residual['slope']),.02,places=9)
            pedestal=next(r for r in slopes if r['dataset']=='pedestal' and r['cellid']=='420008' and r['metric']=='hg_residual' and r['scope']=='run_adjusted' and r['status']=='ok')
            self.assertAlmostEqual(float(pedestal['slope']),.5,delta=.01)
            events(ROOT,base,'bad',[(0,10,20.,100.)],duplicate=True)
            bad=subprocess.run([str(HERE/'CollectTemperature'),'response',str(base/'bad_raw.root'),str(base/'bad_cal.root'),'-',str(base/'bad.root'),str(base/'bad.tsv'),str(base/'bad_count.tsv'),'420008','0','.5','.05'],capture_output=True,text=True)
            self.assertNotEqual(bad.returncode,0)
            self.assertIn('Duplicate full raw CellID',bad.stderr)
            # Existing output cannot be silently overwritten.
            again=subprocess.run(cmd,capture_output=True,text=True)
            self.assertNotEqual(again.returncode,0);self.assertIn('Output exists',again.stderr)


if __name__=='__main__':unittest.main()
