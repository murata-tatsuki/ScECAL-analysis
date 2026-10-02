"""Numerical edge cases that could change the physics conclusion."""
import math
from pathlib import Path
import tempfile
import unittest
import subprocess
from array import array
import numpy as np
from mip_tools import policy_audit, common_weights, load_fit_function, fit_spectrum
from run import energy_audit


def constant(c, value, chi=1, ndf=10):
    return dict(cellid=c, mpv=value, chi2=chi, ndf=ndf, width=10, sigma=5)


class PhysicsChecks(unittest.TestCase):
    def test_legacy_mean_and_clean_fallback_and_runtime_override(self):
        rows, means = policy_audit([constant(0, 400), constant(400000, 100),
                                   constant(400001, -10), constant(400002, 120, 30, 1),
                                   constant(400003, 0), constant(400004, 1)])
        by_id = {r['cellid']: r for r in rows}
        by_policy = {(r['pitch_um'], r['policy']): r for r in means}
        self.assertEqual(by_policy[10, 'legacy']['mean_adc'], 91/4)
        self.assertEqual(by_policy[10, 'clean']['mean_adc'], 101/2)
        self.assertEqual(by_id[400001]['legacy_adc'], -10)
        self.assertEqual(by_id[400001]['clean_adc'], 50.5)
        self.assertEqual(by_id[400002]['digi_adc'], 120)  # digi has no chi2 cut
        self.assertEqual(by_id[400002]['legacy_source'], 'fallback')
        self.assertEqual(by_id[400003]['digi_adc'], 90)   # runtime <=5 override
        self.assertEqual(by_id[400004]['digi_source'], 'fallback')  # sentinel 1

    def test_nonfinite_and_nonpositive_ndf_do_not_poison_clean_pool(self):
        rows, means = policy_audit([constant(0, 400), constant(400000, 100),
                                   constant(400001, 150, 0, 0), constant(400002, 160, 1, -1),
                                   constant(400003, math.nan), constant(400004, 170, math.inf)])
        by_id = {r['cellid']: r for r in rows}
        self.assertEqual(by_id[400001]['clean_adc'], 100)
        self.assertIn('nonpositive_ndf', by_id[400002]['clean_reason'])
        self.assertIn('nonfinite', by_id[400003]['clean_reason'])
        self.assertTrue(math.isnan(by_id[400005]['legacy_adc']))
        self.assertEqual(by_id[400005]['clean_adc'], 100)

    def test_overlap_is_symmetric_and_does_not_match_disjoint_runs(self):
        a = dict(sample=np.array([0]*4+[1]*8+[0, 1]), run=np.array([29]*12+[30, 31]),
                 cellid=np.repeat(400000, 14), longitudinal=np.zeros(14),
                 transverse=np.zeros(14), sx=np.zeros(14), sy=np.zeros(14))
        w, rows = common_weights(a, min_bin=2)
        self.assertEqual(w[:4].sum(), 4)
        self.assertEqual(w[4:12].sum(), 4)
        self.assertEqual(w[-2:].sum(), 0)
        np.testing.assert_array_equal(w[4:12], np.repeat(.5, 8))

    def test_total_energy_freezes_hits_and_invalidates_whole_bad_event(self):
        # First event: negative legacy MPV is retained in baseline, fixed in clean.
        # Second event: raw join missing; a finite subtotal must not become a total.
        p = [dict(cellid=400000, legacy_adc=100., clean_adc=100.),
             dict(cellid=400001, legacy_adc=-10., clean_adc=100.)]
        for row in p:
            row.update(legacy_source='individual', clean_source='individual' if row['cellid']==400000 else 'fallback')
        h = dict(sample=np.array([0, 0, 1]), cellid=np.array([400000, 400001, 400000]),
                 temperature=np.repeat(20., 3), raw_hg=np.array([200., 200., math.nan]),
                 saved_hg=np.array([.305, -3.05, .305]), saved_energy=np.array([.305, -3.05, .305]),
                 raw_ok=np.array([1, 1, 0]), event_index=np.array([0, 0, 1]))
        e = dict(event_index=np.array([0, 1]), sample=np.array([0, 1]), run=np.array([29, 29]),
                 event=np.array([4, 9]), event_time=np.array([5, 9]), n_calib=np.array([2, 1]))
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            energy_audit(None, out, h, e, [p, p], [{400000:100., 400001:100.}]*2)
            import csv
            with (out/'event_energy.tsv').open() as f:
                rows = list(csv.DictReader(f, delimiter='\t'))
            self.assertEqual(rows[0]['valid'], 'True')
            self.assertAlmostEqual(float(rows[0]['legacy_energy_mev']), -2.745)
            self.assertAlmostEqual(float(rows[0]['clean_energy_mev']), .61)
            self.assertEqual(rows[0]['n_hits'], '2')
            self.assertEqual(rows[1]['valid'], 'False')
            self.assertTrue(math.isnan(float(rows[1]['clean_energy_mev'])))

    def test_fitter_recovers_planted_ten_percent_excess(self):
        import ROOT
        ROOT.gROOT.SetBatch(True)
        load_fit_function(ROOT)
        rng = ROOT.TRandom3(513)
        # ROOT Landau location differs from MPV by -0.22278 * width.
        vals = np.array([rng.Landau(100+.22278298*10, 10)+rng.Gaus(0, 15) for _ in range(12000)])
        with tempfile.TemporaryDirectory() as directory:
            f = ROOT.TFile.Open(str(Path(directory)/'fits.root'), 'CREATE')
            data = fit_spectrum(ROOT, f, 'data', vals, np.ones(len(vals)), 100, 30, 450)
            mc = fit_spectrum(ROOT, f, 'mc', 1.1*vals, np.ones(len(vals)), 100, 30, 450)
            f.Close()
        self.assertTrue(data['valid'], data)
        self.assertTrue(mc['valid'], mc)
        self.assertAlmostEqual(data['mpv'], 100, delta=4)
        self.assertAlmostEqual(mc['mpv']/data['mpv'], 1.1, delta=.03)

    def test_collector_joins_shuffled_events_and_rejects_ambiguous_raw_hits(self):
        import ROOT
        ROOT.gROOT.SetBatch(True)
        def make_file(path, name, event_branch, records):
            f = ROOT.TFile.Open(str(path), 'RECREATE')
            tree = ROOT.TTree(name, name)
            event, time = array('i', [0]), array('i', [0])
            tree.Branch(event_branch, event, event_branch+'/I')
            tree.Branch('Event_Time', time, 'Event_Time/I')
            vectors = {}
            length = array('d', [2.])
            if name == 'Raw_Hit':
                specs = {'CellID':'int', 'HitTag':'int', 'HG_Charge':'double'}
            elif name == 'Calib_Hit':
                specs = {k:'double' for k in ('Hit_Energy', 'Hit_HG_Energy', 'NewTemperature', 'Hit_X', 'Hit_Y', 'Hit_Z')}
                specs['CellID'] = 'int'
            else:
                specs = {'hitCellnew':'int', 'trackFitPars':'double'}
                tree.Branch('realLength', length, 'realLength/D')
            for key, typ in specs.items():
                vectors[key] = ROOT.std.vector(typ)()
                tree.Branch(key, vectors[key])
            for event_id, values in records:
                event[0] = event_id; time[0] = event_id+1000
                for key, vector in vectors.items():
                    vector.clear()
                    for value in values[key]:
                        vector.push_back(value)
                tree.Fill()
            tree.Write();f.Close()
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            make_file(d/'raw.root', 'Raw_Hit', 'TriggerID', [
                (101, dict(CellID=[400000,400000], HitTag=[1,1], HG_Charge=[202.,999.])),
                (100, dict(CellID=[400000], HitTag=[1], HG_Charge=[101.]))])
            make_file(d/'cal.root', 'Calib_Hit', 'Event_Num', [
                (i, dict(CellID=[400000], Hit_Energy=[1.], Hit_HG_Energy=[1.], NewTemperature=[20.],
                         Hit_X=[0.], Hit_Y=[0.], Hit_Z=[40.8])) for i in (100,101,102)])
            make_file(d/'track.root', 'T_Event', 'TriggerID', [
                (i, dict(hitCellnew=[l*100000 for l in range(14)],
                         trackFitPars=[0.,.001,0.,.1,0.,.001,0.,.1,1.,10.])) for i in (101,100,102)])
            manifest = d/'manifest.tsv'
            manifest.write_text('\t'.join(['0','29',str(d/'raw.root'),str(d/'cal.root'),str(d/'track.root'),'-'])+'\n')
            subprocess.run([str(Path(__file__).parent/'CollectMuon'),str(manifest),str(d/'flat.root'),
                            str(d/'validation.tsv'),'0','.05','3','0'], check=True, capture_output=True)
            f = ROOT.TFile.Open(str(d/'flat.root'));t = f.Get('Hits')
            found = {}
            for entry in t:
                found[int(entry.event)] = (int(entry.raw_ok), float(entry.raw_hg), int(entry.mip_selected))
            self.assertEqual(found[100], (1,101.,1))
            self.assertEqual(found[101][0], 0)
            self.assertTrue(math.isnan(found[101][1]))
            self.assertNotIn(102, found)
            f.Close()


if __name__ == '__main__':
    unittest.main()
