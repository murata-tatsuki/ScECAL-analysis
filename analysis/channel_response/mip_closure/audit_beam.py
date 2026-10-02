#!/usr/bin/env python3
"""Check all muon primary directions and the saved run-specific position maps."""
import argparse
from pathlib import Path
import re
import numpy as np
from mip_tools import write_tsv
from run import MC, fingerprint


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True, help='new directory')
    p.add_argument('--truth-dir', type=Path, default=MC/'Result_MC/generate/mu-/mu_track/100GeV')
    p.add_argument('--macro-dir', type=Path, default=MC/'CEPCScECAL/generate/mac/mu-/100GeV')
    a = p.parse_args()
    import ROOT
    ROOT.gROOT.SetBatch(True)
    a.output.mkdir(parents=True, exist_ok=False)
    result = []
    for path in sorted(a.truth_dir.glob('ECAL_Run*.root')):
        run = int(re.search(r'Run(\d+)', path.name).group(1))
        f = ROOT.TFile.Open(str(path));t = f.Get('MC_Truth')
        if not t or f.TestBit(ROOT.TFile.kRecovered):
            raise ValueError(f'Invalid/recovered truth: {path}')
        rdf = ROOT.RDataFrame(t)
        for i, key in enumerate(('dx', 'dy', 'dz')):
            rdf = rdf.Define(key, f'PrimaryDirection.size()==3 ? PrimaryDirection[{i}] : std::numeric_limits<double>::quiet_NaN()')
        arr = rdf.AsNumpy(['dx', 'dy', 'dz', 'PrimaryPosX', 'PrimaryPosY', 'PrimaryPosZ', 'PrimaryPDG'])
        n = len(arr['dx'])
        if not n or not np.all(np.isfinite(np.column_stack([arr[k] for k in ('dx', 'dy', 'dz')]))):
            raise ValueError(f'Invalid directions in {path}')
        row = dict(run=run, n_events=n, n_nonmuon=int(np.count_nonzero(np.abs(arr['PrimaryPDG']) != 13)),
                   n_non_axial=int(np.count_nonzero((np.abs(arr['dx'])>1e-12)|(np.abs(arr['dy'])>1e-12)|(np.abs(arr['dz']-1)>1e-12))),
                   truth_path=str(path), truth_size=path.stat().st_size, truth_mtime_ns=path.stat().st_mtime_ns)
        for name in ('dx', 'dy', 'dz', 'PrimaryPosX', 'PrimaryPosY', 'PrimaryPosZ'):
            v = arr[name]
            row.update({name+'_mean': float(v.mean()), name+'_std': float(v.std()),
                        name+'_min': float(v.min()), name+'_max': float(v.max())})
        macro = a.macro_dir/f'ECAL_Run{run}.mac'
        text = macro.read_text() if macro.is_file() else ''
        match = re.search(r'^\s*/mydet/hitmapFile\s+(\S+)', text, re.M)
        row.update(macro=str(macro), macro_sha256=fingerprint(macro) if text else '',
                   macro_direction=' '.join(re.findall(r'^\s*/gps/direction\s+(.+)', text, re.M)),
                   position_map=match[1] if match else '', position_map_valid=False,
                   map_x_mean=np.nan, map_y_mean=np.nan, map_x_std=np.nan, map_y_std=np.nan)
        if text:
            (a.output/macro.name).write_text(text)
        if match and Path(match[1]).is_file():
            mf = ROOT.TFile.Open(match[1]);h = mf.Get('h_hitmap_smooth')
            if h:
                row.update(position_map_valid=True, map_x_mean=h.GetMean(1), map_y_mean=h.GetMean(2),
                           map_x_std=h.GetStdDev(1), map_y_std=h.GetStdDev(2))
            mf.Close()
        result.append(row);f.Close()
        print(f'Run{run}: {n} entries, non-axial={row["n_non_axial"]}', flush=True)
    if not result:
        raise ValueError('No truth files')
    write_tsv(a.output/'beam_audit.tsv', result)
    (a.output/'interpretation.txt').write_text(
        'Current generator_sps_mu.sh fixes /gps/direction 0 0 1.\n'
        'PrimaryGeneratorAction.cc draws (x,y) from h_hitmap_smooth and places it at z=-60 mm.\n'
        'beamsize_track.cc fills that map from reconstructed intercepts at z=0, with uniform pitch smearing.\n'
        'This matches a marginal position distribution, not event-by-event initial hits or position-angle correlations.\n'
        'At exactly axial incidence the transverse position is unchanged between the source and z=0.\n'
        'The table checks saved truth; current macros and position-map moments do not prove historical file identity.\n')


if __name__ == '__main__':
    main()
