#!/usr/bin/env python3
"""Build a PS manifest; missing counterparts are errors, never silently skipped."""
import argparse
import csv
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ps-base',type=Path,default=Path('/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps'))
    p.add_argument('--energies',default='0.5,1,2,3,4,5',help='electron energies in GeV, or none')
    p.add_argument('--muon-energy',default='10',help='muon energy in GeV, or none')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();rows=[]
    for particle,energy,kind in ([('mu-',a.muon_energy,'mip')] if a.muon_energy!='none' else [])+[
            ('e-',e,'response') for e in a.energies.split(',') if e!='none']:
        relative=Path(particle)/(energy+'GeV')
        files=sorted((a.ps_base/'decode'/relative).glob('ECAL_Run*.root'))
        if not files:p.error('No raw files: '+str(relative))
        for raw in files:
            cal=a.ps_base/'calib'/relative/raw.name
            track=a.ps_base/'trackFit3D'/relative/raw.name
            if not cal.is_file() or (kind=='mip' and not track.is_file()):
                p.error('Missing calibrated/track counterpart for '+str(raw))
            rows.append(dict(kind=kind,raw=str(raw.resolve()),calib=str(cal.resolve()),track=str(track.resolve()) if kind=='mip' else '-'))
    if not rows:p.error('No samples selected')
    if len({r['raw'] for r in rows})!=len(rows):p.error('Duplicate input energy/file')
    with a.output.open('x') as f:
        w=csv.DictWriter(f,fieldnames=['kind','raw','calib','track'],delimiter='\t');w.writeheader();w.writerows(rows)
    print(f'Saved {len(rows)} input files to {a.output}')


if __name__=='__main__':main()
