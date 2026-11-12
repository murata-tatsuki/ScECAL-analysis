"""Channel masks from threshold fit quality; never change event or ADC selection."""
import csv
import hashlib
import json
from pathlib import Path
import re
import numpy as np
import uproot

# Orthogonal scans first; combined quality presets are available as follow-ups.
PRESETS = {
    'allChannel': None,
    'fitValid': {},
    'precision': dict(max_rel_threshold_error=.05, max_rel_sigma_error=.30,
                      min_sigma_adc=.6, max_sigma_adc=29.5),
    'chi2_20': dict(max_chi2ndf=20.),
    'chi2_10': dict(max_chi2ndf=10.),
    'chi2_5': dict(max_chi2ndf=5.),
    'chi2_3': dict(max_chi2ndf=3.),
    'chi2_2': dict(max_chi2ndf=2.),
    'quality_loose': dict(max_chi2ndf=10., max_rel_threshold_error=.10,
                         max_rel_sigma_error=.50, min_sigma_adc=.6, max_sigma_adc=29.5),
    'quality_medium': dict(max_chi2ndf=5., max_rel_threshold_error=.05,
                          max_rel_sigma_error=.30, min_sigma_adc=.6, max_sigma_adc=29.5),
    'quality_strict': dict(max_chi2ndf=3., max_rel_threshold_error=.02,
                          max_rel_sigma_error=.20, min_sigma_adc=.6, max_sigma_adc=29.5),
}
BOUNDS = ('max_chi2ndf', 'max_rel_threshold_error', 'max_rel_sigma_error',
          'min_threshold_adc', 'max_threshold_adc', 'min_sigma_adc', 'max_sigma_adc',
          'max_sigma_over_threshold')
FIELDS = ('threshold', 'thresholdErr', 'sigma', 'sigmaErr', 'chi2ndf')
IDS = np.array([l*100000+c*10000+h for l in range(30) for c in range(6)
                for h in range(30 if c==5 else 36)],dtype=np.int32)


def add_arguments(parser):
    parser.add_argument('--threshold-cut', default='allChannel',
                        help='Output name: preset or custom name with explicit numeric bounds')
    parser.add_argument('--threshold-root', type=Path,
                        default=Path(__file__).resolve().parents[1]/'result/threshold/threshold.root')
    parser.add_argument('--threshold-tree', default='expErfThre',
                        help='Current digitizer uses expErfThre; other models select a different quality mask')
    for key in BOUNDS:
        parser.add_argument('--threshold-'+key.replace('_','-'),type=float,default=None)


def settings(args):
    name=args.threshold_cut
    if not re.fullmatch('[A-Za-z0-9_-]+',name): raise RuntimeError('Invalid threshold-cut name')
    supplied={k:getattr(args,'threshold_'+k) for k in BOUNDS if getattr(args,'threshold_'+k) is not None}
    if name=='allChannel':
        if supplied: raise RuntimeError('allChannel cannot have quality bounds; use a different cut name')
        return None
    if name not in PRESETS and not supplied:
        raise RuntimeError('Unknown threshold-cut name; provide numeric bounds or use: '+','.join(PRESETS))
    result=dict(PRESETS.get(name) or {}); result.update(supplied)
    for key,value in result.items():
        if not np.isfinite(value) or value<0 or (key.startswith('max_') and value==0):
            raise RuntimeError('Invalid threshold quality bound: '+key)
    for quantity in ('threshold_adc','sigma_adc'):
        lo=result.get('min_'+quantity);hi=result.get('max_'+quantity)
        if lo is not None and hi is not None and lo>hi:raise RuntimeError('Minimum exceeds maximum: '+quantity)
    return result


def select(args):
    bounds=settings(args)
    catalog={k:np.full(len(IDS),np.nan) for k in FIELDS}
    present=np.zeros(len(IDS),dtype=bool); valid=np.zeros(len(IDS),dtype=bool)
    reasons={};source=None
    if bounds is None:
        selected=np.ones(len(IDS),dtype=bool)
    else:
        p=args.threshold_root.resolve(strict=True);before=p.stat()
        if not args.threshold_tree or ';' in args.threshold_tree:raise RuntimeError('Invalid threshold tree name')
        with uproot.open(p) as f:
            if args.threshold_tree not in f:raise RuntimeError('Missing threshold tree: '+args.threshold_tree)
            tree=f[args.threshold_tree]
            needed=['CellID',*FIELDS]
            if not all(k in tree.keys() for k in needed):raise RuntimeError('Threshold tree lacks quality branches')
            arrays=tree.arrays(needed,library='np')
        after=p.stat()
        if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise RuntimeError('Threshold input changed during read')
        raw=np.asarray(arrays['CellID'],dtype=np.int64)
        if len(np.unique(raw))!=len(raw):raise RuntimeError('Duplicate threshold CellID rows')
        index=np.searchsorted(IDS,raw);connected=(index<len(IDS))
        connected &= IDS[np.minimum(index,len(IDS)-1)]==raw
        for k in FIELDS:catalog[k][index[connected]]=arrays[k][connected]
        present[index[connected]]=True
        valid=present & np.logical_and.reduce([np.isfinite(catalog[k]) for k in FIELDS])
        valid &= (catalog['threshold']>0)&(catalog['sigma']>0)&(catalog['thresholdErr']>=0)&(catalog['sigmaErr']>=0)&(catalog['chi2ndf']>=0)
        selected=valid.copy()
        ratios={
            'max_rel_threshold_error':np.divide(catalog['thresholdErr'],catalog['threshold'],out=np.full(len(IDS),np.nan),where=catalog['threshold']>0),
            'max_rel_sigma_error':np.divide(catalog['sigmaErr'],catalog['sigma'],out=np.full(len(IDS),np.nan),where=catalog['sigma']>0),
            'max_sigma_over_threshold':np.divide(catalog['sigma'],catalog['threshold'],out=np.full(len(IDS),np.nan),where=catalog['threshold']>0),
        }
        checks={'max_chi2ndf':catalog['chi2ndf'],**ratios,
                'min_threshold_adc':catalog['threshold'],'max_threshold_adc':catalog['threshold'],
                'min_sigma_adc':catalog['sigma'],'max_sigma_adc':catalog['sigma']}
        for key,value in bounds.items():
            passing=checks[key]>=value if key.startswith('min_') else checks[key]<=value
            reasons[key]=int((valid & ~passing).sum());selected &= passing
        if not selected.any():raise RuntimeError('No channels pass this threshold quality cut')
        payload={k:np.asarray(arrays[k]).tolist() for k in needed}
        source=dict(path=str(p),bytes=before.st_size,mtime_ns=before.st_mtime_ns,tree=args.threshold_tree,
                    fit_rows=len(raw),unconnected_rows=int((~connected).sum()),
                    fit_table_sha256=hashlib.sha256(json.dumps(payload,sort_keys=True,allow_nan=True).encode()).hexdigest())
    counts=np.bincount(IDS[selected]//100000,minlength=30)
    selected_ids=IDS[selected].tolist()
    info=dict(name=args.threshold_cut,criteria=bounds,source=source,total_channels=len(IDS),
              selected_channels=int(selected.sum()),selected_by_layer=counts.tolist(),
              selected_ids_sha256=hashlib.sha256(np.asarray(selected_ids,dtype='<i4').tobytes()).hexdigest(),
              selected_ids=selected_ids,empty_layers=np.flatnonzero(counts==0).tolist(),
              fit_quality_evaluated=bounds is not None,
              scope='Same physical channel mask for data, MC truth, digitization and calibration',
              event_selection='Keep original axis-valid events and original entrance axes, including zero-hit events after masking')
    if bounds is not None:
        info.update(fit_present_channels=int(present.sum()),valid_fit_channels=int(valid.sum()),
                    missing_fit_channels=int((~present).sum()),invalid_fit_channels=int((present & ~valid).sum()),
                    failed_bounds_among_valid=reasons,
                    success_definition='Producer writes only successful fits. Require finite values, threshold/sigma > 0, errors/chi2ndf >= 0; unavailable fit status/covariance/NDF cannot be independently rechecked.')
    rows=[]
    for i,id in enumerate(IDS):
        failed=[]
        if bounds is not None:
            if not present[i]:failed.append('missing_fit')
            elif not valid[i]:failed.append('invalid_fit')
            else:
                for key,value in bounds.items():
                    ok=checks[key][i]>=value if key.startswith('min_') else checks[key][i]<=value
                    if not ok:failed.append(key)
        rows.append(dict(CellID=int(id),layer=int(id//100000),chip=int(id%100000//10000),channel=int(id%100),
                         selected=int(selected[i]),fit_present=int(present[i]) if bounds is not None else '',valid_fit=int(valid[i]) if bounds is not None else '',
                         **{k:float(catalog[k][i]) for k in FIELDS},rejected_by=','.join(failed)))
    return info,rows


def write_table(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
