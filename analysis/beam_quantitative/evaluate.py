"""Event bootstrap and practical-equivalence evaluation; NumPy/uproot only."""
import csv
import contextlib
import io
import json
import math
from pathlib import Path
import awkward as ak
import numpy as np
import uproot

EDGES = np.linspace(-116.6, 116.6, 45)  # 44 bins, exactly 5.3 mm; shared by data and MC
SSA_LAB_PROFILE_EDGES = (np.arange(43, dtype=float) - 21.) * 5.3
# Entrance subtraction makes positions continuous. Use the same 5.3 mm display
# width, with a bin centered at zero and enough range for the complete detector.
SSA_ALIGNED_PROFILE_EDGES = (np.arange(86, dtype=float) - 42.5) * 5.3
SSA_PROFILE_BINNING = dict(version=1, width_mm=5.3, lab_bins=42, aligned_bins=85,
                           lab_range_mm=[-111.3,111.3], aligned_range_mm=[-225.25,225.25],
                           rule='Lab bins follow strip pitch; aligned bins use the same width around the entrance centroid')
TOLERANCES = dict(mean=0.5, rms_relative=0.05, q90_relative=0.10,
                  tail_fraction=0.01, correlation=0.05, tv=0.05,
                  hits_relative=0.10)

@contextlib.contextmanager
def buffered_root(path):
    # Many small ROOT metadata writes are slow on the shared filesystem.
    # Histograms/metric tables are small: assemble them in RAM, then write once.
    class Buffer(io.BytesIO):
        def close(self):
            pass
    buffer = Buffer()
    try:
        with uproot.recreate(buffer) as root:
            yield root
        Path(path).write_bytes(buffer.getvalue())
    finally:
        io.BytesIO.close(buffer)

def load(path):
    names = ['axisValid', 'axisX', 'axisY', 'layer', 'x', 'y', 'energy',
             'missingParents', 'invalidPre', 'invalidPost']
    with uproot.open(path) as f:
        a = f['events'].arrays(names, library='ak')
    valid = np.asarray(a.axisValid, dtype=bool)
    audit = dict(sampled=len(valid), accepted=int(valid.sum()),
                 axis_invalid=int((~valid).sum()),
                 events_missing_parents=int(ak.sum(a.missingParents > 0)),
                 invalid_pre_hits=int(ak.sum(a.invalidPre)), invalid_post_hits=int(ak.sum(a.invalidPost)))
    a = a[valid]
    n = len(a)
    if not n:
        raise RuntimeError('No events with a valid entrance axis')
    sizes = np.asarray(ak.num(a.layer), dtype=np.int32)
    return dict(n=n, audit=audit, ax=np.asarray(a.axisX), ay=np.asarray(a.axisY),
                event=np.repeat(np.arange(n, dtype=np.int32), sizes),
                layer=np.asarray(ak.flatten(a.layer), dtype=np.int32),
                x=np.asarray(ak.flatten(a.x)), y=np.asarray(ak.flatten(a.y)),
                energy=np.maximum(0., np.asarray(ak.flatten(a.energy))))

def ssa_profile_histograms(sample):
    """Display histograms only; keep exact-position metrics and TV bins unchanged."""
    histograms={}
    for layer in range(30):
        axis='x' if layer % 2 else 'y'
        selected=sample['layer']==layer
        lab=sample[axis][selected]
        aligned=lab-sample['a'+axis][sample['event'][selected]]
        for weight in ('count','positive_energy'):
            values=np.ones(len(lab)) if weight=='count' else sample['energy'][selected]
            for frame,coordinates,edges in (('lab',lab,SSA_LAB_PROFILE_EDGES),('aligned',aligned,SSA_ALIGNED_PROFILE_EDGES)):
                if np.any((coordinates<edges[0]) | (coordinates>=edges[-1])):
                    raise RuntimeError('SSA profile positions outside detector display range')
                histograms[frame+'/'+weight+'/layer_'+str(layer)+'/u']=np.histogram(
                    coordinates,bins=edges,weights=values/sample['n'])
    return histograms


def refresh_ssa_profile_bins(histogram_path, events_path):
    """Refresh cached profiles from the same events without repeating evaluation."""
    with uproot.open(histogram_path) as source:
        if 'binning/ssa_profiles' in source and json.loads(str(source['binning/ssa_profiles']))==SSA_PROFILE_BINNING:
            return False
        retained={}
        for key,classname in source.classnames(cycle=False).items():
            if classname=='TDirectory':continue
            retained[key]=str(source[key]) if classname=='TObjString' else source[key]
    replacements=ssa_profile_histograms(load(events_path))
    with buffered_root(histogram_path) as output:
        for key,value in retained.items():
            if key not in replacements and key!='binning/ssa_profiles':output[key]=value
        for key,value in replacements.items():output[key]=value
        output['binning/ssa_profiles']=json.dumps(SSA_PROFILE_BINNING,sort_keys=True)
    return True


def resamples(n, b, rng):
    # Ordinary nonparametric bootstrap: retain the whole event, including all hits.
    return np.vstack([np.ones(n), rng.multinomial(n, np.full(n, 1./n), size=b)]).astype(float)

def matrix(n, event, index, k, weight):
    return np.bincount(event.astype(np.int64)*k + index, weights=weight,
                       minlength=n*k).reshape(n, k)

def map_index(x, y):
    ix = np.searchsorted(EDGES, x, side='right') - 1
    iy = np.searchsorted(EDGES, y, side='right') - 1
    inside = (ix >= 0) & (ix < 44) & (iy >= 0) & (iy < 44)
    return np.where(inside, ix*44+iy, 44*44)  # retain flow as one explicit category

def tv(a, b):
    sa, sb = a.sum(axis=1), b.sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        return .5*np.abs(a/sa[:, None]-b/sb[:, None]).sum(axis=1)

def bootstrap_tv(a, b, bd, bm):
    # Empty categories contribute zero; omit them before dense multiplication.
    active = (a.sum(axis=0) + b.sum(axis=0)) > 0
    return tv(bd @ a[:, active], bm @ b[:, active])

def quantile(v, w, q=.9):
    if not np.isfinite(w).all() or w.sum() <= 0:
        return np.nan
    order = np.argsort(v, kind='stable')
    j = np.searchsorted(np.cumsum(w[order]), q*w.sum(), side='left')
    return v[order[min(j, len(order)-1)]]

def statistics(v, p):
    total = p.sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        mean = (p @ v)/total
        rms = np.sqrt(np.maximum(0., (p @ (v*v))/total-mean*mean))
    q90 = np.array([quantile(np.abs(v-mu), row) for mu, row in zip(mean, p)])
    return dict(mean=mean, rms=rms, q90=q90)

def relative(a, b):
    with np.errstate(divide='ignore', invalid='ignore'):
        return b/a-1.

def correlation(x, y, w):
    n = w.sum(axis=1)
    mx, my = (w @ x)/n, (w @ y)/n
    vx = np.maximum(0., (w @ (x*x))/n-mx*mx)
    vy = np.maximum(0., (w @ (y*y))/n-my*my)
    with np.errstate(divide='ignore', invalid='ignore'):
        return ((w @ (x*y))/n-mx*my)/np.sqrt(vx*vy)

def json_safe(value):
    if isinstance(value, dict): return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)): return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer): return int(value)
    return value

def evaluate(data_path, mc_path, data_hist, mc_hist, metrics_root, table_path,
             bootstraps=500, seed=20261011, min_events=200, gate_layers=(0,1,2,3), tolerances=None):
    tol = dict(TOLERANCES)
    if tolerances:
        unknown = set(tolerances)-set(tol)
        if unknown: raise RuntimeError('Unknown tolerance keys: '+str(sorted(unknown)))
        tol.update(tolerances)
    if any(not np.isfinite(v) or v <= 0 for v in tol.values()):
        raise RuntimeError('Tolerances must be finite and positive')
    d, m = load(data_path), load(mc_path)
    rng = np.random.default_rng(seed)
    bd, bm = resamples(d['n'], bootstraps, rng), resamples(m['n'], bootstraps, rng)
    rows = []
    sufficient = min(d['n'], m['n']) >= min_events

    def add(group, layer, axis, frame, weight, metric, a, b, delta, gate):
        a, b, delta = np.asarray(a), np.asarray(b), np.asarray(delta)
        samples = delta[1:]
        finite = np.isfinite(samples)
        # Undefined metrics or degenerate resamples must not produce a false PASS.
        complete = bool(np.isfinite(delta[0]) and finite.all())
        lo, hi = np.quantile(samples, [.025,.975]) if complete else (np.nan,np.nan)
        upper = max(abs(lo),abs(hi)) if complete else np.nan
        limit = tol[metric]
        status = 'UNDETERMINED'
        if complete and sufficient:
            if upper <= limit: status = 'PASS'
            elif lo > limit or hi < -limit: status = 'FAIL'
        rows.append(dict(group=group, layer=layer, axis=axis, frame=frame, weight=weight,
                         metric=metric, gate=bool(gate),
                         value_unit=('mm' if metric in ('mean','rms_relative','q90_relative') else 'positions/event' if metric=='hits_relative' else 'dimensionless'),
                         difference_unit=('mm' if metric=='mean' else 'relative' if metric.endswith('_relative') else 'dimensionless'),
                         data_value=float(a[0]), mc_value=float(b[0]),
                         difference=float(delta[0]), ci_low=float(lo), ci_high=float(hi),
                         tolerance=limit, upper_abs=float(upper), score=float(upper/limit),
                         status=status, valid_bootstraps=int(finite.sum()),
                         data_events=d['n'], mc_events=m['n']))

    with buffered_root(data_hist) as dh, buffered_root(mc_hist) as mh:
        for store in (dh,mh):store['binning/ssa_profiles']=json.dumps(SSA_PROFILE_BINNING,sort_keys=True)
        for axis in ('x','y'):
            dv, mv = d['a'+axis], m['a'+axis]
            ds, ms = statistics(dv, bd), statistics(mv, bm)
            for stat, metric in [('mean','mean'),('rms','rms_relative'),('q90','q90_relative')]:
                delta = ms[stat]-ds[stat] if stat == 'mean' else relative(ds[stat],ms[stat])
                add('entrance',-1,axis,'lab','event',metric,ds[stat],ms[stat],delta,True)
            # Both tails use the DATA center, including its bootstrap uncertainty.
            dt = (bd*(np.abs(dv[None,:]-ds['mean'][:,None]) >= 30)).sum(axis=1)/d['n']
            mt = (bm*(np.abs(mv[None,:]-ds['mean'][:,None]) >= 30)).sum(axis=1)/m['n']
            add('entrance',-1,axis,'lab','event','tail_fraction',dt,mt,mt-dt,True)
            dh['entrance/'+axis] = np.histogram(dv,bins=EDGES)
            mh['entrance/'+axis] = np.histogram(mv,bins=EDGES)
        dc, mc = correlation(d['ax'],d['ay'],bd), correlation(m['ax'],m['ay'],bm)
        add('entrance',-1,'xy','lab','event','correlation',dc,mc,mc-dc,True)
        dx = matrix(d['n'],np.arange(d['n']),map_index(d['ax'],d['ay']),1937,np.ones(d['n']))
        mx = matrix(m['n'],np.arange(m['n']),map_index(m['ax'],m['ay']),1937,np.ones(m['n']))
        dist = bootstrap_tv(dx,mx,bd,bm)
        add('entrance',-1,'xy','lab','event','tv',np.zeros_like(dist),dist,dist,True)
        dh['entrance/xy'] = (dx.sum(axis=0)[:-1].reshape(44,44),EDGES,EDGES)
        mh['entrance/xy'] = (mx.sum(axis=0)[:-1].reshape(44,44),EDGES,EDGES)
        dh['entrance/flow'] = (np.array([dx[:,-1].sum()]),np.array([0.,1.]))
        mh['entrance/flow'] = (np.array([mx[:,-1].sum()]),np.array([0.,1.]))
        del dx,mx
        for layer in range(30):
            gate = layer in gate_layers
            axis = 'x' if layer % 2 else 'y'
            di, mi = d['layer'] == layer, m['layer'] == layer
            de, me = d['event'][di], m['event'][mi]
            dv, mv = d[axis][di], m[axis][mi]
            # Fine-axis coordinates are strip centers. Merge insignificant roundoff
            # at 1e-6 mm; quantiles retain the actual discrete geometry, not plot bins.
            dr, mr = np.round(dv,6), np.round(mv,6)
            values = np.unique(np.concatenate([dr,mr]))
            if len(values) > 10000: raise RuntimeError('Unexpectedly many fine-axis positions')
            idxd, idxm = np.searchsorted(values,dr), np.searchsorted(values,mr)
            for weight in ('count','positive_energy'):
                dw = np.ones(len(dv)) if weight == 'count' else d['energy'][di]
                mw = np.ones(len(mv)) if weight == 'count' else m['energy'][mi]
                dp = matrix(d['n'],de,idxd,len(values),dw)
                mp = matrix(m['n'],me,idxm,len(values),mw)
                dsample, msample = bd @ dp, bm @ mp
                ds, ms = statistics(values,dsample), statistics(values,msample)
                for stat, metric in [('mean','mean'),('rms','rms_relative'),('q90','q90_relative')]:
                    delta = ms[stat]-ds[stat] if stat == 'mean' else relative(ds[stat],ms[stat])
                    add('hits',layer,axis,'lab',weight,metric,ds[stat],ms[stat],delta,gate)
                with np.errstate(divide='ignore',invalid='ignore'):
                    dt = (dsample*(np.abs(values[None,:]-ds['mean'][:,None]) >= 30)).sum(axis=1)/dsample.sum(axis=1)
                    mt = (msample*(np.abs(values[None,:]-ds['mean'][:,None]) >= 30)).sum(axis=1)/msample.sum(axis=1)
                add('hits',layer,axis,'lab',weight,'tail_fraction',dt,mt,mt-dt,gate)
                if weight == 'count':
                    dn, mn = dsample.sum(axis=1)/d['n'], msample.sum(axis=1)/m['n']
                    add('hits',layer,'none','lab',weight,'hits_relative',dn,mn,relative(dn,mn),gate)
                # Alignment is a diagnostic, never an input-beam tuning gate.
                def aligned(p, sample, boot, coord):
                    w, s, s2 = p.sum(axis=1),p @ values,p @ (values*values)
                    aa = sample['a'+coord]
                    total = boot @ w
                    with np.errstate(divide='ignore',invalid='ignore'):
                        mean = (boot @ (s-aa*w))/total
                        var = (boot @ (s2-2*aa*s+aa*aa*w))/total-mean*mean
                    return mean,np.sqrt(np.maximum(0.,var))
                da, dar = aligned(dp,d,bd,axis); ma, mar = aligned(mp,m,bm,axis)
                add('hits',layer,axis,'aligned',weight,'mean',da,ma,ma-da,False)
                add('hits',layer,axis,'aligned',weight,'rms_relative',dar,mar,relative(dar,mar),False)
                del dp,mp,dsample,msample
                dx = matrix(d['n'],de,map_index(d['x'][di],d['y'][di]),1937,dw)
                mx = matrix(m['n'],me,map_index(m['x'][mi],m['y'][mi]),1937,mw)
                dist = bootstrap_tv(dx,mx,bd,bm)
                add('hits',layer,'xy','lab',weight,'tv',np.zeros_like(dist),dist,dist,gate)
                for store, sample, sel, w, mat in [(dh,d,di,dw,dx),(mh,m,mi,mw,mx)]:
                    prefix = 'lab/'+weight+'/layer_'+str(layer)
                    store[prefix+'/xy'] = (mat.sum(axis=0)[:-1].reshape(44,44)/sample['n'],EDGES,EDGES)
                    store[prefix+'/flow_per_event'] = (np.array([mat[:,-1].sum()/sample['n']]),np.array([0.,1.]))
                    store[prefix+'/u'] = np.histogram(sample[axis][sel],bins=SSA_LAB_PROFILE_EDGES,weights=w/sample['n'])
                    store['aligned/'+weight+'/layer_'+str(layer)+'/u'] = np.histogram(sample[axis][sel]-sample['a'+axis][sample['event'][sel]],bins=SSA_ALIGNED_PROFILE_EDGES,weights=w/sample['n'])
                del dx,mx
            print(f'Quantitative evaluation: layer {layer}, events data={d["n"]}, MC={m["n"]}',flush=True)
    with table_path.open('w',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    with buffered_root(metrics_root) as f:
        schema = {k: ('string' if isinstance(rows[0][k],str) else 'bool' if isinstance(rows[0][k],bool) else 'int64' if isinstance(rows[0][k],int) else 'float64') for k in rows[0]}
        tree = f.mktree('metrics',schema)
        tree.extend({k: ([r[k] for r in rows] if schema[k]=='string' else np.array([r[k] for r in rows])) for k in schema})
    gates = [r for r in rows if r['gate']]
    def decision(rr):
        if any(r['status']=='FAIL' for r in rr): return 'FAIL'
        if not rr or any(r['status']!='PASS' for r in rr): return 'UNDETERMINED'
        return 'PASS'
    score = max(r['score'] for r in gates) if all(np.isfinite(r['score']) for r in gates) else np.nan
    worst = sorted(gates,key=lambda r:r['score'] if np.isfinite(r['score']) else float('inf'),reverse=True)[:10]
    return json_safe(dict(beam_status=decision(gates),all_layer_lab_status=decision([r for r in rows if r['frame']=='lab']),
                          aligned_diagnostic_status=decision([r for r in rows if r['frame']=='aligned']),
                          score=score,tolerances=tol,gate_layers=list(gate_layers),bootstrap_replicates=bootstraps,
                          confidence_level=.95,confidence_method='two-sided percentile event bootstrap; approximate',
                          sufficient_events=sufficient,min_events=min_events,data=d['audit'],mc=m['audit'],
                          worst_metrics=worst,metric_count=len(rows),
                          interpretation='A PASS is conditional on the fixed MC/material/selection model; not proof of true beam distribution.'))
