"""Numerical diagnostics. Slopes never mix physical channels or ADC bins."""
import math
import numpy as np


def factor(t, cell, reference=20., alpha10=1.6/135, alpha15=3.5/230):
    layer = np.asarray(cell)//100000
    alpha = np.where((layer >= 4) & (layer <= 27), alpha10, alpha15)
    return 1-alpha*(np.asarray(t)-reference)


def energies(h, l, t, cell, ph, pl, mip, gain, offset, **kwargs):
    f = factor(t, cell, **kwargs)
    return dict(raw=(h-ph)*.305/mip,
                hg=(h-ph)*f*.305/mip,
                lg=((l-pl)*f*gain+offset)*.305/mip,
                lg_offset_before=((l-pl)*gain+offset)*f*.305/mip,
                lg_raw=((l-pl)*gain+offset)*.305/mip)


def moments(values):
    v = np.asarray(values)
    v = v[np.isfinite(v)]
    return dict(n=len(v), mean=float(v.mean()) if len(v) else math.nan,
                sem=float(v.std(ddof=1)/np.sqrt(len(v))) if len(v)>1 else math.nan,
                median=float(np.median(v)) if len(v) else math.nan,
                q16=float(np.quantile(v, .16)) if len(v) else math.nan,
                q84=float(np.quantile(v, .84)) if len(v) else math.nan)


def slope(rows, value, error, reference=20., min_span=.5, fixed_run=False):
    """WLS with independent run intercepts optionally; run and T confounding is explicit.

    Errors are formal conditional errors, inflated if reduced chi2 exceeds one.
    A constant-temperature run cannot identify a within-run slope.
    """
    rows = [r for r in rows if all(math.isfinite(float(r[k])) for k in ('temperature', value, error))
            and float(r[error]) > 0]
    n = len(rows)
    result = dict(n_points=n, n_runs=len({r['run'] for r in rows}), temperature_span=0.,
                  slope=math.nan, slope_err=math.nan, slope_significance=math.nan,
                  value_at_reference=math.nan, chi2=math.nan, ndf=0, status='insufficient_points')
    if n < 3:
        return result
    t, y, err = (np.array([float(r[k]) for r in rows]) for k in ('temperature', value, error))
    span = float(np.ptp(t))
    result['temperature_span'] = span
    if span < min_span:
        result['status'] = 'insufficient_temperature_span'
        return result
    runs = np.array([r['run'] for r in rows])
    groups = np.unique(runs)
    if fixed_run:
        within_span = max(np.ptp(t[runs == g]) for g in groups)
        if within_span < min_span:
            result['status'] = 'run_temperature_confounded'
            return result
        x = np.column_stack([t-reference, *[(runs == g).astype(float) for g in groups]])
    else:
        x = np.column_stack([t-reference, np.ones(n)])
    xw, yw = x/err[:, None], y/err
    if n <= x.shape[1] or np.linalg.matrix_rank(xw) < x.shape[1]:
        result['status'] = 'rank_deficient_or_no_ndf'
        return result
    beta = np.linalg.lstsq(xw, yw, rcond=None)[0]
    chi2 = float(np.square((y-x@beta)/err).sum())
    ndf = n-x.shape[1]
    cov = np.linalg.inv(xw.T@xw)*max(1., chi2/ndf)
    se = float(np.sqrt(cov[0, 0]))
    result.update(slope=float(beta[0]), slope_err=se,
                  slope_significance=float(beta[0]/se) if se else math.nan,
                  value_at_reference=float(beta[1]) if not fixed_run else math.nan,
                  chi2=chi2, ndf=ndf, status='ok')
    return result


def grouped(rows, keys):
    result = {}
    for r in rows:
        result.setdefault(tuple(r[k] for k in keys), []).append(r)
    return result


def trends(rows, keys, value, error, reference, span):
    output = []
    for key, group in grouped(rows, keys).items():
        base = dict(zip(keys, key))
        scopes = [('pooled_runs', 'all', group, False), ('run_adjusted', 'all', group, True)]
        scopes += [('within_run', run, g, False) for (run,), g in grouped(group, ['run']).items()]
        for scope, run, points, fixed in scopes:
            output.append(dict(base, scope=scope, run=run, metric=value,
                               **slope(points, value, error, reference, span, fixed)))
    return output
