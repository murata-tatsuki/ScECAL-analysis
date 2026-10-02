"""MIP policies, common phase-space weights, and a shared Landau-Gaussian fit."""
import csv
import math
from pathlib import Path

import numpy as np


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    with Path(path).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def physical(c):
    c = int(c)
    if c < 0 or c // 100000 >= 32 or c % 100000 // 10000 >= 6 or c % 100 >= 36:
        raise ValueError(f'Invalid CellID {c}')
    return c // 100000 * 100000 + c % 100000 // 10000 * 10000 + c % 100


def pitch(c):
    return 10 if 4 <= c // 100000 <= 27 else 15


def ieee_div(a, b):
    with np.errstate(divide='ignore', invalid='ignore'):
        return float(np.divide(float(a), float(b)))


def policy_audit(source_rows):
    """Reproduce legacy C++ including NaN comparisons and the value==1 sentinel.

    The clean policy keeps existing chi2 cuts and additionally rejects nonpositive
    MPV/NDF, negative chi2 and nonfinite fit fields. Physical connected channels
    are the clean fallback pool; legacy pools reproduce all source rows < layer30.
    """
    source = {}
    for row in source_rows:
        c = physical(row['cellid'])
        if c // 100000 >= 30:
            continue
        if c in source:
            raise ValueError(f'Duplicate MIP constant for {c}')
        r = dict(row, cellid=c, pitch_um=pitch(c))
        r['chi2_ndf'] = ieee_div(r['chi2'], r['ndf'])
        cut = 2 if pitch(c) == 10 else 1.7
        r['legacy_accept'] = not r['chi2_ndf'] > cut
        r['digi_accept'] = not r['mpv'] < 0  # zero and NaN pass the original code
        reasons = []
        if not all(math.isfinite(r[k]) for k in ('mpv', 'chi2', 'ndf', 'width', 'sigma')):
            reasons.append('nonfinite')
        if r['mpv'] <= 0:
            reasons.append('nonpositive_mpv')
        if r['ndf'] <= 0:
            reasons.append('nonpositive_ndf')
        if r['chi2'] < 0:
            reasons.append('negative_chi2')
        if math.isfinite(r['chi2_ndf']) and r['chi2_ndf'] > cut:
            reasons.append('chi2_cut')
        if c % 100000 // 10000 == 5 and c % 100 >= 30:
            reasons.append('unconnected')
        r['clean_accept'] = not reasons
        r['clean_reason'] = ','.join(reasons) or 'accepted'
        source[c] = r
    means = {}
    for p in (10, 15):
        for policy in ('legacy', 'digi', 'clean'):
            values = [r['mpv'] for r in source.values() if r['pitch_um'] == p and r[policy + '_accept']]
            if not values:
                raise ValueError(f'No fallback donors for {policy}, {p} um')
            means[p, policy] = {'pitch_um': p, 'policy': policy, 'n_donors': len(values),
                               'nonpositive_donors': sum(v <= 0 for v in values),
                               'nonfinite_donors': sum(not math.isfinite(v) for v in values),
                               'mean_adc': sum(values) / len(values)}
    rows = []
    for layer in range(30):
        for chip in range(6):
            for channel in range(30 if chip == 5 else 36):
                c = layer * 100000 + chip * 10000 + channel
                r = source.get(c)
                out = dict(cellid=c, pitch_um=pitch(c), input_mpv_adc=r['mpv'] if r else math.nan,
                           ndf=r['ndf'] if r else 0, chi2_ndf=r['chi2_ndf'] if r else math.nan,
                           clean_reason=r['clean_reason'] if r else 'missing')
                for policy in ('legacy', 'digi', 'clean'):
                    accepted = bool(r and r[policy + '_accept'])
                    sentinel = policy != 'clean' and accepted and r['mpv'] == 1
                    fallback = not accepted or sentinel
                    value = means[pitch(c), policy]['mean_adc'] if fallback else r['mpv']
                    out[policy + '_source'] = 'fallback' if fallback else 'individual'
                    out[policy + '_adc'] = value
                out['digi_runtime_override'] = out['digi_adc'] <= 5
                if out['digi_runtime_override']:
                    out['digi_adc'] = 90. if pitch(c) == 10 else 400.
                    out['digi_source'] = 'runtime_default'
                out['changed'] = not (out['legacy_adc'] == out['clean_adc'])
                rows.append(out)
    return rows, list(means.values())


def read_mip(ROOT, path):
    f = ROOT.TFile.Open(str(path)); t = f.Get('MIP_Fit') if f else None
    if not t:
        raise ValueError(f'Missing MIP_Fit in {path}')
    sigma_branch = 'GauSigma' if t.GetBranch('GauSigma') else 'GausSigma'
    rows = [dict(cellid=int(e.CellID), mpv=float(e.LandauMPV), chi2=float(e.ChiSquare),
                 ndf=int(e.NDF), width=float(e.LandauWidth), sigma=float(getattr(e, sigma_branch))) for e in t]
    f.Close()
    return policy_audit(rows)


def read_pedestal(ROOT, path):
    f = ROOT.TFile.Open(str(path)); t = f.Get('ChnLevel') if f else None
    if not t:
        raise ValueError(f'Missing ChnLevel in {path}')
    result = {}
    for e in t:
        if len(e.CellID) != len(e.PedHighMean):
            raise ValueError('Pedestal vector size mismatch')
        for c, v in zip(e.CellID, e.PedHighMean):
            c, v = physical(c), float(np.float32(v))  # production stores float arrays
            if not math.isfinite(v) or (c in result and result[c] != v):
                raise ValueError(f'Invalid/conflicting pedestal for {c}')
            result[c] = v
    f.Close()
    return result


def read_threshold(ROOT, path):
    f = ROOT.TFile.Open(str(path)); t = f.Get('expErfThre') if f else None
    if not t:
        raise ValueError(f'Missing expErfThre in {path}')
    result, pools = {}, {}
    for e in t:
        c, threshold, sigma = physical(e.CellID), float(e.threshold), float(e.sigma)
        if c // 100000 >= 30 or not math.isfinite(threshold + sigma) or sigma <= 0:
            continue
        if c in result:
            raise ValueError(f'Duplicate threshold for {c}')
        result[c] = (threshold, sigma, 'individual')
        pools.setdefault(c // 10000, []).append((threshold, sigma))
    f.Close()
    for c in (l*100000+ch*10000+i for l in range(30) for ch in range(6) for i in range(30 if ch == 5 else 36)):
        if c not in result and c // 10000 in pools:
            a = np.mean(pools[c // 10000], axis=0)
            result[c] = (float(a[0]), float(a[1]), 'chip_mean')
    return result


def common_weights(a, long_bin=5., transverse_bin=1.5, slope_bin=.01, min_bin=3):
    """Symmetric overlap target min(Ndata,Nmc), at fixed run/channel/position/slopes.

    Weights depend on geometry and counts, never on ADC, energy, or either policy.
    A run/geometry cell missing in either sample contributes zero to both fits.
    """
    keys = np.column_stack([a['run'], a['cellid'],
                            np.floor(a['longitudinal']/long_bin),
                            np.floor(a['transverse']/transverse_bin),
                            np.floor(a['sx']/slope_bin + .5),
                            np.floor(a['sy']/slope_bin + .5)]).astype(np.int64)
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    counts = [np.bincount(inverse[a['sample'] == s], minlength=len(unique)) for s in (0, 1)]
    target = np.minimum(*counts)
    target[target < min_bin] = 0
    weight = np.zeros(len(inverse), dtype=float)
    for s in (0, 1):
        mask = a['sample'] == s
        weight[mask] = target[inverse[mask]] / counts[s][inverse[mask]]
    rows = [dict(run=int(k[0]), cellid=int(k[1]), long_bin=int(k[2]), transverse_bin=int(k[3]),
                 sx_bin=int(k[4]), sy_bin=int(k[5]), n_data=int(counts[0][i]), n_mc=int(counts[1][i]),
                 target=int(target[i])) for i, k in enumerate(unique)]
    return weight, rows


def load_fit_function(ROOT):
    ROOT.gInterpreter.Declare(r'''
    #include <TMath.h>
    double muonLangau(double*x,double*p) {
        const int n=100; const double step=10.*p[3]/n;
        if(p[0]<=0||p[3]<=0)return 0;
        const double location=p[1]+0.22278298*p[0];
        double sum=0;
        for(int i=0;i<n;++i) {
            double y=x[0]-5*p[3]+(i+.5)*step;
            sum+=TMath::Landau(y,location,p[0],true)*TMath::Gaus(x[0],y,p[3],true);
        }
        return p[2]*step*sum;
    }
    ''')


def fit_spectrum(ROOT, rootout, name, values, weights, scale, low, high, min_effective=100):
    """Fit in scaled coordinates; outputs use the original units.

    Plateau-only fit avoids inventing efficiency below the hardware threshold.
    Reject boundary/covariance failures and an MPV outside the fitted interval.
    """
    finite = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[finite], weights[finite]
    neff = weights.sum()**2 / np.square(weights).sum() if len(weights) else 0.
    result = dict(n=len(values), sumw=float(weights.sum()), effective_n=float(neff),
                  fit_low=float(low), fit_high=float(high), mpv=math.nan, mpv_err=math.nan,
                  convolution_mode=math.nan, width=math.nan, sigma=math.nan,
                  chi2=math.nan, ndf=0, status=-1, covariance=-1, valid=False, reason='low_statistics')
    # 0.04 nominal MIP per bin; identical bin edges for data and MC and every stage.
    rootout.cd()
    hist = ROOT.TH1D(name, name, 225, -1., 8.)
    hist.Sumw2()
    if len(values):
        hist.FillN(len(values), np.ascontiguousarray(values/scale), np.ascontiguousarray(weights))
    result['underflow'] = float(hist.GetBinContent(0))
    result['overflow'] = float(hist.GetBinContent(hist.GetNbinsX()+1))
    in_fit = (values >= low) & (values <= high)
    wfit = weights[in_fit]
    fit_neff = wfit.sum()**2 / np.square(wfit).sum() if len(wfit) else 0.
    result['fit_effective_n'] = float(fit_neff)
    if fit_neff < min_effective or not 0 <= low < high <= 8*scale:
        hist.Write(); hist.SetDirectory(0)
        return result
    lo, hi = low/scale, high/scale
    # Align the interval with bin centers so that threshold-edge bins are not included.
    first = max(1, hist.FindBin(lo))
    while first <= hist.GetNbinsX() and hist.GetBinLowEdge(first) < lo:
        first += 1
    lo = hist.GetBinLowEdge(first)
    trials = []
    for seed in (1., max(lo+.08, min(hi-.08, float(np.median(values[in_fit])/scale)*.75))):
        func = ROOT.TF1(name+'_fit', ROOT.muonLangau, lo, hi, 4)
        func.SetParNames('LandauWidth', 'LandauMPV', 'Area', 'GaussianSigma')
        func.SetParameters(.12, seed, hist.Integral()*.04, .15)
        func.SetParLimits(0, .005, 1.5)
        func.SetParLimits(1, .01, hi)
        func.SetParLimits(2, 1.e-8, max(1., hist.Integral())*20.)
        func.SetParLimits(3, .005, 1.5)
        fit = hist.Fit(func, 'QRS0N')
        trials.append((int(fit), int(fit.CovMatrixStatus()), float(func.GetChisquare()), func))
    status, covariance, chi2, func = min(trials, key=lambda t: (t[0] != 0 or t[1] < 2, t[2]))
    mpv = float(func.GetParameter(1))*scale
    result.update(mpv=mpv, mpv_err=float(func.GetParError(1))*scale,
                  convolution_mode=float(func.GetMaximumX(lo, hi))*scale,
                  width=float(func.GetParameter(0))*scale, sigma=float(func.GetParameter(3))*scale,
                  chi2=chi2, ndf=int(func.GetNDF()), status=status, covariance=covariance)
    reasons = []
    if status != 0 or covariance < 2:
        reasons.append('fit_failure')
    if result['ndf'] <= 0:
        reasons.append('nonpositive_ndf')
    if not all(math.isfinite(result[k]) for k in ('mpv', 'mpv_err', 'width', 'sigma', 'chi2')) or mpv <= 0:
        reasons.append('invalid_parameter')
    if not low < mpv < high:
        reasons.append('mpv_outside_plateau_fit')
    if any(func.GetParameter(p) < .0051 or func.GetParameter(p) > 1.499 for p in (0, 3)):
        reasons.append('parameter_boundary')
    if result['ndf'] > 0 and chi2/result['ndf'] > 3:
        reasons.append('chi2_ndf_gt_3')
    result.update(valid=not reasons, reason=','.join(reasons) or 'ok')
    hist.GetListOfFunctions().Add(func.Clone(name+'_model'))
    hist.Write(); hist.SetDirectory(0)
    return result
