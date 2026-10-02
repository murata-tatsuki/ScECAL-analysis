#!/usr/bin/env python3
"""PS HG-only MIP temperature closure, then HG/LG and pedestal diagnostics."""
import argparse
import csv
import fcntl
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'mip_closure'))
from mip_tools import load_fit_function, fit_spectrum, read_threshold
from temperature_tools import factor, energies, moments, grouped, trends


def table(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows, fields=None):
    with Path(path).open('w') as f:
        writer = csv.DictWriter(f, fieldnames=fields or (list(rows[0]) if rows else ['status']), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True, help='TSV header: kind,raw,calib,track; kind=mip/response/pedestal')
    for key in ('pedestal', 'mip', 'hl'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--threshold', type=Path, help='expErfThre; required for MIP fits')
    p.add_argument('--output', type=Path, required=True, help='new output directory')
    p.add_argument('--channels', default='all', help='physical CellID CSV')
    p.add_argument('--events-per-file', type=int, default=0, help='systematic midpoint sample of raw entries; 0=all')
    p.add_argument('--temp-width', type=float, default=.5)
    p.add_argument('--min-temp-span', type=float, default=.5)
    p.add_argument('--min-fit-hits', type=int, default=100)
    p.add_argument('--min-hits', type=int, default=30)
    p.add_argument('--max-angle', type=float, default=.05, help='projected track angle in rad')
    p.add_argument('--reference-temperature', type=float, default=20.)
    p.add_argument('--coefficient-10', type=float, default=1.6/135)
    p.add_argument('--coefficient-15', type=float, default=3.5/230)
    p.add_argument('--mip-measurement-temperature', type=float, help='documented acquisition temperature; never inferred from current samples')
    p.add_argument('--mip-constant-temperature', type=float, help='temperature represented by constants, after any correction during MIP production')
    p.add_argument('--mip-temperature-evidence', type=Path, help='production log/config documenting the two temperatures')
    p.add_argument('--plateau-sigmas', type=float, default=3.)
    p.add_argument('--fit-high-mip', type=float, default=3.)
    p.add_argument('--hg-max', type=float, default=2200., help='MIP low-signal raw HG upper cut, <=2600')
    p.add_argument('--overlap-min', type=float, default=800., help='raw HG ADC')
    p.add_argument('--overlap-max', type=float, default=2200.)
    p.add_argument('--adc-width', type=float, default=100.)
    p.add_argument('--compare-models', action='store_true', help='also fit sign-flipped and reciprocal factors as diagnostic alternatives')
    p.add_argument('--no-plots', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    numeric = [v for v in vars(a).values() if isinstance(v, float)]
    if not all(math.isfinite(v) for v in numeric):
        p.error('All numeric options must be finite')
    if min(a.temp_width, a.min_temp_span, a.max_angle, a.adc_width) <= 0 or a.events_per_file < 0 or a.min_fit_hits < 20 or a.min_hits < 3:
        p.error('Invalid widths, statistics or sample limit')
    if not 0 < a.hg_max <= 2600 or not 0 <= a.overlap_min < a.overlap_max <= 2600 or not 0 < a.fit_high_mip <= 8 or a.plateau_sigmas < 0:
        p.error('Invalid HG/fit/overlap range')
    if a.channels != 'all':
        from mip_tools import physical
        for c in map(int, a.channels.split(',')):
            if physical(c) != c or c//100000 >= 30 or (c%100000//10000 == 5 and c%100 >= 30):
                p.error('channels must be connected physical CellIDs in layers 0-29')
    return a


def inputs(a):
    rows = table(a.manifest)
    if not rows or not {'kind', 'raw', 'calib', 'track'} <= rows[0].keys():
        raise ValueError('Manifest requires kind,raw,calib,track header and at least one row')
    seen = set()
    for r in rows:
        if r['kind'] not in ('mip', 'response', 'pedestal'):
            raise ValueError('Unknown input kind '+r['kind'])
        for key in ('raw', 'calib', 'track'):
            needed = key == 'raw' or (key == 'calib' and r['kind'] != 'pedestal') or (key == 'track' and r['kind'] == 'mip')
            if not needed:
                r[key] = '-'
                continue
            path = Path(r[key]).expanduser()
            if not path.is_absolute():
                path = a.manifest.resolve().parent/path
            if not path.is_file():
                raise ValueError(f'Missing {key}: {path}')
            r[key] = str(path.resolve())
        # Pedestal samples must not be counted twice through several roles.
        if r['raw'] in seen:
            raise ValueError('Duplicate raw file: '+r['raw']+'; each mip/response input also audits its pedestal')
        seen.add(r['raw'])
    if any(r['kind'] == 'mip' for r in rows) and not a.threshold:
        raise ValueError('--threshold is required for MIP plateau fits')
    for key in ('pedestal', 'mip', 'hl', 'threshold', 'mip_temperature_evidence'):
        path = getattr(a, key)
        if path and not path.is_file():
            raise ValueError(f'Missing {key}: {path}')
    return rows


def snapshot(a, rows):
    dest = a.output/'inputs'
    dest.mkdir()
    (dest/'configuration.json').write_text(json.dumps(vars(a), default=str, indent=2)+'\n')
    write(dest/'manifest.tsv', rows)
    hashes = []
    paths = [*HERE.glob('*.py'), *HERE.glob('*.cc'), *HERE.glob('*.hh'), HERE/'Makefile', HERE/'run.sh',
             HERE.parent/'adc_energy/FastPng.hh',
             HERE.parent/'mip_closure/mip_tools.py',
             HERE.parent/'calibration_residual/CalibrationAudit.hh', HERE.parent/'calibration_residual/RootInput.hh']
    paths += [getattr(a, k) for k in ('pedestal', 'mip', 'hl', 'threshold', 'mip_temperature_evidence') if getattr(a, k)]
    production = Path('/megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/src')
    # Current source snapshots are evidence of the inspected implementation,
    # not a claim that these files produced the input ROOT files.
    paths += [production/name for name in ('Calibration.cxx','MIPFit.cxx','EBUdecode.cxx') if (production/name).is_file()]
    for i, path in enumerate(paths):
        hashes.append(dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        if path.suffix != '.root':
            shutil.copyfile(path, dest/f'{i}_{path.name}')
    write(dest/'checksums.tsv', hashes)
    metadata = []
    for i, r in enumerate(rows):
        for key in ('raw', 'calib', 'track'):
            if r[key] != '-':
                s = Path(r[key]).stat()
                metadata.append(dict(file_index=i, role=key, path=r[key], size=s.st_size, mtime_ns=s.st_mtime_ns))
    write(dest/'data_files.tsv', metadata)
    ref, measured, represented = a.reference_temperature, a.mip_measurement_temperature, a.mip_constant_temperature
    status = 'unconfirmed' if represented is None else ('declared_match' if abs(represented-ref)<.01 else 'declared_mismatch')
    write(a.output/'reference_audit.tsv', [dict(reference_C=ref, mip_measurement_C=measured,
          mip_constant_C=represented, status=status, evidence=str(a.mip_temperature_evidence or 'UNCONFIRMED'),
          acquisition_matches_reference='unknown' if measured is None else abs(measured-ref)<.01,
          factor_29C_10um=float(factor(29,400000,ref,a.coefficient_10,a.coefficient_15)),
          factor_29C_15um=float(factor(29,0,ref,a.coefficient_10,a.coefficient_15)),
          note='Declarations/evidence are not automatic verification of saved-file production provenance')])


def arrays(ROOT, path):
    return ROOT.RDataFrame('Hits', str(path)).AsNumpy()


def indices(h, tw):
    if not len(h['cellid']):
        return []
    valid = (h['temperature_valid'] == 1) & np.isfinite(h['temperature'])
    tb = np.full(len(valid), -2147483648, dtype=np.int64)
    tb[valid] = np.floor(h['temperature'][valid]/tw).astype(np.int64)
    keys = np.column_stack([h['run'], h['cellid'], tb])
    unique, inv = np.unique(keys, axis=0, return_inverse=True)
    order = np.argsort(inv, kind='stable')
    bounds = np.searchsorted(inv[order], np.arange(len(unique)+1))
    return [(tuple(map(int, key)), order[bounds[i]:bounds[i+1]]) for i, key in enumerate(unique)]


def process(ROOT, a, fi, role, h, params, threshold, rootout):
    mips, residuals, closures, excluded = [], [], [], []
    kw = dict(reference=a.reference_temperature, alpha10=a.coefficient_10, alpha15=a.coefficient_15)
    for (run, cell, tb), idx in indices(h, a.temp_width):
        p = params[cell]
        t, rh, rl = (h[k][idx] for k in ('temperature', 'raw_hg', 'raw_lg'))
        base = dict(file_index=fi, kind=role, run=run, cellid=cell, temp_bin=tb,
                    temperature=float(np.mean(t)) if tb != -2147483648 else math.nan)
        reasons = []
        if not p.pedPresent:
            reasons.append('missing_pedestal')
        if not math.isfinite(p.mip):
            reasons.append('nonfinite_mip')
        elif p.mip <= 0:
            reasons.append('nonpositive_mip')
        if reasons:
            # Audit legacy constants without replacing them or aborting other channels.
            # No division/fit with an invalid MIP, even when it passed legacy chi2 cuts.
            excluded.append(dict(base, n=len(idx), reason=','.join(reasons),
                                 mip_constant=float(p.mip), ped_present=bool(p.pedPresent)))
            continue
        en = energies(rh, rl, t, cell, p.ph, p.pl, p.mip, p.gain, p.offset, **kw)
        f = factor(t, cell, **kw)
        dh, dl = en['hg']-h['saved_hg'][idx], en['lg']-h['saved_lg'][idx]
        # Closure always uses the exact current production formula (20 C), even during a hypothesis scan.
        prod = energies(rh, rl, t, cell, p.ph, p.pl, p.mip, p.gain, p.offset)
        cdh, cdl = prod['hg']-h['saved_hg'][idx], prod['lg']-h['saved_lg'][idx]
        def absmax(v):
            v = v[np.isfinite(v)]
            return float(np.max(np.abs(v))) if len(v) else math.nan
        closures.append(dict(base, n=len(idx), temperature_valid=tb != -2147483648,
            production_hg_absmax=absmax(cdh), production_lg_absmax=absmax(cdl),
            production_hg_mean=moments(cdh)['mean'], production_lg_mean=moments(cdl)['mean'],
            production_hg_mismatch=int(np.count_nonzero(~np.isfinite(cdh) | (np.abs(cdh)>1e-7+1e-6*np.abs(h['saved_hg'][idx])))),
            production_lg_mismatch=int(np.count_nonzero(~np.isfinite(cdl) | (np.abs(cdl)>1e-7+1e-6*np.abs(h['saved_lg'][idx])))),
            hypothesis_hg_mean=moments(dh)['mean'], hypothesis_lg_mean=moments(dl)['mean'],
            saved_minus_sensor_mean=moments(t-h['sensor_temperature'][idx])['mean'],
            saved_minus_sensor_absmax=absmax(t-h['sensor_temperature'][idx])))
        if role == 'mip':
            # Freeze the HG-only hit set across hypotheses; no cuts on corrected MIP/energy.
            mask = np.isfinite(rh) & (rh < a.hg_max) & np.isfinite(f) & (f > 0)
            v = (rh-p.ph)*h['path_factor'][idx]
            stages = dict(raw=v, corrected=v*f, saved=h['saved_hg'][idx]/.305*p.mip*h['path_factor'][idx])
            factors = dict(raw=np.ones(len(t)), corrected=f, saved=f)
            if a.compare_models:
                stages.update(sign_flipped=v*(2-f), reciprocal=v/f)
                factors.update(sign_flipped=2-f, reciprocal=1/f)
            status = 'ok'
            if tb == -2147483648:
                status = 'missing_temperature_sensors'
            elif cell not in threshold:
                status = 'missing_threshold'
            low, high = math.nan, math.nan
            if status == 'ok' and np.any(mask):
                # threshold ADC is pedestal-subtracted (same convention as mip_closure).
                th, sigma, source = threshold[cell]
                # Common, conservative range for ALL stages, fully above the hardware turn-on
                # and fully below the HG cut after every transformation/path correction.
                production_factor = factor(t, cell)
                lower = max(10/float(np.min(production_factor[mask])), th+a.plateau_sigmas*sigma)
                low = max(lower*float(np.max(q[mask]*h['path_factor'][idx][mask])) for q in factors.values())
                high = min(a.fit_high_mip*p.mip, *(float(np.min(q[mask]*h['path_factor'][idx][mask]))*(a.hg_max-p.ph) for q in factors.values()))
                if not 0 <= low < high or high > 8*p.mip:
                    status = 'no_common_fit_interval'
            for stage, values in stages.items():
                name = f'f{fi}_r{run}_c{cell}_t{tb}_{stage}'
                result = dict(n=int(mask.sum()), mpv=math.nan, mpv_err=math.nan, valid=False, reason=status if status != 'ok' else 'empty_selection')
                if status == 'ok' and np.any(mask):
                    result = fit_spectrum(ROOT, rootout, name, values[mask], np.ones(mask.sum()), p.mip, low, high, a.min_fit_hits)
                # Stable schema even for missing-temperature and fit-failure rows.
                mips.append(dict(base, stage=stage, n=result['n'], fit_n=result.get('fit_effective_n',0),
                    fit_low=low, fit_high=high, mpv=result['mpv'], mpv_err=result['mpv_err'],
                    mpv_MIP=result['mpv']/p.mip, mpv_MIP_err=result['mpv_err']/p.mip,
                    valid=result['valid'], status=result['reason'], chi2=result.get('chi2',math.nan),
                    ndf=result.get('ndf',0), covariance=result.get('covariance',-1),
                    mip_constant=p.mip, mip_fallback=bool(p.mipFallback)))
        if role == 'response':
            good = np.isfinite(rh) & np.isfinite(rl) & np.isfinite(f) & (f>0) & (rh>=a.overlap_min) & (rh<a.overlap_max)
            good &= np.isfinite(h['saved_hg'][idx]) & np.isfinite(h['saved_lg'][idx])
            ab = np.zeros(len(rh), dtype=np.int64)
            ab[good] = np.floor((rh[good]-a.overlap_min)/a.adc_width).astype(np.int64)
            for b in np.unique(ab[good]):
                m = good & (ab == b)
                row = dict(base, adc_bin=int(b), raw_hg_low=a.overlap_min+b*a.adc_width,
                           raw_hg_high=min(a.overlap_max,a.overlap_min+(b+1)*a.adc_width),
                           mean_raw_hg=float(rh[m].mean()), n=int(m.sum()))
                for key, values in dict(saved_delta=h['saved_lg'][idx]-h['saved_hg'][idx],
                                        corrected_delta=en['lg']-en['hg'], raw_delta=en['lg_raw']-en['raw'],
                                        offset_before_delta=en['lg_offset_before']-en['hg']).items():
                    for stat, value in moments(values[m]).items():
                        row[key+'_'+stat] = value
                row['status'] = 'missing_temperature_sensors' if tb == -2147483648 else ('ok' if row['n']>=a.min_hits else 'low_statistics')
                residuals.append(row)
    return mips, residuals, closures, excluded


def pedestal_row(p, fi, role, param, min_hits):
    """Raw pedestal means remain available even if no calibration reference exists."""
    q={k:(int(v) if k in ('run','cellid','memory','temp_bin','n','missing_temperature') else float(v)) for k,v in p.items()}
    q.update(file_index=fi,kind=role,ped_present=bool(param.pedPresent),
             hg_residual=q['hg_mean']-param.ph if param.pedPresent else math.nan,
             lg_residual=q['lg_mean']-param.pl if param.pedPresent else math.nan)
    q['status']=('missing_pedestal' if not param.pedPresent else
                 'missing_temperature_sensors' if q['missing_temperature'] else
                 'ok' if q['n']>=min_hits else 'low_statistics')
    return q


def plots(ROOT, a, mips, residuals, peds):
    """ROOT graphs saved individually as well as eight-panel PNGs by physical channel."""
    if ROOT.gSystem.Load('libASImage') < 0:
        raise RuntimeError('Cannot load libASImage for fast PNG output')
    if not ROOT.gInterpreter.Declare('#include '+json.dumps(str(HERE.parent/'adc_energy/FastPng.hh'))):
        raise RuntimeError('Cannot load FastPng.hh')
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetLabelSize(.04, 'XYZ'); ROOT.gStyle.SetTitleSize(.04, 'XYZ')
    ROOT.gStyle.SetPadLeftMargin(.16); ROOT.gStyle.SetPadBottomMargin(.15)
    out = ROOT.TFile(str(a.output/'trends.root'), 'CREATE')
    fig = a.output/'figures';fig.mkdir()
    gm, gr, gp = (grouped(v,['cellid']) for v in (mips,residuals,peds))
    timings = []
    for (cell,) in sorted(set(gm)|set(gr)|set(gp)):
        canvas = ROOT.TCanvas(f'temperature_canvas_{cell}','',2400,1100); canvas.Divide(4,2)
        keep = []
        for panel, rows, xkey, ykey, ekey, groupkeys, title in [
            (1,gm.get((cell,),[]),'temperature','mpv_MIP','mpv_MIP_err',['stage','run'],'HG MIP MPV;Temperature [C];MPV [nominal MIP]'),
            (2,gm.get((cell,),[]),'run','mpv_MIP','mpv_MIP_err',['stage'],'HG MIP MPV;Run;MPV [nominal MIP]'),
            (3,gr.get((cell,),[]),'temperature','saved_delta_mean','saved_delta_sem',['run','adc_bin'],'Saved LG-HG: separate ADC bins;Temperature [C];LG-HG [MeV]'),
            (4,gr.get((cell,),[]),'run','saved_delta_mean','saved_delta_sem',['adc_bin','temp_bin'],'Saved LG-HG;Run;LG-HG [MeV]'),
            (5,gp.get((cell,),[]),'temperature','hg_residual','hg_sem',['run','memory'],'HitTag=0 HG: separate memory cells;Sensor temperature [C];HG - fixed pedestal [ADC]'),
            (6,gp.get((cell,),[]),'temperature','lg_residual','lg_sem',['run','memory'],'HitTag=0 LG: separate memory cells;Sensor temperature [C];LG - fixed pedestal [ADC]'),
            (7,gp.get((cell,),[]),'run','hg_residual','hg_sem',['memory','temp_bin'],'HitTag=0 HG;Run;HG - fixed pedestal [ADC]'),
            (8,gp.get((cell,),[]),'run','lg_residual','lg_sem',['memory','temp_bin'],'HitTag=0 LG;Run;LG - fixed pedestal [ADC]')]:
            canvas.cd(panel); ROOT.gPad.SetGrid()
            valid = [r for r in rows if r.get('valid',True) and (r.get('status','ok')=='ok' or (panel>=7 and r.get('status')=='missing_temperature_sensors' and r['n']>=a.min_hits))
                     and all(math.isfinite(float(r[k])) for k in (xkey,ykey,ekey))]
            multi=ROOT.TMultiGraph();multi.SetTitle(f'CellID {cell}: {title}')
            legend=ROOT.TLegend(.17,.70,.88,.91);legend.SetNColumns(2);legend.SetTextSize(.022)
            for i,(key,values) in enumerate(grouped(valid,groupkeys).items()):
                graph=ROOT.TGraphErrors(len(values));graph.SetName(f'c{cell}_panel{panel}_series{i}')
                label=', '.join(f'{k}={v}' for k,v in zip(groupkeys,key));graph.SetTitle(label)
                for j,r in enumerate(sorted(values,key=lambda v:v[xkey])):
                    graph.SetPoint(j,float(r[xkey]),float(r[ykey]));graph.SetPointError(j,0,float(r[ekey]))
                graph.SetMarkerColor(1+i%8);graph.SetLineColor(1+i%8);graph.SetMarkerStyle(20+i%5)
                multi.Add(graph,'P');ROOT.SetOwnership(graph,False)
                out.cd();graph.Write();keep.append(graph)
                if i<8:legend.AddEntry(graph,label,'p')
            if valid:
                multi.Draw('A');legend.Draw();keep.extend([multi,legend])
            else:
                label,xtitle,ytitle=title.split(';')
                ROOT.gPad.DrawFrame(0,0,1,1,f'CellID {cell}: {label} (insufficient);{xtitle};{ytitle}')
        png = fig/f'channel_{cell}.png'
        timing = ROOT.fastPng(canvas, str(png))
        # ROOT WriteImage reports some I/O errors without throwing an exception.
        # Confirm that a PNG was written before marking the analysis complete.
        with png.open('rb') as stream:
            if stream.read(8) != b'\x89PNG\r\n\x1a\n':
                raise RuntimeError(f'Invalid PNG output: {png}')
        timings.append(dict(cellid=cell, raster_seconds=float(timing.raster),
                            encode_seconds=float(timing.encode), bytes=png.stat().st_size))
        canvas.Close()
        if len(timings)%100 == 0:
            print(f'  PNG saved: {len(timings)}', flush=True)
    out.Close()
    write(a.output/'png_timings.tsv', timings,
          ['cellid','raster_seconds','encode_seconds','bytes'])
    print(f'PNG saved={len(timings)}, raster_seconds={sum(r["raster_seconds"] for r in timings):.3f}, '
          f'encode_seconds={sum(r["encode_seconds"] for r in timings):.3f}', flush=True)


def main():
    a = arguments(); rows = inputs(a)
    if a.output.exists():
        raise ValueError('Output exists; choose a new directory')
    if a.dry_run:
        for row in rows: print(row)
        print('Output:',a.output);return
    import ROOT
    ROOT.gROOT.SetBatch(True);ROOT.TH1.AddDirectory(False)
    ROOT.RDataFrame(1).Define('temperature_init', '0.').AsNumpy()
    # Force lazy graphics/function classes to initialize before declaring helpers.
    _classes = (ROOT.TF1, ROOT.TCanvas, ROOT.TGraphErrors)
    # Reuse the exact legacy sentinel/fallback policy, recording it as candidate provenance.
    ROOT.gInterpreter.Declare('#include '+json.dumps(str(HERE.parent/'calibration_residual/CalibrationAudit.hh')))
    cpp = ROOT.readCalibration(str(a.pedestal),str(a.mip),str(a.hl))
    params = {int(pair.first):pair.second for pair in cpp}
    threshold = read_threshold(ROOT,a.threshold) if a.threshold else {}
    with (HERE/'.build.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        subprocess.run(['make','-C',str(HERE),'all'],check=True)
    a.output.mkdir(parents=True)
    snapshot(a,rows);ROOT.writeCalibration(cpp,str(a.output/'calibration_channels.tsv'))
    load_fit_function(ROOT)
    _fit_function = ROOT.muonLangau
    rootout=ROOT.TFile(str(a.output/'mip_spectra.root'),'CREATE')
    if not rootout or rootout.IsZombie():raise ValueError('Cannot create spectra file')
    mips,residuals,closures,peds,counters,excluded=[],[],[],[],[],[]
    excluded_fields=['file_index','kind','run','cellid','temp_bin','temperature','n',
                     'reason','mip_constant','ped_present']
    for fi,r in enumerate(rows):
        part=a.output/f'file_{fi}';part.mkdir()
        command=[str(HERE/'CollectTemperature'),r['kind'],r['raw'],r['calib'],r['track'],str(part/'hits.root'),
                 str(part/'pedestal.tsv'),str(part/'counters.tsv'),a.channels,str(a.events_per_file),str(a.temp_width),str(a.max_angle)]
        print(f'File {fi}: {r["kind"]} {r["raw"]}',flush=True)
        with (part/'collector.log').open('w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        h=arrays(ROOT,part/'hits.root')
        mm,rr,cc,ee=process(ROOT,a,fi,r['kind'],h,params,threshold,rootout)
        mips.extend(mm);residuals.extend(rr);closures.extend(cc)
        excluded.extend(ee)
        # Flush exclusions after each file, so a later input error retains the audit.
        write(a.output/'excluded_channels.tsv',excluded,excluded_fields)
        if ee:
            print(f'  WARNING: excluded {len({v["cellid"] for v in ee})} channels / '
                  f'{sum(v["n"] for v in ee)} signal hits with invalid constants; '
                  'see excluded_channels.tsv (pedestal audit continues)',flush=True)
        del h
        for p in table(part/'pedestal.tsv'):
            peds.append(pedestal_row(p,fi,r['kind'],params[int(p['cellid'])],a.min_hits))
        counters.extend(dict(file_index=fi,kind=r['kind'],**c) for c in table(part/'counters.tsv'))
        print(f'  MIP fit rows={len(mm)}, residual bins={len(rr)}, closure bins={len(cc)}',flush=True)
    rootout.Close()
    write(a.output/'mip_fits.tsv',mips);write(a.output/'residual_bins.tsv',residuals)
    write(a.output/'closure.tsv',closures);write(a.output/'pedestal.tsv',peds);write(a.output/'counters.tsv',counters)
    trend_rows=[]
    for dataset,values,keys,metrics in [
        ('mip',[r for r in mips if r['valid']],['cellid','stage'],[('mpv_MIP','mpv_MIP_err')]),
        ('residual',[r for r in residuals if r['status']=='ok'],['cellid','adc_bin'],[(k+'_mean',k+'_sem') for k in ('saved_delta','corrected_delta','raw_delta','offset_before_delta')]),
        ('pedestal',[r for r in peds if r['status']=='ok'],['cellid','memory'],[('hg_residual','hg_sem'),('lg_residual','lg_sem')])]:
        for value,error in metrics:
            for r in trends(values,keys,value,error,a.reference_temperature,a.min_temp_span):
                row=dict(dataset=dataset,stage='',adc_bin='',memory='')
                row.update(r);trend_rows.append(row)
    # Different datasets have different keys; write their union in a stable order.
    fields=list(dict.fromkeys(k for r in trend_rows for k in r))
    write(a.output/'temperature_slopes.tsv',trend_rows,fields or None)
    if not a.no_plots:plots(ROOT,a,mips,residuals,peds)
    summary=dict(mip_inputs=sum(r['kind']=='mip' for r in rows),response_inputs=sum(r['kind']=='response' for r in rows),
                 valid_mip_fits=sum(r['valid'] for r in mips),
                 missing_sensor_signal_hits=sum(int(r['missing_sensor_hits']) for r in counters),
                 production_hg_mismatch=sum(r['production_hg_mismatch'] for r in closures),
                 production_lg_mismatch=sum(r['production_lg_mismatch'] for r in closures),
                 excluded_signal_channels=len({r['cellid'] for r in excluded}),
                 excluded_signal_hits=sum(r['n'] for r in excluded),
                 evaluated_signal_hits=sum(r['n'] for r in closures),
                 pedestal_bins_missing_reference=sum(not r['ped_present'] for r in peds),
                 signal_status=('no_valid_signal_constants' if excluded and not closures else
                                'evaluated_with_exclusions' if excluded else 'no_constant_exclusions'),
                 mip_temperature_provenance='unconfirmed' if a.mip_constant_temperature is None else 'user_declared',
                 interpretation='Check HG MIP slopes first, then fixed-ADC LG-HG residuals and per-memory pedestal. No data/MC pass/fail inference.')
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (a.output/'COMPLETE').touch();print('Saved',a.output)


if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError,OSError,subprocess.CalledProcessError) as e:
        print('ERROR:',e,file=sys.stderr);sys.exit(1)
