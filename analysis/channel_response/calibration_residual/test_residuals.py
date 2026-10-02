"""Synthetic end-to-end regression checks; requires PyROOT only for `make check`."""
from array import array
import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import ROOT

ROOT.gROOT.SetBatch(True)
EXE = str(Path(sys.argv[1]).resolve())


def table(path):
    with open(path) as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_constants(base):
    f = ROOT.TFile(str(base / "ped.root"), "RECREATE")
    t = ROOT.TTree("ChnLevel", "")
    vectors = {"CellID": ROOT.std.vector("int")(),
               "PedHighMean": ROOT.std.vector("double")(),
               "PedLowMean": ROOT.std.vector("double")()}
    for name, vector in vectors.items():
        t.Branch(name, vector)
    for cell in [0, 420008, 420009]:
        vectors["CellID"].push_back(cell)
        vectors["PedHighMean"].push_back(500)
        vectors["PedLowMean"].push_back(400)
    t.Fill()
    t.Write()
    f.Close()
    for filename, treename, spec, rows in [
        ("mip.root", "MIP_Fit", [("CellID", "i"), ("LandauMPV", "d"), ("ChiSquare", "d"), ("NDF", "i")],
         [(0, 61, 1, 10), (420008, 61, 1, 10), (420009, 100, 50, 10)]),
        ("hl.root", "InterCalib", [("CellID", "i"), ("Slope", "d"), ("Intercept", "d")],
         [(0, .025, 0), (420008, .025, 0), (420009, .1, 0)]),
    ]:
        f = ROOT.TFile(str(base / filename), "RECREATE")
        t = ROOT.TTree(treename, "")
        holders = {name: array(kind, [0]) for name, kind in spec}
        for name, kind in spec:
            t.Branch(name, holders[name], f"{name}/{'I' if kind == 'i' else 'D'}")
        for row in rows:
            for (name, _), value in zip(spec, row):
                holders[name][0] = value
            t.Fill()
        t.Write()
        f.Close()


def write_events(path, events, calibrated=False, bad_hit=False, repeated_keys=False):
    f = ROOT.TFile(str(path), "RECREATE")
    t = ROOT.TTree("Calib_Hit" if calibrated else "Raw_Hit", "")
    scalar_names = ["Run_Num", "Event_Time", "Event_Num" if calibrated else "TriggerID"]
    scalars = {name: array("i", [0]) for name in scalar_names}
    for name, value in scalars.items():
        t.Branch(name, value, f"{name}/I")
    spec = {"CellID": "int"}
    spec.update({name: "double" for name in
                 (["Hit_HG_Energy", "Hit_LG_Energy", "Hit_Energy", "NewTemperature"] if calibrated else ["HG_Charge", "LG_Charge"])})
    if not calibrated:
        spec["HitTag"] = "int"
    vectors = {name: ROOT.std.vector(kind)() for name, kind in spec.items()}
    for name, vector in vectors.items():
        t.Branch(name, vector)
    for i, run, temperature, hg in events:
        scalars["Run_Num"][0] = run
        scalars["Event_Time"][0] = 1000 if repeated_keys else 1000 + i
        scalars[scalar_names[2]][0] = 1 if repeated_keys else i
        for vector in vectors.values():
            vector.clear()
        # Two memory cells per physical channel, reversed hit order in calibration.
        for memory in ([1, 0] if calibrated else [0, 1]):
            cell = 420008 + memory * 100
            correction = 1 - (temperature - 20) * 1.6 / 135
            eh = .005 * (hg - 500) * correction
            residual = .75 + .0001 * (hg - 2600) + (run - 10) * .1 + (temperature - 20) * .2
            el = eh + residual
            values = {"CellID": cell}
            if calibrated:
                values.update(Hit_HG_Energy=eh, Hit_LG_Energy=el,
                              Hit_Energy=eh if hg < 2600 else el, NewTemperature=temperature)
                if bad_hit and i == 0:
                    values["CellID"] = 420009
            else:
                values.update(HitTag=1, HG_Charge=hg, LG_Charge=400 + el / (.005 * correction * 40))
            for name, value in values.items():
                vectors[name].push_back(value)
        t.Fill()
    t.Write()
    f.Close()


def run(base, name, raw, cal, succeeds=True, plots=0, extra_args=()):
    manifest = base / f"{name}.tsv"
    manifest.write_text(f"{raw}\t{cal}\n")
    out = base / name
    args = [EXE, "--manifest", str(manifest), "--pedestal", str(base / "ped.root"),
            "--mip", str(base / "mip.root"), "--hl", str(base / "hl.root"),
            "--out", str(out), "--sample", "data", "--channels", "420008,420009",
            "--min-hits", "5", "--plots", str(plots)]
    args.extend(extra_args)
    result = subprocess.run(args, capture_output=True, text=True)
    assert (result.returncode == 0) == succeeds, result.stdout + result.stderr
    if succeeds:
        assert (out / "residuals.root").exists()
    else:
        assert not (out / "residuals.root").exists()
    return out, result.stderr


with tempfile.TemporaryDirectory(prefix="scecal_residual_test_") as folder:
    base = Path(folder)
    write_constants(base)
    events = []
    for run_number in [10, 11]:
        for temp in [20.0, 21.0]:
            for hg in range(2350, 2850, 5):
                events.append((len(events), run_number, temp, hg))
    write_events(base / "raw.root", events)
    write_events(base / "cal.root", events[::-1], calibrated=True)
    out, _ = run(base, "ordered_by_key", base / "raw.root", base / "cal.root", plots=1)
    rows = table(out / "channels.tsv")
    assert len(rows) == 4  # run x temperature, memory cells aggregated only after matching.
    for row in rows:
        expected = .75 + (int(row["run"]) - 10) * .1 + (float(row["temp_bin_low_C"]) - 20) * .2
        assert row["switch_status"] == "ok"
        assert int(row["hits"]) == 200
        for field in ["selected_jump_OLS_MeV", "selected_jump_median_bins_MeV", "paired_delta_at_switch_OLS_MeV"]:
            assert math.isclose(float(row[field]), expected, abs_tol=1e-10), (field, row[field], expected)
        assert int(row["selection_mismatch"]) == 0
        assert float(row["hg_closure_abs_max"]) < 1e-10
        assert float(row["lg_closure_abs_max"]) < 1e-10
    b = next(x for x in table(out / "residual_bins.tsv")
             if x["run"] == "10" and x["temp_bin_low_C"] == "20" and x["adc_bin_low"] == "1900")
    # HG 2400..2445 (5 ADC), each twice. Exact quantile interpolation, not histogram centers.
    values = sorted(.75 + .0001 * (h - 2600) for h in range(2400, 2450, 5) for _ in range(2))
    def quantile(p):
        rank = (len(values) - 1) * p
        k = int(rank)
        return values[k] + (values[min(k + 1, len(values) - 1)] - values[k]) * (rank - k)
    for column, p in [("q16_MeV", .16), ("median_MeV", .5), ("q84_MeV", .84)]:
        assert math.isclose(float(b[column]), quantile(p), abs_tol=1e-12)
    audit = {x["cellid"]: x for x in table(out / "calibration_channels.tsv")}
    assert audit["420009"]["mip_fallback"] == "1"
    assert audit["420009"]["hl_gain_fallback"] == "1"
    assert audit["420009"]["hl_offset_fallback"] == "0"
    assert audit["420009"]["hl_offset_unset_zero"] == "1"
    f = ROOT.TFile(str(out / "residuals.root"))
    assert f.Get("file_0/run_10/temperature_bin_40/channel_420008/delta_mean_sem")
    f.Close()
    assert len(list((out / "figures").rglob("*.png"))) == 4
    assert table(out / "files.tsv")[0]["matching"] == "unique_event_key"
    # Duplicate keys may use entry order only after the full sequence is validated.
    write_events(base / "dupraw.root", events, repeated_keys=True)
    write_events(base / "dupcal.root", events, calibrated=True, repeated_keys=True)
    dup, _ = run(base, "repeated", base / "dupraw.root", base / "dupcal.root")
    assert table(dup / "files.tsv")[0]["matching"] == "verified_entry_order"
    write_events(base / "dupcal_bad.root", events[::-1], calibrated=True, repeated_keys=True)
    _, error = run(base, "ambiguous", base / "dupraw.root", base / "dupcal_bad.root", succeeds=False)
    assert "Ambiguous duplicate" in error
    write_events(base / "badcal.root", events, calibrated=True, bad_hit=True)
    _, error = run(base, "missing_hit", base / "raw.root", base / "badcal.root", succeeds=False)
    assert "missing raw HitTag=1" in error
    low = [e for e in events if e[3] < 2600]
    write_events(base / "lowraw.root", low)
    write_events(base / "lowcal.root", low, calibrated=True)
    low_out, _ = run(base, "one_sided", base / "lowraw.root", base / "lowcal.root")
    assert all(x["switch_status"] == "insufficient_support" and math.isnan(float(x["selected_jump_OLS_MeV"])) for x in table(low_out / "channels.tsv"))
    # Total hits exceed min-hits, but each ADC bin is below it: no drawable points.
    sparse = [(i, 12, 20.0, 800 + i * 100) for i in range(10)]
    write_events(base / "sparseraw.root", sparse)
    write_events(base / "sparsecal.root", sparse, calibrated=True)
    sparse_out, _ = run(base, "empty_plot", base / "sparseraw.root", base / "sparsecal.root", plots=1)
    assert not (sparse_out / "figures").exists()
    assert int(table(sparse_out / "channels.tsv")[0]["hits"]) == 20
    assert len(table(sparse_out / "residual_bins.tsv")) == 10
    f = ROOT.TFile(str(sparse_out / "residuals.root"))
    assert f.Get("file_0/run_12/temperature_bin_40/channel_420008/delta_mean_sem").GetN() == 0
    f.Close()
    # Residual points alone must still produce a PNG when the switch window is empty.
    residual_only = [(i, 12, 20.0, 800 + i) for i in range(5)]
    write_events(base / "residualraw.root", residual_only)
    write_events(base / "residualcal.root", residual_only, calibrated=True)
    residual_out, _ = run(base, "residual_only", base / "residualraw.root", base / "residualcal.root", plots=1)
    assert len(list((residual_out / "figures").rglob("*.png"))) == 1
    assert not table(residual_out / "switch_bins.tsv")
    # Switch points alone must still produce PNGs when the residual display range excludes all hits.
    switch_out, _ = run(base, "switch_only", base / "raw.root", base / "cal.root", plots=1,
                        extra_args=("--adc-max", "100"))
    assert not table(switch_out / "residual_bins.tsv")
    assert len(list((switch_out / "figures").rglob("*.png"))) == 4
print("PASS: run/temperature separation, exact quantiles, injected switch jumps, saved selection, fallback audit, joins, missing hits, sparse support, ROOT/PNG output, empty PNG suppression")
