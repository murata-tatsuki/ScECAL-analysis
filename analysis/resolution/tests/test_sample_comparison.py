#!/usr/bin/env python3
"""Integration checks using synthetic SingleEnergy ROOTs (requires matching PyROOT).

Run after make MultiEnergyAnalysis. --channels also renders all 360 channel plots.
"""
import argparse
import json
import math
from pathlib import Path
import shlex
import subprocess
import tempfile

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ForceStyle(False)  # Inspect stored styles even if .rootlogon enables ForceStyle.
SOURCE = Path(__file__).resolve().parents[1]
ENERGIES = [0.5, 10, 20, 50, 100, 250]
LABELS = ["Calibration A (reference)", "Temperature correction", "MC alternative"]


def run(args, expected=None):
    result = subprocess.run([str(x) for x in args], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if expected is None:
        assert result.returncode == 0, result.stdout[-8000:]
    else:
        assert result.returncode != 0 and expected in result.stdout, result.stdout[-8000:]
    return result.stdout


def fixture(path, energy, sample):
    path.parent.mkdir(parents=True, exist_ok=True)
    output = ROOT.TFile(str(path), "RECREATE")
    output.mkdir("beforeEventCut").cd()
    names = [f"edep_{int(energy)}GeV", f"nhits_{int(energy)}GeV", "edep_1hit_10um", "edep_1hit_15um"]
    names += [f"{kind}_{int(energy)}GeV_layer{layer}" for layer in range(30) for kind in ("edep", "nhits")]
    for name in names:
        hist = ROOT.TH1F(name, name, 400, 0, 12000 if "1hit" not in name else 12)
        for value in range(100):
            hist.Fill(600 + value if "1hit" not in name else 0.1 + value * 0.005)
        hist.Write()
        hist.SetDirectory(0)
    for name in ("fit_gaus", "fit_gaus_nhit"):
        fit = ROOT.TF1(name, "gaus", 0, 12000)
        mean = energy * 60
        fit.SetParameters(100, mean, mean * math.sqrt((0.18 + sample * 0.02) ** 2 / energy + 0.01 ** 2))
        fit.SetParError(1, mean * 0.001)
        fit.SetParError(2, mean * 0.001)
        fit.Write()
    output.Close()


def add_channels(path):
    output = ROOT.TFile(str(path), "UPDATE")
    output.mkdir("10GeV").cd()
    canvas = ROOT.TCanvas("raw_fixture", "raw fixture", 600, 600)
    canvas.Divide(6, 6)
    histograms = []
    for channel in range(36):
        hist = ROOT.TH1F(f"fixture_{channel}", "", 200, 0, 2)
        hist.SetDirectory(0)
        if channel != 0:  # Exercise the empty-channel normalization path too.
            for value in range(100):
                hist.Fill(0.1 + value * 0.005)
        canvas.cd(channel + 1)
        hist.Draw()
        histograms.append(hist)
    for prefix in ("", "cut_"):
        for layer in range(30):
            for chip in range(6):
                for channel, hist in enumerate(histograms):
                    hist.SetName(f"edep_channel_{layer}_{chip}_{channel}")
                canvas.Write(f"{prefix}edep_Layer{layer}_Chip{chip}")
    canvas.Close()
    output.Close()


def legends(canvas):
    return [str(entry.GetLabel()) for obj in canvas.GetListOfPrimitives()
            if obj.InheritsFrom("TLegend") for entry in obj.GetListOfPrimitives()]



def check_fit_curves(pad):
    for obj in pad.GetListOfPrimitives():
        if obj.InheritsFrom("TPad"):
            check_fit_curves(obj)
        functions = ([obj] if obj.InheritsFrom("TF1") else
                     list(obj.GetListOfFunctions()) if obj.InheritsFrom("TGraph") else [])
        for fit in functions:
            if not fit.InheritsFrom("TF1") or not fit.GetName().startswith("res_func_low_"):
                continue
            for energy in (.5, 2, 30, 80, 100, 250):
                p = [fit.GetParameter(i) for i in range(fit.GetNpar())]
                expected = math.sqrt(p[0]**2/energy + (p[1]**2/energy**2+p[2]**2 if len(p)==3 else p[1]**2))
                assert abs(fit.Eval(energy)-expected)<1e-12, (fit.GetName(),energy)


def check_output(path, labels, generic=True):
    output = ROOT.TFile(str(path))
    for name in ("resolution", "resolution_with_noise"):
      for log in (False, True):
        canvas = output.Get(name + ("_loglog" if log else ""))
        assert canvas, name
        check_fit_curves(canvas)
        pad = canvas if generic else canvas.GetPad(1)
        assert pad.GetLogx() == int(log) and pad.GetLogy() == int(log)
        entries = legends(pad)
        assert all(any(entry.startswith(label) for entry in entries) for label in labels), entries
        graphs = [obj for obj in pad.GetListOfPrimitives() if obj.InheritsFrom("TGraphErrors")]
        assert len(graphs) == len(labels), len(graphs)
        if generic:
            assert len({graph.GetLineColor() for graph in graphs}) == len(labels)
        for sample, graph in enumerate(graphs):
            assert graph.GetN() == len(ENERGIES)
            axis = graph.GetXaxis()
            if log:
                assert 0 < axis.GetXmin() < min(ENERGIES) and axis.GetXmax() > max(ENERGIES)
                assert graph.GetMinimum() > 0
            else:
                assert axis.GetXmin() == 0 and axis.GetXmax() == 130
                assert graph.GetMinimum() == 0 and graph.GetMaximum() == .3
            for index, energy in enumerate(ENERGIES):
                assert graph.GetPointX(index) == energy
                expected = math.sqrt((0.18 + sample * 0.02) ** 2 / energy + 0.01 ** 2)
                assert abs(graph.GetPointY(index) - expected) < 1e-12
        if not generic:
            for index in (3, 4):
                assert canvas.GetPad(index).GetLogx() == int(log)
                assert canvas.GetPad(index).GetLogy() == int(log)
    output.Close()


def checks(work, channels):
    inputs = []
    for sample in range(3):
        files = [work / f"sample {sample}" / f"{energy:g}GeV_result.root" for energy in ENERGIES]
        for energy, path in zip(ENERGIES, files):
            fixture(path, energy, sample)
        inputs.append(files)
    output, figures = work / "comparison.root", work / "figures"
    config = work / "compare.sh"
    config.write_text("\n".join([
        f"multi_output={shlex.quote(str(output))}",
        f"multi_figures={shlex.quote(str(figures))}",
        "multi_options=(--skip-channel-plots)",
        *["multi_sample " + shlex.join([label, *map(str, reversed(files))]) for label, files in zip(LABELS, inputs)],
    ]))
    run(["bash", SOURCE / "run_samples.sh", "--keep-tail-events", "multi", config])
    check_output(output, LABELS)
    assert (figures / "resolution.png").is_file()
    assert (figures / "resolution_with_noise.png").is_file()
    assert (figures / "resolution_loglog.png").is_file()
    assert (figures / "resolution_with_noise_loglog.png").is_file()
    print("PASS: shell-selected three samples, labels with spaces, unordered and fractional energies, resolution values/styles", flush=True)

    binary = SOURCE / "MultiEnergyAnalysis"
    def command(groups, labels=LABELS):
        return [binary, work / "invalid.root", len(groups), len(groups[0]),
                *[path for group in groups for path in group], work / "invalid_figs", "--labels", *labels, "--skip-channel-plots", "--keep-tail-events"]
    mismatch = [*inputs[:2], [*inputs[2][:-1], work / "200GeV_result.root"]]
    run(command(mismatch), "Beam energies do not match for sample 3")
    duplicate = [*inputs[:2], [*inputs[2][:-1], inputs[2][0]]]
    run(command(duplicate), "Duplicate beam energy")
    run(command(inputs, LABELS[:2]), "Labels must be non-empty")
    missing = [[work / "absent" / file.name for file in inputs[0]], *inputs[1:]]
    run(command(missing), "cannot read")
    empty = work / "empty" / "0.5GeV_result.root"
    empty.parent.mkdir()
    ROOT.TFile(str(empty), "RECREATE").Close()
    run(command([[empty, *inputs[0][1:]], *inputs[1:]]), "Missing TH1F")
    print("PASS: reject third-sample energy mismatch, duplicate energies, missing labels/files/objects", flush=True)

    legacy = []
    for sample, directory in enumerate(("data", "simulation/threshold")):
        files = []
        for energy, source in zip(ENERGIES, inputs[sample]):
            dest = work / directory / f"{energy:g}GeV_20mm.root"
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.symlink_to(source)
            files.append(dest)
        legacy.extend(files)
    legacy_out = work / "legacy.root"
    run([binary, legacy_out, 2, len(ENERGIES), *legacy, work / "legacy_figures", "--skip-channel-plots", "--keep-tail-events"])
    check_output(legacy_out, ["data 20mm", "sim 20mm"], generic=False)
    print("PASS: legacy positional data/simulation interface and legends", flush=True)

    # Capture SingleEnergy invocations instead of running a costly detector analysis.
    recorder = work / "record.py"
    capture = work / "calls.jsonl"
    recorder.write_text("import json,sys\nwith open(sys.argv[1], 'a') as f: f.write(json.dumps(sys.argv[2:])+'\\n')\n")
    calib_dirs = [work / "calib PS", work / "calib SPS"]
    for index, base in enumerate(calib_dirs):
        directory = base / "10GeV"
        directory.mkdir(parents=True)
        (directory / f"run {index}.root").touch()
        (directory / "notes.txt").touch()
    single_config = work / "single.sh"
    single_config.write_text("\n".join([
        "single_runner=(" + shlex.join(["python3", str(recorder), str(capture)]) + ")",
        "single_sample " + shlex.join([str(work / "chosen output"), str(work / "chosen figures"), "20", "1", "*.root", *map(str, calib_dirs)]),
    ]))
    run(["bash", SOURCE / "run_samples.sh", "--dry-run", "single", single_config])
    assert not capture.exists()
    run(["bash", SOURCE / "run_samples.sh", "single", single_config])
    args = json.loads(capture.read_text())
    assert args[1] == str(work / "chosen output/10GeV_20mm.root")
    assert args[2:5] == ["1", "10GeV", "2"]
    assert all(path.endswith(".root") for path in args[5:7])
    assert args[7:10] == ["20", "1", str(work / "chosen figures")]
    assert args[10] == "--exclude-tail-events" and args[11] == "--tail-results"
    print("PASS: SingleEnergy output selection, merged calibration dirs, glob count, paths with spaces, dry run", flush=True)

    if channels:
        raw_input = inputs[0][1]
        add_channels(raw_input)
        raw_output = work / "channels.root"
        run([binary, raw_output, 3, 1, raw_input, raw_input, raw_input, work / "channel_figures", "--labels", *LABELS, "--keep-tail-events"])
        raw = ROOT.TFile(str(raw_output))
        directory = raw.GetDirectory("samples/10GeV/raw")
        assert directory
        for prefix in ("", "cut_"):
            for layer in range(30):
                for chip in range(6):
                    assert directory.GetKey(f"{prefix}edep_Layer{layer}_Chip{chip}")
            canvas = raw.Get(f"samples/10GeV/raw/{prefix}edep_Layer0_Chip0")
            assert legends(canvas) == LABELS
            histograms = [obj for obj in canvas.GetPad(2).GetListOfPrimitives() if obj.InheritsFrom("TH1")]
            assert len(histograms) == 3
            assert len({hist.GetLineColor() for hist in histograms}) == 3
            assert all(abs(hist.Integral() - 1) < 1e-6 for hist in histograms)
            empty_histograms = [obj for obj in canvas.GetPad(1).GetListOfPrimitives() if obj.InheritsFrom("TH1")]
            assert len(empty_histograms) == 3 and all(hist.Integral() == 0 for hist in empty_histograms)
        raw.Close()
        assert len(list((work / "channel_figures/samples/10GeV/raw").glob("*.png"))) == 360
        print("PASS: all 360 generic channel canvases/PNGs, three legends/colors, normalization and empty channels", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--channels", action="store_true")
    parser.add_argument("--work-dir", type=Path)
    options = parser.parse_args()
    if options.work_dir:
        options.work_dir.mkdir(parents=True, exist_ok=True)
        checks(options.work_dir.resolve(), options.channels)
    else:
        with tempfile.TemporaryDirectory(prefix="resolution-comparison-") as directory:
            checks(Path(directory), options.channels)
