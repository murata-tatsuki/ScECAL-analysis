"""Readable overview, detector map, and paired-channel plots; no refitting.

Usage: python3.12 visual_report.py RESULT_DIRECTORY
"""
import csv
import hashlib
import json
import math
from array import array
from pathlib import Path


def read_channels(path):
    with path.open() as f:
        channels = list(csv.DictReader(f, delimiter='\t'))
    for row in channels:
        for key, value in row.items():
            if key in ('cellid', 'pitch_um'):
                row[key] = int(value)
            elif key.endswith('_valid'):
                row[key] = value == 'True'
            elif not key.endswith('_source'):
                row[key] = float(value)
    return channels


def draw(ROOT, out, channels):
    ROOT.gROOT.SetBatch(True)
    style = ROOT.TStyle('mip_readable', 'MIP closure presentation')
    style.SetOptStat(0)
    style.SetOptTitle(0)
    style.SetCanvasColor(0)
    style.SetPadColor(0)
    style.SetFrameFillColor(0)
    style.SetTextFont(42)
    style.SetLabelFont(42, 'XYZ')
    style.SetTitleFont(42, 'XYZ')
    style.SetLabelSize(.035, 'XYZ')
    style.SetTitleSize(.04, 'XYZ')
    style.SetTitleOffset(1.2, 'X')
    style.SetTitleOffset(1.25, 'Y')
    style.SetPadLeftMargin(.12)
    style.SetPadRightMargin(.04)
    style.SetPadBottomMargin(.15)
    style.SetPadTopMargin(.16)
    style.SetLegendBorderSize(0)
    style.SetLegendFont(42)
    style.SetHistLineWidth(3)
    style.cd()
    alive = [style]
    valid = [r for r in channels if r['mc_adc_fit_valid']
             and math.isfinite(r['mc_adc_mpv_over_input'])]
    colors = {'individual': ROOT.TColor.GetColor('#0072B2'),
              'fallback': ROOT.TColor.GetColor('#D55E00')}

    def label(x, y, text, size=.035, ndc=True):
        t = ROOT.TLatex(x, y, text)
        t.SetTextFont(42); t.SetTextSize(size)
        if ndc:
            t.SetNDC()
        t.Draw(); alive.append(t)

    def save(canvas, name):
        for extension in ('png', 'pdf', 'svg'):
            canvas.SaveAs(str(out/f'{name}.{extension}'))

    canvas = ROOT.TCanvas('response_distribution', '', 1600, 750)
    canvas.Divide(2, 1)
    for pad, pitch in enumerate((10, 15), 1):
        canvas.cd(pad)
        rs = [r for r in valid if r['pitch_um'] == pitch]
        values = [r['mc_adc_mpv_over_input'] for r in rs]
        low = min(.7, math.floor(min(values, default=.7)/.05)*.05)
        high = max(1.4, (math.floor(max(values, default=1.4)/.05)+1)*.05)
        bins = max(1, round((high-low)/.05))
        hists = []
        for group in colors:
            h = ROOT.TH1D(f'distribution_{pitch}_{group}', '', bins, low, high)
            h.SetDirectory(0)
            for row in rs:
                if row['mc_legacy_source'] == group:
                    h.Fill(row['mc_adc_mpv_over_input'])
            h.SetLineColor(colors[group]); h.SetLineStyle(1 if group == 'individual' else 2)
            hists.append(h)
        maximum = max(1, max(h.GetMaximum() for h in hists))*1.55
        frame = ROOT.TH2D(f'distribution_frame_{pitch}', ';MC ADC MPV / input MIP;Channels / 0.05',
                          10, low, high, 10, 0, maximum)
        frame.SetDirectory(0); frame.Draw(); alive.append(frame)
        for h in hists:
            h.Draw('HIST SAME'); alive.append(h)
        target = ROOT.TLine(1, 0, 1, maximum*.73)
        target.SetLineStyle(3); target.SetLineWidth(2); target.Draw(); alive.append(target)
        label(.12, .93, f'{pitch} #mum  |  MC response', .05)
        label(.12, .875, f'{len(rs)} valid channels; reference = 1', .032)
        legend = ROOT.TLegend(.16, .64, .93, .81)
        legend.SetTextSize(.031); legend.SetFillStyle(0)
        for h, group in zip(hists, colors):
            vs = sorted(r['mc_adc_mpv_over_input'] for r in rs if r['mc_legacy_source'] == group)
            n = len(vs)
            median = (vs[(n-1)//2]+vs[n//2])/2 if n else math.nan
            legend.AddEntry(h, f'{group}: n={n}, median={median:.3f}', 'l')
        legend.Draw(); alive.append(legend)
        label(.12, .035, 'Groups: calibration policy (legacy). Each channel counts once.', .026)
    save(canvas, 'closure_summary')

    # A layer is 210 connected channels: chip 5 stops at channel 29.
    layers = max(30, max((r['cellid']//100000+1 for r in channels), default=30))
    canvas = ROOT.TCanvas('response_map', '', 1700, 1050)
    canvas.SetLeftMargin(.09); canvas.SetRightMargin(.14)
    canvas.SetTopMargin(.14); canvas.SetBottomMargin(.18)
    ROOT.TColor.CreateGradientColorTable(3, array('d', [0, .5, 1]),
        array('d', [.13, 1, .70]), array('d', [.40, 1, .09]), array('d', [.67, 1, .17]), 101)
    style.SetNumberContours(101)
    frame = ROOT.TH2D('channel_map_frame', ';Channel within layer (36 #times chip + channel);Layer',
                      210, -.5, 209.5, layers, -.5, layers-.5)
    frame.SetDirectory(0)
    for x in range(1, 211):
        for y in range(1, layers+1):
            frame.SetBinContent(x, y, 1)
    frame.SetMinimum(.5); frame.SetMaximum(1.5)
    frame.GetZaxis().SetTitle('MC ADC MPV / input MIP')
    frame.GetZaxis().SetTitleOffset(1.25)
    for axis in (frame.GetXaxis(), frame.GetYaxis(), frame.GetZaxis()):
        axis.SetTitleSize(.028); axis.SetLabelSize(.023)
    frame.GetYaxis().SetTitleOffset(1.2)
    frame.GetYaxis().SetTickLength(0)
    for layer in range(layers):
        frame.GetYaxis().SetBinLabel(layer+1, str(layer))
    frame.GetYaxis().SetLabelSize(.020)
    frame.GetXaxis().SetNdivisions(510)
    frame.Draw('COLZ'); alive.append(frame)
    background = ROOT.TBox(-.5, -.5, 209.5, layers-.5)
    background.SetFillColor(ROOT.TColor.GetColor('#E0E0E0'))
    background.SetLineColor(0); background.Draw(); alive.append(background)
    for row in valid:
        layer, rest = divmod(row['cellid'], 100000)
        chip, ch = divmod(rest, 10000)
        if not (0 <= chip < 6 and 0 <= ch < (30 if chip == 5 else 36)):
            raise ValueError(f"Invalid connected CellID {row['cellid']}")
        x = 36*chip+ch
        index = round(min(1, max(0, row['mc_adc_mpv_over_input']-.5))*100)
        box = ROOT.TBox(x-.5, layer-.5, x+.5, layer+.5)
        box.SetFillColor(ROOT.TColor.GetColorPalette(index)); box.SetLineWidth(0)
        box.Draw(); alive.append(box)
    for chip in range(1, 6):
        x = chip*36-.5
        line = ROOT.TLine(x, -.5, x, layers-.5)
        line.SetLineColor(ROOT.kWhite); line.SetLineWidth(2); line.Draw(); alive.append(line)
    for boundary in (3.5, 27.5):
        line = ROOT.TLine(-.5, boundary, 209.5, boundary)
        line.SetLineColor(ROOT.kBlack); line.SetLineStyle(2); line.Draw(); alive.append(line)
    ROOT.gPad.RedrawAxis()
    label(.08, .955, 'MC response by layer and connected channel', .035)
    label(.08, .905, 'Blue: below 1    White: near 1    Red: above 1    Gray: no valid MC ADC fit', .027)
    clipped = sum(not .5 <= r['mc_adc_mpv_over_input'] <= 1.5 for r in valid)
    label(.09, .050, f'Chip 0-4: 36 channels; chip 5: 30. Color saturates at 0.5 / 1.5 ({clipped} channels outside).', .023)
    label(.09, .020, 'Layers 4-27: 10 #mum; other layers: 15 #mum. Dashed lines mark the sensor-type boundaries.', .023)
    save(canvas, 'mc_channel_map')

    paired = sorted((r for r in channels if r['mc_adc_fit_valid'] and r['data_fit_valid']
                     and r['mc_fit_valid'] and math.isfinite(r['clean_data_over_mc'])),
                    key=lambda r: (r['pitch_um'], r['cellid']))
    pages = [paired[i:i+20] for i in range(0, len(paired), 20)] or [[]]
    for page_index, page_rows in enumerate(pages):
        canvas = ROOT.TCanvas(f'paired_channels_{page_index}', '', 1450, max(650, 220+38*len(page_rows)))
        canvas.SetLeftMargin(.30); canvas.SetRightMargin(.06)
        coordinates = []
        for r in page_rows:
            d, m = r['data_clean_mip_mpv'], r['mc_clean_mip_mpv']
            ed, em = r['data_clean_mip_mpv_err'], r['mc_clean_mip_mpv_err']
            err = math.sqrt((ed/m)**2+(d*em/m**2)**2)
            coordinates.append((r['clean_data_over_mc'], err))
        low = min([.8]+[v-e for v, e in coordinates]); high = max([1.2]+[v+e for v, e in coordinates])
        margin = .08*(high-low)
        n = max(1, len(page_rows))
        frame = ROOT.TH2D('pairs_frame', ';Calibrated MIP MPV: data / MC (clean policy);',
                         10, low-margin, high+margin, n, -.5, n-.5)
        frame.SetDirectory(0); frame.GetYaxis().SetLabelSize(.045); frame.GetYaxis().SetTickLength(0)
        for i, r in enumerate(page_rows):
            layer, rest = divmod(r['cellid'], 100000); chip, ch = divmod(rest, 10000)
            frame.GetYaxis().SetBinLabel(n-i, f"{r['pitch_um']} #mum   L{layer:02d} / chip {chip} / ch {ch:02d}")
        frame.Draw(); alive.append(frame)
        line = ROOT.TLine(1, -.5, 1, n-.5)
        line.SetLineStyle(2); line.SetLineColor(ROOT.kGray+2); line.Draw(); alive.append(line)
        for i, (r, (v, e)) in enumerate(zip(page_rows, coordinates)):
            graph = ROOT.TGraphErrors(1)
            graph.SetPoint(0, v, n-1-i); graph.SetPointError(0, e, 0)
            graph.SetMarkerStyle(20 if r['pitch_um'] == 10 else 21)
            graph.SetMarkerColor(colors['individual'] if r['pitch_um'] == 10 else colors['fallback'])
            graph.SetLineColor(graph.GetMarkerColor()); graph.SetMarkerSize(1.2)
            graph.Draw('P SAME'); alive.append(graph)
        label(.08, .95, f'Data / MC closure: {len(paired)} channels  |  page {page_index+1}/{len(pages)}', .035)
        label(.08, .90, 'One row per channel; reference = 1. Only channels passing both fits are shown.', .026)
        label(.08, .05, 'Bars: propagated fit statistical errors (data / MC assumed independent).', .022)
        label(.08, .02, 'Calibration and selection systematics are not included. Sparse coverage is not a detector-wide validation.', .022)
        if not page_rows:
            label(.35, .5, 'No valid paired channels', .04)
        if len(pages) == 1:
            save(canvas, 'data_mc_pairs')
        else:
            for extension in ('png', 'svg'):
                canvas.SaveAs(str(out/f'data_mc_pairs_page_{page_index+1:03d}.{extension}'))
                if page_index == 0:
                    canvas.SaveAs(str(out/f'data_mc_pairs.{extension}'))
            suffix = '(' if page_index == 0 else ')' if page_index == len(pages)-1 else ''
            canvas.Print(str(out/'data_mc_pairs.pdf')+suffix)
    manifest = dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    input_sha256=hashlib.sha256((out/'channel_closure.tsv').read_bytes()).hexdigest(),
                    valid_mc_channels=len(valid), paired_channels=len(paired), paired_pages=len(pages),
                    plots=['closure_summary', 'mc_channel_map', 'data_mc_pairs'],
                    formats=['png', 'pdf', 'svg'], refitted=False)
    (out/'visualization.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    import argparse
    import ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result_directory', type=Path)
    args = parser.parse_args()
    draw(ROOT, args.result_directory, read_channels(args.result_directory/'channel_closure.tsv'))
