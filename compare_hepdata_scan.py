#!/usr/bin/env python3
"""Compare a saved profile scan with HEPData's observed 4-lepton NSBI curve.

Reads downloaded HEPData table JSON (record/data endpoint) and optionally
cross-checks another table and the CSV download. No fit or model is changed.
"""

import argparse
from array import array
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def series_points(table, series="Obs NSBI"):
    if table["x_count"] != 1 or "4l only" not in table["name"]:
        raise ValueError("Expected a one-dimensional, four-lepton-only table")
    labels = [header["name"] for header in table["headers"][1:]]
    group = labels.index(series)
    points = {}
    for row in table["values"]:
        mu = float(row["x"][0]["value"])
        value = float(next(y["value"] for y in row["y"] if y["group"] == group))
        if mu in points or not math.isfinite(mu) or not math.isfinite(value):
            raise ValueError("Duplicate or non-finite reference point")
        points[mu] = value
    return points


def csv_series_points(path, series="Obs NSBI"):
    points = {}
    selected = False
    with path.open(newline="") as stream:
        for row in csv.reader(stream):
            if not row or row[0].startswith("#"):
                continue
            if len(row) == 2 and row[1] == series:
                selected = True
                continue
            if selected:
                try:
                    mu, value = map(float, row)
                except ValueError:
                    break  # The next series starts here.
                if mu in points:
                    raise ValueError("Duplicate CSV point")
                points[mu] = value
    if not points:
        raise ValueError(f"CSV does not contain the {series} series")
    return points


def make_plot(scan, rows, summary, output):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    low, high = summary["comparison_range"]
    canvas = ROOT.TCanvas("hepdata_comparison", "Observed scan comparison", 1000, 850)
    upper = ROOT.TPad("scan_overlay", "", 0, 0.31, 1, 1)
    lower = ROOT.TPad("scan_residual", "", 0, 0, 1, 0.31)
    for pad in (upper, lower):
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.04)
        pad.SetTicks(1, 1)
        pad.Draw()
    upper.SetBottomMargin(0.025)
    upper.SetTopMargin(0.07)
    lower.SetTopMargin(0.035)
    lower.SetBottomMargin(0.31)

    def graph(name, points, color, marker=20):
        g = ROOT.TGraph(len(points), array("d", [p[0] for p in points]),
                        array("d", [p[1] for p in points]))
        g.SetName(name)
        g.SetLineColor(color)
        g.SetMarkerColor(color)
        g.SetMarkerStyle(marker)
        g.SetMarkerSize(0.9)
        g.SetLineWidth(3)
        return g

    upper.cd()
    ymax = 1.12 * max(max(scan.values()), max(r["q_hepdata"] for r in rows))
    frame = upper.DrawFrame(low, 0, high, ymax)
    frame.GetYaxis().SetTitle("t_{#mu} = -2 ln #lambda(#mu)")
    frame.GetYaxis().SetTitleSize(0.052)
    frame.GetYaxis().SetTitleOffset(1.20)
    frame.GetYaxis().SetLabelSize(0.043)
    frame.GetXaxis().SetLabelSize(0)
    local = graph("workspace_scan", sorted(scan.items()), ROOT.kAzure + 2)
    reference = graph("hepdata_observed", [(r["mu"], r["q_hepdata"]) for r in rows], ROOT.kBlack)
    local.Draw("L SAME")
    reference.Draw("P SAME")
    legend = ROOT.TLegend(0.18, 0.68, 0.68, 0.84)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.036)
    legend.AddEntry(local, "Supplied workspace: RooFit profile", "l")
    legend.AddEntry(reference, "HEPData: observed NSBI, stat+syst", "p")
    legend.Draw()
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextSize(0.042)
    label.DrawLatex(0.18, 0.88, "H* #rightarrow ZZ #rightarrow 4l")
    label.SetTextSize(0.032)
    label.DrawLatex(0.18, 0.62, f"#hat{{#mu}} = {summary['mu_hat']:.5f}; 127 systematics + 3 qqZZ factors profiled")
    label.DrawLatex(0.18, 0.56, "HEPData 156805 v1, Table 2 / Figure 11(b)")

    lower.cd()
    residual_limit = max(0.7, 1.25 * 1000 * summary["max_abs_delta_q"])
    residual_frame = lower.DrawFrame(low, -residual_limit, high, residual_limit)
    residual_frame.GetXaxis().SetTitle("#mu_{off-shell}")
    residual_frame.GetXaxis().SetTitleSize(0.115)
    residual_frame.GetXaxis().SetLabelSize(0.09)
    residual_frame.GetXaxis().SetTitleOffset(1.05)
    residual_frame.GetYaxis().SetTitle("10^{3} #Delta t_{#mu}")
    residual_frame.GetYaxis().SetTitleSize(0.10)
    residual_frame.GetYaxis().SetTitleOffset(0.61)
    residual_frame.GetYaxis().SetLabelSize(0.085)
    residual_frame.GetYaxis().SetNdivisions(505)
    rounding = ROOT.TBox(low, -0.5, high, 0.5)
    rounding.SetFillColor(ROOT.kGray)
    rounding.SetFillStyle(1001)
    rounding.Draw()
    zero = ROOT.TLine(low, 0, high, 0)
    zero.SetLineStyle(2)
    zero.Draw()
    residual = graph("local_minus_hepdata", [(r["mu"], 1000 * r["delta_q"]) for r in rows], ROOT.kAzure + 2)
    residual.Draw("P SAME")
    label.SetTextSize(0.068)
    label.DrawLatex(0.17, 0.87, "Local - HEPData; gray: #pm0.0005 rounding interval")
    lower.RedrawAxis()
    canvas.cd()
    for extension in ("png", "pdf"):
        canvas.SaveAs(str(output / f"observed_scan_comparison.{extension}"))
    root_file = ROOT.TFile(str(output / "observed_scan_comparison.root"), "RECREATE")
    for obj in (local, reference, residual, canvas):
        obj.Write()
    root_file.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan", type=Path, required=True)
    parser.add_argument("--bestfit", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--crosscheck-table", type=Path)
    parser.add_argument("--reference-csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    table = json.loads(args.reference.read_text())
    if table["doi"] != "10.17182/hepdata.156805.v1/t2":
        raise ValueError("This comparison expects HEPData 156805 v1 Table 2")
    reference = series_points(table)
    checks = {}
    if args.crosscheck_table:
        other = json.loads(args.crosscheck_table.read_text())
        if series_points(other) != reference:
            raise ValueError("The duplicate observed series disagree")
        checks["duplicate_table"] = other["doi"]
    if args.reference_csv:
        if csv_series_points(args.reference_csv) != reference:
            raise ValueError("The JSON and CSV observed values disagree")
        checks["csv_matches_json"] = True
    scan = {}
    with args.scan.open(newline="") as stream:
        for row in csv.DictReader(stream):
            mu, q = float(row["mu"]), float(row["q"])
            if int(row["migrad_status"]) or not math.isfinite(q) or mu in scan:
                raise ValueError("Invalid, unconverged, or duplicate local scan point")
            scan[mu] = q
    # Both published and local grids contain exact decimal coordinates. Compare
    # only shared coordinates: no interpolation, rescaling, or minimum shift.
    common = sorted(scan.keys() & reference.keys())
    if len(common) < 2 or common[0] != min(scan) or common[-1] != max(scan):
        raise ValueError("Reference and local scan must share range endpoints")
    rows = [{"mu": mu, "q_local": scan[mu], "q_hepdata": reference[mu],
             "delta_q": scan[mu] - reference[mu]} for mu in common]
    largest = max(rows, key=lambda row: abs(row["delta_q"]))
    best = json.loads(args.bestfit.read_text())
    inputs = {name: {"path": str(path.resolve()), "sha256": digest(path)}
              for name, path in vars(args).items() if name != "output" and path}
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "reference_doi": table["doi"],
        "reference_url": "https://www.hepdata.net/record/156805?table="
                         "%24%5Cmu_%7B%5Cmathrm%7Boff-shell%7D%7D%24%20NSBI%20"
                         "stat%2Bsyst%20vs%20stat-only%2C%204l%20only",
        "reference_data_url": "https://www.hepdata.net/record/data/156805/1759753/1",
        "reference_series": "Obs NSBI",
        "reference_figure": table["location"],
        "reference_points": len(reference),
        "reference_range": [min(reference), max(reference)],
        "local_points": len(scan),
        "comparison_points": len(rows),
        "comparison_range": [common[0], common[-1]],
        "definition": "q = t_mu = 2 * (NLL_profile(mu) - NLL_global_min)",
        "method": "Exact shared mu coordinates; no interpolation or curve shifts",
        "mu_hat": best["mu_hat"],
        "reference_grid_minimum": {"mu": min(reference, key=reference.get),
                                   "q": min(reference.values())},
        "max_abs_delta_q": abs(largest["delta_q"]),
        "largest_residual": largest,
        "rms_delta_q": math.sqrt(sum(r["delta_q"]**2 for r in rows) / len(rows)),
        "mean_delta_q": sum(r["delta_q"] for r in rows) / len(rows),
        "reference_rounding_half_step": 0.0005,
        "points_within_rounding_half_step": sum(abs(r["delta_q"]) <= 0.0005 for r in rows),
        "max_excess_over_rounding_half_step": max(0, abs(largest["delta_q"]) - 0.0005),
        "rounding_note": "Three decimal places in the published values; this is not a statistical uncertainty",
        "crosschecks": checks,
        "inputs": inputs,
    }
    with (args.output / "comparison_points.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_plot(scan, rows, summary, args.output)
    print(json.dumps({k: v for k, v in summary.items() if k != "inputs"}, indent=2))


if __name__ == "__main__":
    main()
