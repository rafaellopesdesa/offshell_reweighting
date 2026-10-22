#!/usr/bin/env python3
"""Compare two weighted Asimov profiles with the published expected NSBI scan."""

import argparse
from array import array
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path

from compare_hepdata_scan import series_points, csv_series_points, digest
from root_workspace_fit import require_minimum, write_json


def metrics(rows):
    largest = max(rows, key=lambda row: abs(row["delta_q"]))
    return {"points": len(rows), "max_abs_delta_q": abs(largest["delta_q"]),
            "largest_residual": largest,
            "rms_delta_q": math.sqrt(sum(r["delta_q"]**2 for r in rows)/len(rows)),
            "mean_delta_q": sum(r["delta_q"] for r in rows)/len(rows)}


def make_plot(scans, bestfits, reference, comparisons, output, stem="expected_scan_comparison"):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    canvas = ROOT.TCanvas(stem, "Asimov profile comparison", 1100, 900)
    upper = ROOT.TPad(stem + "_profiles", "", 0, .31, 1, 1)
    lower = ROOT.TPad(stem + "_residuals", "", 0, 0, 1, .31)
    for pad in (upper, lower):
        pad.SetLeftMargin(.13)
        pad.SetRightMargin(.04)
        pad.SetTicks(1, 1)
        pad.Draw()
    upper.SetTopMargin(.07)
    upper.SetBottomMargin(.025)
    lower.SetTopMargin(.035)
    lower.SetBottomMargin(.30)
    objects = []

    def graph(name, pairs, color, marker):
        pairs = sorted(pairs)
        obj = ROOT.TGraph(len(pairs), array("d", [x for x,y in pairs]), array("d", [y for x,y in pairs]))
        obj.SetName(name)
        obj.SetLineColor(color)
        obj.SetMarkerColor(color)
        obj.SetLineWidth(3)
        obj.SetMarkerStyle(marker)
        obj.SetMarkerSize(.65)
        objects.append(obj)
        return obj

    upper.cd()
    published = {x:y for x,y in reference.items() if 0 <= x <= 1.5}
    ymax = 1.2*max(list(published.values()) + [y for scan in scans.values() for y in scan.values()])
    frame = upper.DrawFrame(0, 0, 1.5, ymax)
    frame.GetYaxis().SetTitle("t_{#mu} = -2 ln #lambda(#mu)")
    frame.GetYaxis().SetTitleSize(.052)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.GetYaxis().SetLabelSize(.044)
    frame.GetXaxis().SetLabelSize(0)
    ref = graph("hepdata_expected", published.items(), ROOT.kBlack, 20)
    ref.SetLineWidth(2)
    ref.Draw("LP SAME")
    colors = {4096: ROOT.TColor.GetColor("#2166ac"), 8192: ROOT.TColor.GetColor("#d6604d")}
    legend = ROOT.TLegend(.39, .67, .93, .87)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(.036)
    legend.AddEntry(ref, "HEPData: expected NSBI, stat+syst", "lp")
    for size, scan in scans.items():
        pairs = list(scan.items())
        mu_hat = bestfits[size]["mu_hat"]
        if min(scan) <= mu_hat <= max(scan):
            pairs.append((mu_hat, 0))
        local = graph(f"asimov_{size}", pairs, colors[size], 24 if size == 4096 else 25)
        local.Draw("LP SAME")
        legend.AddEntry(local, f"{size} reference events; #hat{{#mu}} = {mu_hat:.4f}", "lp")
    legend.Draw()
    text = ROOT.TLatex()
    text.SetNDC(True)
    text.SetTextFont(42)
    text.SetTextSize(.04)
    text.DrawLatex(.39, .89, "H* #rightarrow ZZ #rightarrow 4l Asimov scans")
    text.SetTextSize(.031)
    text.DrawLatex(.39, .61, "127 systematics + 3 qqZZ factors profiled")
    text.DrawLatex(.39, .56, "Independent same-sample normalization; seed 42")
    text.DrawLatex(.39, .51, "Local endpoint: #mu = 0.01 (no fit at zero)")
    upper.RedrawAxis()
    lower.cd()
    limit = max(.001, 1.3*max(abs(r["delta_q"]) for rows in comparisons.values() for r in rows))
    frame2 = lower.DrawFrame(0, -limit, 1.5, limit)
    frame2.GetXaxis().SetTitle("#mu_{off-shell}")
    frame2.GetXaxis().SetTitleSize(.115)
    frame2.GetXaxis().SetLabelSize(.09)
    frame2.GetYaxis().SetTitle("#Delta t_{#mu}")
    frame2.GetYaxis().SetTitleSize(.11)
    frame2.GetYaxis().SetTitleOffset(.5)
    frame2.GetYaxis().SetLabelSize(.085)
    frame2.GetYaxis().SetNdivisions(505)
    band = ROOT.TBox(.05, -.0005, 1.5, .0005)
    band.SetFillColor(ROOT.kGray)
    band.Draw()
    zero = ROOT.TLine(0, 0, 1.5, 0)
    zero.SetLineStyle(2)
    zero.Draw()
    for size, rows in comparisons.items():
        residual = graph(f"residual_{size}", [(r["mu"], r["delta_q"]) for r in rows], colors[size], 24 if size == 4096 else 25)
        residual.Draw("LP SAME")
    text.SetTextSize(.067)
    text.DrawLatex(.17, .86, "Local - HEPData at shared #mu #geq 0.05; no shifts or rescaling")
    lower.RedrawAxis()
    canvas.cd()
    for extension in ("png", "pdf"):
        canvas.SaveAs(str(output/f"{stem}.{extension}"))
    root_file = ROOT.TFile.Open(str(output/f"{stem}.root"), "RECREATE")
    for obj in objects + [canvas]:
        obj.Write()
    root_file.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits", type=Path, required=True)
    parser.add_argument("--sample-subdirectory", default="",
                        help="Optional fit directory within each n_SIZE directory")
    parser.add_argument("--normalization-mode", choices=("nominal", "anchors"), default="nominal")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--crosscheck-table", type=Path, required=True)
    parser.add_argument("--reference-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    table = json.loads(args.reference.read_text())
    if table["doi"] != "10.17182/hepdata.156805.v1/t2":
        raise ValueError("Expected HEPData 156805 v1, Table 2")
    reference = series_points(table, "Exp NSBI")
    if series_points(json.loads(args.crosscheck_table.read_text()), "Exp NSBI") != reference:
        raise ValueError("Expected series differs between Tables 1 and 2")
    if csv_series_points(args.reference_csv, "Exp NSBI") != reference:
        raise ValueError("Expected series differs between JSON and CSV")
    scans, bestfits, comparisons, reports = {}, {}, {}, {}
    inputs = {str(p.resolve()): digest(p) for p in
              (args.reference, args.crosscheck_table, args.reference_csv)}
    for size in (4096, 8192):
        directory = args.fits/f"n_{size}"/args.sample_subdirectory
        best_path, scan_path = directory/"bestfit.json", directory/"scan_points.json"
        best = json.loads(best_path.read_text())
        validation_path = directory/"likelihood_validation.json"
        validation = json.loads(validation_path.read_text())
        if not validation["passed"] or len(validation.get("scan_checks", [])) != 31:
            raise ValueError("Require independent likelihood validation of all 31 profile points")
        if not best["complete"] or best["fit_status"] or best["covariance_quality"] != 3:
            raise ValueError("Require successful global fits with accurate covariance")
        scan = {}
        for row in json.loads(scan_path.read_text()):
            require_minimum(row)
            if abs(row["q"]-2*(row["nll"]-best["nll_min"])) > 1e-8:
                raise ValueError("Scan has a different global-minimum baseline")
            if row["q"] < -2e-4:
                raise ValueError("Profile improves the global minimum; refit before comparing")
            mu = round(row["mu"], 10)
            if mu in scan:
                raise ValueError("Duplicate scan point")
            scan[mu] = row["q"]
        expected_grid = [.01]+[i/20 for i in range(1,31)]
        if sorted(scan) != expected_grid:
            raise ValueError("Require 0.01 followed by 0.05, 0.10, ..., 1.50")
        common = sorted(scan.keys() & reference.keys())
        rows = [{"mu": mu, "q_local": scan[mu], "q_hepdata": reference[mu],
                 "delta_q": scan[mu]-reference[mu]} for mu in common]
        scans[size], bestfits[size], comparisons[size] = scan, best, rows
        reports[size] = {"mu_hat": best["mu_hat"], "q_at_0p01": scan[.01],
                         "comparison": metrics(rows),
                         "fit_edm": best["edm"], "covariance_quality": best["covariance_quality"]}
        with (args.output/f"comparison_{size}.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        inputs.update({str(p.resolve()): digest(p) for p in (best_path, scan_path, validation_path)})
    pairwise = [{"mu": mu, "q_4096": scans[4096][mu], "q_8192": scans[8192][mu],
                 "delta_q": scans[4096][mu]-scans[8192][mu]} for mu in sorted(scans[4096])]
    with (args.output/"sample_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(pairwise[0]))
        writer.writeheader()
        writer.writerows(pairwise)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "reference_doi": table["doi"], "reference_series": "Exp NSBI",
        "reference_description": table["description"], "reference_figure": table["location"],
        "reference_grid_minimum": {"mu": min(reference, key=reference.get), "q": min(reference.values())},
        "reference_q_at_zero": reference[0],
        "evaluated_range": [.01, 1.5], "plot_range": [0, 1.5],
        "comparison_range": [.05, 1.5],
        "endpoint_note": "0.01 is a local zero proxy, not the published zero coordinate; excluded from reference residuals",
        "definition": "2 * (NLL_profile(mu) - NLL_continuous_global_min)",
        "comparison_method": "Shared coordinates only; no interpolation, rescaling, or minimum shift",
        "reference_rounding_half_step": .0005,
        "reference_crosschecks": {"tables_1_and_2_identical": True, "csv_matches_json": True},
        "results": reports, "sample_difference_4096_minus_8192": metrics(pairwise),
        "normalization_mode": args.normalization_mode,
        "model_caveat": ("Finite reference samples with nominal and up/down anchor same-sample normalization to supplied physical rates; nonlinear interpolation is not normalized at every nuisance point. NumPy inputs differ from the observed ROOT template."
                         if args.normalization_mode == "anchors" else
                         "Finite reference samples with nominal-only same-sample normalization; supplied systematics and NumPy yields retained. These inputs differ from the observed ROOT template."),
        "inputs_sha256": inputs,
    }
    write_json(args.output/"comparison_summary.json", summary)
    make_plot(scans, bestfits, reference, comparisons, args.output)
    for size in scans:
        make_plot({size: scans[size]}, {size: bestfits[size]}, reference,
                  {size: comparisons[size]}, args.output,
                  stem=f"expected_scan_comparison_{size}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
