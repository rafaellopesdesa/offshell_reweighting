#!/usr/bin/env python3
"""Fit and profile the supplied data workspace without changing the input file.

Run within a Slurm allocation on Unity. All outputs belong in persistent /work
storage outside this public source repository.
"""

import argparse
from array import array
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import time

import ROOT


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def parameter_record(variable):
    return {
        "name": variable.GetName(),
        "value": variable.getVal(),
        "error": variable.getError(),
        "constant": bool(variable.isConstant()),
        "min": variable.getMin() if math.isfinite(variable.getMin()) else None,
        "max": variable.getMax() if math.isfinite(variable.getMax()) else None,
    }


def restore(workspace, records):
    for record in records:
        var = workspace.var(record["name"])
        var.setVal(record["value"])
        var.setError(record["error"])
        var.setConstant(record["constant"])


def load_model(args):
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    source = ROOT.TFile.Open(str(args.input), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot read {args.input}")
    workspace = source.Get(args.workspace)
    if not workspace:
        raise RuntimeError(f"Missing workspace {args.workspace}")
    configuration = workspace.obj(args.model_config)
    data = workspace.data(args.data)
    if not configuration or not data:
        raise RuntimeError("Missing ModelConfig or observed dataset")
    pdf = configuration.GetPdf()
    parameters = pdf.getParameters(data)
    floating = [v for v in parameters if v.InheritsFrom("RooRealVar") and not v.isConstant()]
    names = {v.GetName() for v in floating}
    expected = {v.GetName() for v in configuration.GetNuisanceParameters()}
    expected.update(v.GetName() for v in configuration.GetParametersOfInterest())
    if names != expected:
        raise RuntimeError(f"Unexpected floating parameter set: {names ^ expected}")
    asimov = getattr(args, "weighted_asimov", False)
    if asimov:
        metadata = source.Get("asimov_metadata")
        if not metadata or json.loads(str(metadata.GetString()))["data_kind"] != "weighted_reference_asimov":
            raise RuntimeError("Weighted-Asimov mode requires construction metadata")
        if not workspace.getSnapshot("asimov_generation"):
            raise RuntimeError("Missing Asimov generation snapshot")
    elif data.isNonPoissonWeighted():
        raise RuntimeError("This data driver expects unit events or integer count weights")
    options = [
        ROOT.RooFit.Extended(True),
        ROOT.RooFit.Constrain(configuration.GetNuisanceParameters()),
        ROOT.RooFit.GlobalObservables(configuration.GetGlobalObservables()),
        ROOT.RooFit.GlobalObservablesSource("model"),
        ROOT.RooFit.EvalBackend("cpu"),
        ROOT.RooFit.Offset("initial"),
    ]
    conditional = configuration.GetConditionalObservables()
    external = configuration.GetExternalConstraints()
    if conditional and conditional.getSize():
        options.append(ROOT.RooFit.ConditionalObservables(conditional))
    if external and external.getSize():
        options.append(ROOT.RooFit.ExternalConstraints(external))
    nll = pdf.createNLL(data, *options)
    initial = nll.getVal()
    if not math.isfinite(initial):
        raise RuntimeError("Non-finite initial likelihood")
    print(f"Loaded {data.numEntries()} rows, count sum {data.sumEntries()}, "
          f"{len(floating)} floating parameters; initial NLL {initial:.12f}", flush=True)
    return source, workspace, configuration, data, parameters, floating, nll


def minimizer(nll, strategy, verbose=False):
    result = ROOT.RooMinimizer(nll)
    result.setMinimizerType("Minuit2")
    result.setStrategy(strategy)
    result.setEps(0.01)
    result.setErrorLevel(0.5)
    result.setMaxFunctionCalls(100000)
    result.setMaxIterations(10000)
    result.setPrintLevel(1 if verbose else -1)
    result.setPrintEvalErrors(3)
    result.setOffsetting(True)
    return result


def diagnostics(result, status, nll, elapsed):
    return {
        "migrad_status": int(status),
        "fit_status": int(result.status()),
        "covariance_quality": int(result.covQual()),
        "edm": result.edm(),
        "nll": nll.getVal(),
        "invalid_nll_evaluations": int(result.numInvalidNLL()),
        "seconds": elapsed,
    }


def require_minimum(record):
    if record["migrad_status"] != 0 or not math.isfinite(record["nll"]) or record["edm"] > 1e-4:
        raise RuntimeError(f"Unconverged minimization: {record}")


def fit(args):
    source, w, mc, data, parameters, floating, nll = load_model(args)
    initial_records = [parameter_record(v) for v in floating]
    initial_nll = nll.getVal()
    attempts = []
    best_records, best_value = None, math.inf
    best_result = None
    for start_mu in args.starts:
        restore(w, initial_records)
        w.var("mu").setVal(start_mu)
        print(f"Starting global fit at mu={start_mu}", flush=True)
        minim = minimizer(nll, args.fit_strategy, verbose=True)
        t0 = time.monotonic()
        status = minim.minimize("Minuit2", "Migrad")
        result = minim.save()
        diagnostic = diagnostics(result, status, nll, time.monotonic() - t0)
        diagnostic["start_mu"] = start_mu
        diagnostic["mu_hat"] = w.var("mu").getVal()
        diagnostic["strategy"] = args.fit_strategy
        attempts.append(diagnostic)
        write_json(args.output / "global_fit_attempts.json", attempts)
        write_json(args.output / f"global_fit_attempt_{len(attempts)-1}.json", {
            "input_sha256": sha256(args.input), "diagnostics": diagnostic,
            "parameters": [parameter_record(v) for v in floating],
        })
        print("GLOBAL FIT", json.dumps(diagnostic), flush=True)
        require_minimum(diagnostic)
        if diagnostic["nll"] < best_value:
            best_value = diagnostic["nll"]
            best_records = [parameter_record(v) for v in floating]
            best_result = result.Clone("bestfit_result")
        del minim
    restore(w, best_records)
    # Publish a checkpoint before the potentially expensive full Hesse matrix.
    summary = {
        "input": str(args.input.resolve()), "input_sha256": sha256(args.input),
        "root_version": ROOT.gROOT.GetVersion(), "workspace": args.workspace,
        "model_config": args.model_config, "dataset": args.data,
        "data_kind": "weighted_reference_asimov" if args.weighted_asimov else "observed",
        "weight_treatment": "Native signed weighted extended NLL; no SumW2 correction",
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "mu_hat": w.var("mu").getVal(), "nll_min": nll.getVal(),
        "initial_nll": initial_nll, "attempts": attempts,
        "parameters": best_records, "complete": False,
    }
    write_json(args.output / "bestfit_checkpoint.json", summary)
    if args.hesse:
        minim = minimizer(nll, args.fit_strategy, verbose=True)
        status = minim.minimize("Minuit2", "Migrad")
        require_minimum(diagnostics(minim.save(), status, nll, 0.0))
        print("Starting full Hesse covariance", flush=True)
        minim.setStrategy(1)
        t0 = time.monotonic()
        hesse_status = minim.hesse()
        best_result = minim.save("bestfit_result", "Workspace global fit")
        summary["hesse_status"] = int(hesse_status)
        summary["hesse_seconds"] = time.monotonic() - t0
        if hesse_status != 0 or best_result.covQual() != 3:
            write_json(args.output / "failed_hesse.json", summary)
            raise RuntimeError("Hesse did not produce an accurate positive-definite covariance")
        require_minimum(diagnostics(best_result, status, nll, 0.0))
    summary.update({
        "mu_hat": w.var("mu").getVal(), "mu_hesse_error": w.var("mu").getError(),
        "nll_min": nll.getVal(), "fit_status": int(best_result.status()),
        "edm": best_result.edm(), "covariance_quality": int(best_result.covQual()),
        "parameters": [parameter_record(v) for v in floating], "complete": True,
        "alpha_hat": {v.GetName(): v.getVal() for v in mc.GetNuisanceParameters()},
        "qqzz_effective_scales": {name: w.function("f_qqZZ_" + name).getVal()
                                  for name in ("0", "1", "2")},
    })
    write_json(args.output / "bestfit.json", summary)
    with (args.output / "bestfit_parameters.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary["parameters"][0]))
        writer.writeheader()
        writer.writerows(summary["parameters"])
    output = ROOT.TFile.Open(str(args.output / "bestfit.root"), "RECREATE")
    best_result.Write("fit_result")
    parameters.snapshot().Write("bestfit_parameters")
    best_result.covarianceMatrix().Write("covariance")
    best_result.correlationMatrix().Write("correlation")
    output.Close()
    w.saveSnapshot("bestfit", parameters, True)
    w.writeToFile(str(args.output / "postfit_workspace.root"))
    print(f"BEST FIT mu={summary['mu_hat']:.10f}, NLL={summary['nll_min']:.12f}, "
          f"EDM={summary['edm']:.3g}, covQual={summary['covariance_quality']}", flush=True)
    source.Close()


def scan_grid(args):
    if args.grid:
        return sorted(args.grid)
    lower, upper = args.range
    return [lower + (upper-lower)*i/(args.points-1) for i in range(args.points)]


def valid_profile_start(workspace, floating, nll, nominal):
    """Move an invalid warm start toward nominal parameters at the fixed mu.

    This changes initialization only. No intensity, parameter bound, or
    likelihood term is clipped or altered; all nuisance parameters are then
    profiled normally. The nominal reference point must itself be physical.
    """
    current = {v.GetName(): v.getVal() for v in floating if v.GetName() != "mu"}
    for fraction in (1.0, 0.5, 0.25, 0.125, 0.0):
        if fraction != 1.0:
            for name, value in current.items():
                workspace.var(name).setVal(nominal[name] + fraction*(value-nominal[name]))
        if math.isfinite(nll.getVal()):
            return fraction
    raise RuntimeError("No finite profile start, including the nominal nuisance point")


def scan(args):
    source, w, mc, data, parameters, floating, nll = load_model(args)
    nominal = {v.GetName(): v.getVal() for v in floating}
    nominal_errors = {v.GetName(): v.getError() for v in floating}
    best = json.loads(args.bestfit.read_text())
    if best["input_sha256"] != sha256(args.input):
        raise RuntimeError("Best-fit checkpoint is for a different input file")
    restore(w, best["parameters"])
    if abs(nll.getVal() - best["nll_min"]) > 1e-6:
        raise RuntimeError("Restored best-fit likelihood does not match saved value")
    grid = scan_grid(args)
    lower, upper = min(grid), max(grid)
    scan_path = args.output / "scan_points.json"
    records = json.loads(scan_path.read_text()) if args.resume and scan_path.exists() else []
    done = {round(row["mu"], 10) for row in records}
    branches = [sorted((x for x in grid if x < best["mu_hat"]), reverse=True),
                sorted(x for x in grid if x >= best["mu_hat"])]
    for branch in branches:
        restore(w, best["parameters"])
        w.var("mu").setConstant(True)
        minim = minimizer(nll, 0)
        for value in branch:
            if round(value, 10) in done:
                continue
            w.var("mu").setVal(value)
            start_fraction = valid_profile_start(w, floating, nll, nominal)
            if start_fraction != 1.0:
                print(f"mu={value}: valid start at {start_fraction:g} of warm-start displacement from nominal", flush=True)
                del minim
                for var in floating:
                    if var.GetName() != "mu":
                        var.setError(nominal_errors[var.GetName()])
                minim = minimizer(nll, 0)
            t0 = time.monotonic()
            status = minim.minimize("Minuit2", "Migrad")
            result = minim.save()
            row = diagnostics(result, status, nll, time.monotonic() - t0)
            if status != 0 or row["edm"] > 1e-4:
                minim.setStrategy(1)
                status = minim.minimize("Minuit2", "Migrad")
                result = minim.save()
                row = diagnostics(result, status, nll, time.monotonic() - t0)
                minim.setStrategy(0)
            row.update({"mu": value, "q": 2.0 * (row["nll"] - best["nll_min"]),
                        "warm_start_fraction": start_fraction,
                        "parameters": {v.GetName(): v.getVal() for v in floating}})
            records.append(row)
            write_json(scan_path, sorted(records, key=lambda r: r["mu"]))
            print("SCAN", json.dumps({k: v for k, v in row.items() if k != "parameters"}), flush=True)
            require_minimum(row)
            if row["q"] < -2e-4:
                raise RuntimeError("Profile found a lower minimum than the global fit")
        del minim
    with (args.output / "likelihood_scan.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["mu", "q", "nll", "migrad_status", "edm", "seconds"])
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames}
                         for row in sorted(records, key=lambda r: r["mu"]))
    write_json(args.output / "scan_summary.json", {
        "range": [lower, upper], "points": len(records), "mu_hat": best["mu_hat"],
        "baseline_nll": best["nll_min"], "definition": "2 * (NLL_profile(mu) - NLL_global_min)",
        "profiled": 130, "all_converged": all(r["migrad_status"] == 0 and r["edm"] <= 1e-4 for r in records),
        "max_edm": max(r["edm"] for r in records),
    })
    source.Close()


def plot(args):
    ROOT.gROOT.SetBatch(True)
    best = json.loads(args.bestfit.read_text())
    points = json.loads((args.output / "scan_points.json").read_text())
    for row in points:
        require_minimum(row)
    xy = sorted([(r["mu"], r["q"]) for r in points] + [(best["mu_hat"], 0.0)])
    graph = ROOT.TGraph(len(xy), array("d", [x for x, y in xy]), array("d", [y for x, y in xy]))
    graph.SetName("profile_likelihood")
    graph.SetTitle(";Off-shell signal strength #mu;-2 ln #Lambda(#mu)")
    graph.SetLineColor(ROOT.TColor.GetColor("#2166ac"))
    graph.SetLineWidth(3)
    graph.SetMarkerColor(graph.GetLineColor())
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(0.45)
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas("likelihood_scan", "Workspace profile likelihood", 1000, 750)
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.09)
    graph.Draw("ALP")
    graph.GetXaxis().SetLimits(*args.range)
    graph.GetXaxis().SetTitleSize(0.045)
    graph.GetYaxis().SetTitleSize(0.045)
    graph.SetMinimum(0.0)
    graph.SetMaximum(max(1.6, max(y for x, y in xy) * 1.5))
    reference = ROOT.TLine(args.range[0], 1.0, args.range[1], 1.0)
    reference.SetLineStyle(2)
    reference.SetLineColor(ROOT.kGray + 2)
    reference.Draw()
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.034)
    kind = "Asimov" if args.weighted_asimov else "Observed"
    text.DrawLatex(0.16, 0.84, kind + " H* #rightarrow ZZ #rightarrow 4l workspace")
    text.DrawLatex(0.16, 0.79, "Signal region + three control regions")
    text.DrawLatex(0.16, 0.74, "127 systematic + 3 qqZZ parameters profiled")
    text.DrawLatex(0.16, 0.68, f"#hat{{#mu}} = {best['mu_hat']:.5f}")
    canvas.RedrawAxis()
    for extension in ("png", "pdf"):
        canvas.SaveAs(str(args.output / ("likelihood_scan." + extension)))
    output = ROOT.TFile.Open(str(args.output / "likelihood_scan.root"), "RECREATE")
    graph.Write()
    canvas.Write()
    output.Close()


def finalize(args):
    """Normalize a completed scan to the final, possibly refined global fit."""
    best = json.loads(args.bestfit.read_text())
    if not best["complete"] or best["input_sha256"] != sha256(args.input):
        raise RuntimeError("Require a complete best fit for this exact input file")
    paths = args.scan_parts or [args.output / "scan_points.json"]
    by_mu, overlaps = {}, []
    for path in paths:
        for row in json.loads(path.read_text()):
            require_minimum(row)
            key = round(row["mu"], 10)
            if key in by_mu:
                difference = 2 * abs(row["nll"] - by_mu[key]["nll"])
                overlaps.append({"mu": row["mu"], "absolute_q_difference": difference})
                if difference > 1e-4:
                    raise RuntimeError(f"Scan chunks disagree at mu={row['mu']}: delta q={difference}")
                if row["nll"] >= by_mu[key]["nll"]:
                    continue
            by_mu[key] = row
    points = list(by_mu.values())
    lower, upper = args.range
    expected = scan_grid(args)
    actual = sorted(row["mu"] for row in points)
    if len(actual) != len(expected) or any(abs(a-b) > 1e-9 for a, b in zip(actual, expected)):
        raise RuntimeError("The scan does not cover the complete requested grid")
    for row in points:
        require_minimum(row)
        row["q"] = 2.0 * (row["nll"] - best["nll_min"])
        if row["q"] < -2e-4:
            raise RuntimeError("A profile point is below the final global minimum")
    points.sort(key=lambda row: row["mu"])
    write_json(args.output / "scan_points.json", points)
    with (args.output / "likelihood_scan.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["mu", "q", "nll", "migrad_status", "edm", "seconds"])
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in points)
    summary_path = args.output / "scan_summary.json"
    previous = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    write_json(summary_path, {
        "range": [lower, upper], "points": len(points), "mu_hat": best["mu_hat"],
        "evaluated_range": [min(expected), max(expected)], "grid": expected,
        "input_sha256": best["input_sha256"], "bestfit": str(args.bestfit.resolve()),
        "baseline_nll": best["nll_min"],
        "initial_scan_baseline_nll": previous.get("initial_scan_baseline_nll", previous.get("baseline_nll")),
        "definition": "2 * (NLL_profile(mu) - NLL_global_min)",
        "profiled": len(best["parameters"]) - 1, "all_converged": True,
        "max_edm": max(row["edm"] for row in points),
        "q_min_grid": min(row["q"] for row in points),
        "q_max_grid": max(row["q"] for row in points),
        "scan_parts": [str(path.resolve()) for path in paths],
        "overlap_checks": overlaps,
    })
    plot(args)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["fit", "scan", "plot", "finalize"])
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workspace", default="combined")
    parser.add_argument("--model-config", default="ModelConfig_4l_offshell")
    parser.add_argument("--data", default="obsData_4l_offshell")
    parser.add_argument("--starts", type=float, nargs="+", default=[1.0, 1.4])
    parser.add_argument("--hesse", action="store_true")
    parser.add_argument("--fit-strategy", type=int, choices=(0, 1, 2), default=1,
                        help="Migrad strategy for global starts/refinement; explicit Hesse always uses 1")
    parser.add_argument("--weighted-asimov", action="store_true",
                        help="Accept signed Asimov quadrature weights with construction metadata")
    parser.add_argument("--bestfit", type=Path)
    parser.add_argument("--range", type=float, nargs=2, default=[0.5, 1.5])
    parser.add_argument("--points", type=int, default=41)
    parser.add_argument("--grid", type=float, nargs="+",
                        help="Explicit evaluated mu coordinates; --range sets the plot axis")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--scan-parts", type=Path, nargs="+",
                        help="Optional scan_points.json chunks to combine during finalization")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    args.bestfit = args.bestfit or args.output / "bestfit.json"
    if args.mode in ("fit", "scan") and not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run this fit/scan inside a Slurm allocation")
    if args.points < 2 or args.range[0] >= args.range[1]:
        parser.error("Require at least two points and an increasing scan range")
    if args.mode in ("scan", "finalize"):
        grid = scan_grid(args)
        if len(set(grid)) != len(grid) or any(not math.isfinite(x) or x <= 0 for x in grid):
            parser.error("Scan coordinates must be unique, finite, and positive; use 0.01 as the zero proxy")
    {"fit": fit, "scan": scan, "plot": plot, "finalize": finalize}[args.mode](args)


if __name__ == "__main__":
    main()
