#!/usr/bin/env python3
"""Refine a saved optimum and evaluate its full observed Hesse covariance.

Run through Slurm. The workspace, original fit and scan baseline remain read-only.
"""
import argparse
import json
import os
from pathlib import Path
import time

import numpy as np
import ROOT

from root_workspace_fit import load_model, minimizer, restore, sha256, write_json, diagnostics, require_minimum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--bestfit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workspace", default="combined")
    parser.add_argument("--model-config", default="ModelConfig_4l_offshell")
    parser.add_argument("--data", default="obsData_4l_offshell")
    parser.add_argument("--weighted-asimov", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    best = json.loads(args.bestfit.read_text())
    assert best["input_sha256"] == sha256(args.input)
    source, w, mc, data, parameters, floating, nll = load_model(args)
    restore(w, best["parameters"])
    baseline = nll.getVal()
    assert abs(baseline-best["nll_min"]) < 1e-6
    # RooMinimizer requires Migrad before Hesse. This independent refinement
    # initializes Minuit; it does not replace the published scan baseline.
    minim = minimizer(nll, 0, verbose=True)
    minim.setEps(0.001)
    refine_started = time.monotonic()
    refine_status = minim.minimize("Minuit2", "Migrad")
    refinement = diagnostics(minim.save(), refine_status, nll, time.monotonic()-refine_started)
    require_minimum(refinement)
    write_json(args.output/"refinement.json", {"diagnostics": refinement,
        "parameters": [{"name": v.GetName(), "value": v.getVal()} for v in floating]})
    minim.setStrategy(1)
    started = time.monotonic()
    status = minim.hesse()
    elapsed = time.monotonic()-started
    result = minim.save("fit_result", "Observed Hesse after independent refinement of the saved optimum")
    rows = list(result.floatParsFinal())
    names = [v.GetName() for v in rows]
    raw = result.covarianceMatrix()
    cov = np.array([[raw[i][j] for j in range(len(rows))] for i in range(len(rows))])
    eigenvalues = np.linalg.eigvalsh(cov)
    values = {r["name"]: r["value"] for r in best["parameters"]}
    drift = max(abs(v.getVal()-values[v.GetName()]) for v in floating)
    record = {
        "input": str(args.input.resolve()), "input_sha256": best["input_sha256"],
        "bestfit": str(args.bestfit.resolve()), "bestfit_sha256": sha256(args.bestfit),
        "information": "Observed Hessian of the native extended NLL at the saved global optimum",
        "minimization_performed": True, "refinement": refinement, "hesse_status": int(status),
        "fit_status": int(result.status()), "covariance_quality": int(result.covQual()),
        "minimum_eigenvalue": float(eigenvalues[0]), "condition_number": float(eigenvalues[-1]/eigenvalues[0]),
        "edm": result.edm(), "nll": nll.getVal(), "baseline_nll": baseline,
        "maximum_parameter_drift": drift, "seconds": elapsed,
        "likelihood_evaluations": int(minim.evalCounter()),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "root_version": ROOT.gROOT.GetVersion(),
        "parameter_order": names, "covariance": cov.tolist(),
        "parameters": [{"name": v.GetName(), "value": v.getVal(), "error": v.getError()} for v in rows],
        "globals": {v.GetName(): v.getVal() for v in mc.GetGlobalObservables()},
        "passed": bool(status == 0 and result.covQual() == 3 and eigenvalues[0] > 0
                       and drift < 0.01 and -1e-5 < nll.getVal()-baseline < 1e-6),
    }
    output = ROOT.TFile.Open(str(args.output/"hesse.root"), "RECREATE")
    result.Write("fit_result")
    raw.Write("covariance")
    output.Close()
    write_json(args.output/"covariance.json", record)
    print(json.dumps({k: v for k, v in record.items() if k not in
                      ("covariance", "parameters", "parameter_order", "globals")}), flush=True)
    source.Close()
    if not record["passed"]:
        raise RuntimeError("Full observed-Hesse validation failed")


if __name__ == "__main__":
    main()
