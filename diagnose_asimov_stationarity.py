#!/usr/bin/env python3
"""Read-only diagnostic of the generation point of a reference Asimov model.

Run in a Slurm allocation. This checks control-region counts, native NLL
derivatives, and fits with systematic nuisance parameters fixed at zero.
It never writes a workspace or changes the full-profile fit prescription.
"""

import argparse
import json
import os
from pathlib import Path
import time

import numpy as np
import ROOT

from root_workspace_fit import diagnostics, load_model, minimizer, sha256, write_json
from root_workspace_validation import CPP, intensity_check


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workspace", default="combined")
    parser.add_argument("--model-config", default="ModelConfig_4l_offshell")
    parser.add_argument("--data", default="obsData_4l_offshell")
    parser.add_argument("--step", type=float, default=1e-4)
    parser.add_argument("--top", type=int, default=5)
    parser.set_defaults(weighted_asimov=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run this diagnostic inside a Slurm allocation")
    source, w, mc, data, parameters, floating, nll = load_model(args)
    w.loadSnapshot("asimov_generation")
    nominal = {v.GetName(): v.getVal() for v in floating}
    globals_ = {v.GetName(): v.getVal() for v in mc.GetGlobalObservables()}
    assert nominal["mu"] == 1.0
    assert all(value == (0.0 if name.startswith("alpha_") else 1.0)
               for name, value in nominal.items())
    assert all(value == 0.0 for value in globals_.values())
    baseline = nll.getVal()
    counts = {}
    for i in range(data.numEntries()):
        row = data.get(i)
        channel = row.getCatLabel("cat")
        counts[channel] = counts.get(channel, 0.0) + data.weight()
    cr = []
    for c in range(3):
        channel = f"CR{c}"
        expected = sum(v.getVal() for v in w.pdf("pdf_" + channel).coefList())
        separate = w.data("data_" + channel).sumEntries()
        cr.append({"channel": channel, "combined_data_count": counts[channel],
                   "separate_data_count": separate, "model_expected_count": expected,
                   "data_minus_model": counts[channel] - expected,
                   "poisson_score_prefactor": 1.0 - counts[channel] / expected})
    report = {
        "input_sha256": sha256(args.input), "root_version": ROOT.gROOT.GetVersion(),
        "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "generation_parameters": nominal, "global_observables": globals_,
        "nominal_nll": baseline, "channel_counts": counts, "control_regions": cr,
        "derivative_definition": "d NLL / d parameter, central finite difference",
        "nominal_scores": [], "quadrature_score_checks": [], "fixed_systematics_fits": [],
        "complete": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, report)
    print("CONTROL REGIONS", json.dumps(cr), flush=True)

    def derivative(variable, step):
        value = nominal[variable.GetName()]
        variable.setVal(value + step)
        plus = nll.getVal()
        variable.setVal(value - step)
        minus = nll.getVal()
        variable.setVal(value)
        return (plus - minus) / (2 * step)

    for variable in floating:
        score = derivative(variable, args.step)
        report["nominal_scores"].append({"parameter": variable.GetName(),
                                          "step": args.step, "score": score})
    report["nominal_scores"].sort(key=lambda row: -abs(row["score"]))
    top = [row for row in report["nominal_scores"]
           if row["parameter"].startswith("alpha_")][:args.top]
    checked = top + [row for row in report["nominal_scores"]
                     if not row["parameter"].startswith("alpha_")]
    for row in checked:
        row["check_step"] = args.step * 10
        row["check_score"] = derivative(w.var(row["parameter"]), args.step * 10)
    write_json(args.output, report)
    print("NOMINAL SCORES", json.dumps(checked), flush=True)

    ROOT.gInterpreter.Declare(CPP)
    weights_path = args.input.parent / "reference_weights_normalized.npy"
    weights = np.load(weights_path, mmap_mode="r", allow_pickle=False)
    reference_weights = ROOT.std.vector("double")(weights.tolist())
    report["reference_weights_sha256"] = sha256(weights_path)
    _, report["nominal_quadrature"] = intensity_check(w, data, [0.01, 1.5], reference_weights)
    for row in top:
        variable = w.var(row["parameter"])
        step = args.step
        variable.setVal(step)
        _, plus = intensity_check(w, data, [0.01, 1.5], reference_weights)
        variable.setVal(-step)
        _, minus = intensity_check(w, data, [0.01, 1.5], reference_weights)
        variable.setVal(0.0)
        rate_derivative = (plus["sr_model_expected_yield"] - minus["sr_model_expected_yield"]) / (2 * step)
        integral_derivative = (plus["sr_reference_integral"] - minus["sr_reference_integral"]) / (2 * step)
        check = {"parameter": row["parameter"], "step": step,
                 "sr_expected_yield_derivative": rate_derivative,
                 "sr_reference_integral_derivative": integral_derivative,
                 "predicted_nll_score": rate_derivative - integral_derivative,
                 "native_nll_score": row["score"],
                 "score_disagreement": row["score"] - rate_derivative + integral_derivative}
        report["quadrature_score_checks"].append(check)
        write_json(args.output, report)
        print("QUADRATURE SCORE", json.dumps(check), flush=True)

    for variable in floating:
        variable.setVal(nominal[variable.GetName()])
        if variable.GetName().startswith("alpha_"):
            variable.setConstant(True)
    for start in (0.8, 1.2):
        for variable in floating:
            variable.setVal(nominal[variable.GetName()])
        w.var("mu").setVal(start)
        minim = minimizer(nll, 1)
        minim.setEps(1e-7)
        started = time.monotonic()
        status = minim.minimize("Minuit2", "Migrad")
        result = minim.save()
        record = diagnostics(result, status, nll, time.monotonic() - started)
        record.update({"start_mu": start, "parameters": {
            v.GetName(): v.getVal() for v in floating if not v.isConstant()},
            "delta_nll_from_nominal": nll.getVal() - baseline,
            "fixed_systematic_count": sum(v.GetName().startswith("alpha_") for v in floating)})
        report["fixed_systematics_fits"].append(record)
        write_json(args.output, report)
        print("FIXED SYSTEMATICS FIT", json.dumps(record), flush=True)
        if status != 0 or record["edm"] > 1e-8:
            raise RuntimeError("Fixed-systematics diagnostic did not converge")
        del minim
    report["complete"] = True
    write_json(args.output, report)
    source.Close()


if __name__ == "__main__":
    main()
