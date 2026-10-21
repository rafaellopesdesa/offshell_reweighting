#!/usr/bin/env python3
"""Validate a serialized Asimov workspace against its staged NumPy inputs."""

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np
import ROOT

from root_asimov_inputs import PROCESSES, SR_PROCESSES, sha256, write_json


def coefficients(mu, qq=(1.0, 1.0, 1.0)):
    t, root10 = math.sqrt(mu), math.sqrt(10)
    denominator = -10 + root10
    return {"S": mu - t, "SBI": t, "B": 1 - t,
            "EWB": ((1-root10)*mu + 9*t - 10 + root10) / denominator,
            "EWSBI": (root10*mu - 10*t) / denominator,
            "EWSBI10": (-mu + t) / denominator, "ttV": 1.0,
            "qqZZ_0": qq[0], "qqZZ_1": qq[0]*qq[1], "qqZZ_2": qq[0]*qq[1]*qq[2]}


def evaluate_numpy(report, ratios, weights, staged, mu, qq=(1, 1, 1), variation=None):
    factors = coefficients(mu, qq)
    sr_rate, intensity = 0.0, np.zeros_like(weights)
    cr_rate = np.zeros(3)
    for process in SR_PROCESSES:
        base = "qqZZ" if process.startswith("qqZZ_") else process
        event_factor, rate_factor = 1.0, 1.0
        if variation:
            index, direction = variation
            event_factor = staged[f"arr_ratio_{direction}_{base}.npy"][index, :len(weights)]
            rate_factor = staged[f"list_weight_{direction}_{process}.npy"][index]
        intensity += factors[process] * report["sr_component_yields"][process] * ratios[process] * event_factor
        sr_rate += factors[process] * report["sr_component_yields"][process] * rate_factor
    for process in PROCESSES:
        scale = np.array([factors[f"qqZZ_{c}"] for c in range(3)]) if process == "qqZZ" else factors[process]
        variation_factor = 1.0
        if variation:
            index, direction = variation
            variation_factor = staged[f"arr_ratio_{direction}_{process}_CR.npy"][index]
        cr_rate += scale * np.array(report["cr_component_yields"][process]) * variation_factor
    if np.any(intensity <= 0) or np.any(cr_rate <= 0) or sr_rate <= 0:
        raise ValueError(f"Nonpositive intensity/rate at mu={mu}, variation={variation}")
    nll = sr_rate - np.dot(weights, np.log(intensity))
    nll += cr_rate.sum() - np.dot(report["cr_asimov_counts"], np.log(cr_rate))
    if variation:
        nll += 0.5  # A single Gaussian-constrained alpha is +/-1.
    return float(nll), intensity, float(sr_rate), cr_rate


def continuous_domain(report, ratios, low=0.5, high=1.5):
    lo, hi = math.sqrt(low), math.sqrt(high)
    def minimum(values):
        a = (values[2] - 2*values[1] + values[0]) / 2
        b = values[1] - values[0] - a
        c = values[0]
        result = np.minimum((a*lo+b)*lo+c, (a*hi+b)*hi+c)
        vertex = np.zeros_like(a)
        np.divide(-b, 2*a, out=vertex, where=a > 0)
        use = (a > 0) & (vertex > lo) & (vertex < hi)
        result[use] = np.minimum(result[use], (a[use]*vertex[use]+b[use])*vertex[use]+c[use])
        return result
    values, rates = [], []
    for mu in (0, 1, 4):
        factors = coefficients(mu)
        values.append(sum(factors[p] * report["sr_component_yields"][p] * ratios[p] for p in SR_PROCESSES))
        sr = sum(factors[p] * report["sr_component_yields"][p] for p in SR_PROCESSES)
        cr = [sum((1 if p == "qqZZ" else factors[p]) * report["cr_component_yields"][p][c]
                  for p in PROCESSES) for c in range(3)]
        rates.append(np.array([sr, *cr]))
    minima, rate_minima = minimum(values), minimum(rates)
    if np.any(minima <= 0) or np.any(rate_minima <= 0):
        raise ValueError("Nominal intensities or rates fail the continuous interval check")
    return {"mu_interval": [low, high], "minimum_intensity": float(minima.min()),
            "minimum_channel_yields": dict(zip(("SR", "CR0", "CR1", "CR2"), rate_minima.tolist())),
            "nonpositive_events": int(np.sum(minima <= 0)),
            "scope": "Nominal systematic and qqZZ parameters; quadratic endpoints and interior minima"}


def validate(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run numerical validation inside a Slurm allocation")
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    directory = args.input.parent
    report = json.loads((directory / "construction.json").read_text())
    if sha256(args.input) != report["workspace_sha256"]:
        raise ValueError("Workspace differs from the construction manifest")
    stage_manifest = json.loads((args.staged / "stage_manifest.json").read_text())
    staged = {name: np.load(args.staged / name, mmap_mode="r", allow_pickle=False)
              for name in stage_manifest["outputs"]}
    size = report["reference_events"]
    ratios_archive = np.load(directory / "nominal_ratios.npz", allow_pickle=False)
    ratios = {p: ratios_archive[p] for p in SR_PROCESSES}
    weights = np.load(directory / "asimov_weights.npy", mmap_mode="r")
    u = np.load(directory / "reference_weights_normalized.npy", mmap_mode="r")
    source = ROOT.TFile.Open(str(args.input), "READ")
    if not source or source.IsZombie():
        raise RuntimeError("Cannot reopen Asimov workspace")
    embedded = json.loads(str(source.Get("asimov_metadata").GetString()))
    if any(report[key] != value for key, value in embedded.items()):
        raise ValueError("ROOT and JSON construction metadata differ")
    w = source.Get("combined")
    mc, data = w.obj("ModelConfig_4l_offshell"), w.data("obsData_4l_offshell")
    if not w.getSnapshot("asimov_generation"):
        raise ValueError("Missing generation snapshot")
    w.loadSnapshot("asimov_generation")
    arrays = data.to_numpy(copy=False)
    np.testing.assert_array_equal(arrays["total_weight"][3:], weights)
    np.testing.assert_array_equal(arrays["total_weight"][:3], report["cr_asimov_counts"])
    for process in SR_PROCESSES:
        np.testing.assert_array_equal(arrays[f"r_{process}"][3:], ratios[process])
        np.testing.assert_allclose(np.dot(u, ratios[process]), 1, rtol=0, atol=2e-14)
        base = "qqZZ" if process.startswith("qqZZ_") else process
        for direction in ("down", "up"):
            original = staged[f"arr_ratio_{direction}_{base}.npy"]
            for item in report["systematic_mapping"]:
                name = f"var_{item['workspace_name'][6:]}_ratio_{process}_{direction}"
                np.testing.assert_array_equal(arrays[name][3:], original[item["array_index"], :size])
    channel_rows, channel_counts = {}, {}
    for label in ("CR0", "CR1", "CR2", "SR"):
        w.cat("cat").setLabel(label)
        mask = arrays["cat"] == w.cat("cat").getIndex()
        channel_rows[label] = int(np.sum(mask))
        channel_counts[label] = float(np.sum(arrays["total_weight"][mask]))
    assert channel_rows == {"CR0": 1, "CR1": 1, "CR2": 1, "SR": size}
    for c in range(3):
        np.testing.assert_allclose(w.data(f"data_CR{c}").sumEntries(), report["cr_asimov_counts"][c], rtol=1e-14)
    parameters = mc.GetPdf().getParameters(data)
    floating = {v.GetName() for v in parameters if v.InheritsFrom("RooRealVar") and not v.isConstant()}
    expected = {v.GetName() for v in mc.GetNuisanceParameters()} | {v.GetName() for v in mc.GetParametersOfInterest()}
    assert floating == expected and len(floating) == 131
    nll = mc.GetPdf().createNLL(data, ROOT.RooFit.Extended(True),
                              ROOT.RooFit.Constrain(mc.GetNuisanceParameters()),
                              ROOT.RooFit.GlobalObservables(mc.GetGlobalObservables()),
                              ROOT.RooFit.GlobalObservablesSource("model"),
                              ROOT.RooFit.EvalBackend("cpu"), ROOT.RooFit.Offset("initial"))
    nominal_numpy, nominal_intensity, _, _ = evaluate_numpy(report, ratios, weights, staged, 1)
    nominal_root = nll.getVal()
    np.testing.assert_allclose(weights, u * nominal_intensity, rtol=3e-14, atol=2e-14)
    tests = [("nominal", 1.0, (1, 1, 1), None),
             ("mu_low", 0.5, (1, 1, 1), None), ("mu_high", 1.5, (1, 1, 1), None),
             ("qqZZ_cumulative", 1.2, (1.03, 0.97, 1.04), None)]
    # mu=1.3 makes every signal/interference basis coefficient nonzero, testing
    # variations even for components absent at the SM point.
    tests += [(item["workspace_name"] + "_" + direction, 1.3, (1, 1, 1),
               (item["array_index"], direction))
              for item in report["systematic_mapping"] for direction in ("down", "up")]
    checks = []
    for label, mu, qq, variation in tests:
        w.var("mu").setVal(mu)
        for name, value in zip(("mu_qqZZ", "mu_qqZZ_1", "mu_qqZZ_2"), qq):
            w.var(name).setVal(value)
        if variation:
            index, direction = variation
            alpha = w.var(report["systematic_mapping"][index]["workspace_name"])
            alpha.setVal(-1 if direction == "down" else 1)
        manual, intensity, rate, cr = evaluate_numpy(report, ratios, weights, staged, mu, qq, variation)
        native = nll.getVal()
        residual = (native - nominal_root) - (manual - nominal_numpy)
        if not math.isfinite(residual) or abs(residual) > 1e-6:
            raise ValueError(f"Independent likelihood mismatch at {label}: {residual}")
        # Check rate and shape separately so a compensation cannot hide a bug.
        root_sr_rate = sum(w.function(f"nfFactor_{p}").getVal() for p in SR_PROCESSES)
        np.testing.assert_allclose(root_sr_rate, rate, rtol=1e-12)
        for c in range(3):
            root_cr = sum(w.function(f"nu_CR{c}_{p}").getVal() for p in PROCESSES)
            np.testing.assert_allclose(root_cr, cr[c], rtol=1e-12)
        checks.append({"name": label, "mu": mu, "delta_nll_root": native-nominal_root,
                       "delta_nll_numpy": manual-nominal_numpy, "difference": residual,
                       "minimum_intensity": float(intensity.min())})
        if variation:
            alpha.setVal(0)
        if len(checks) % 25 == 0:
            print(f"Validated {len(checks)}/{len(tests)} likelihood settings", flush=True)
    w.loadSnapshot("asimov_generation")
    # Exact nominal quadrature implies zero score for mu and the free qqZZ rates.
    step = 1e-4
    scores = {}
    for parameter in ("mu", "mu_qqZZ", "mu_qqZZ_1", "mu_qqZZ_2"):
        values = []
        for delta in (-step, step):
            w.var(parameter).setVal(1+delta)
            values.append(nll.getVal())
        w.var(parameter).setVal(1)
        scores[parameter] = (values[1]-values[0])/(2*step)
    if max(abs(v) for v in scores.values()) > 2e-4:
        raise ValueError(f"Nominal score closure failed: {scores}")
    result = {"passed": True, "workspace_sha256": report["workspace_sha256"],
              "dataset_rows": data.numEntries(), "floating_parameters": len(floating),
              "channel_rows": channel_rows, "channel_counts": channel_counts,
              "negative_sr_weights": int(np.sum(weights < 0)),
              "all_nominal_and_systematic_columns_match_staged_inputs": True,
              "native_likelihood_settings_checked": len(checks),
              "maximum_absolute_delta_nll_difference": max(abs(x["difference"]) for x in checks),
              "nominal_score_finite_difference": scores,
              "continuous_domain": continuous_domain(report, ratios),
              "checks": checks,
              "scope": "No fit or profile scan; validates serialization, input mapping, normalization, weighted likelihood, rates, and nominal score closure. Systematic variations are not subset-renormalized."}
    write_json(args.output, result)
    print(json.dumps({k: v for k, v in result.items() if k != "checks"}, indent=2), flush=True)
    source.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--staged", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    validate(parser.parse_args())


if __name__ == "__main__":
    main()
