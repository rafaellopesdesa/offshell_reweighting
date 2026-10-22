#!/usr/bin/env python3
"""Check exact restriction, dataset preservation, closure and fit recovery."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import ROOT

from build_shifted_asimov import CPP, sample_shapes, cr_expectations
from prune_data_like_asimov import dataset_digest
from root_asimov_inputs import SR_PROCESSES
from root_workspace_fit import load_model, sha256, write_json
from root_workspace_validation import validate_saved_fit
from validate_shifted_asimov import gradients


def main(args):
    started = time.monotonic()
    source, w, mc, data, parameters, floating, nll = load_model(args)
    construction = json.loads((args.input.parent / "construction.json").read_text())
    assert construction["workspace_sha256"] == sha256(args.input)
    generation = construction["generation_parameters"]
    active = construction["active_generation_parameters"]
    pruning = construction["pruning"]
    assert set(v.GetName() for v in floating) == set(active)
    assert mc.GetNuisanceParameters().getSize() == pruning["constrained_count"]
    assert mc.GetGlobalObservables().getSize() == pruning["constrained_count"]
    for name, value in generation.items():
        assert w.var(name).getVal() == value
        assert bool(w.var(name).isConstant()) == (name in pruning["fixed_constrained"])
    for v in mc.GetGlobalObservables():
        assert v.isConstant() and v.getVal() == generation["alpha_" + v.GetName()[4:]]
    assert {d.GetName(): dataset_digest(d) for d in w.allData()} == construction["dataset_sha256"]
    generation_nll = nll.getVal()
    report = {"passed": False, "input_sha256": sha256(args.input), "generation_nll": generation_nll,
              "floating_count": len(floating), "dataset_roundtrip_exact": True,
              "fixed_parameters_at_generation": pruning["fixed_values"], "checks": []}
    if args.bestfit:
        best = json.loads(args.bestfit.read_text())
        assert best["input_sha256"] == sha256(args.input)
        assert set(r["name"] for r in best["parameters"]) == set(active)
        report["saved_artifacts"] = validate_saved_fit(args.bestfit)
        report["recovery"] = {"mu_hat": best["mu_hat"], "mu_difference": best["mu_hat"]-1,
            "nll_minus_generation": best["nll_min"]-generation_nll,
            "maximum_parameter_difference": max(abs(r["value"]-generation[r["name"]]) for r in best["parameters"])}
        assert abs(report["recovery"]["mu_difference"]) < 1e-4
        assert abs(report["recovery"]["nll_minus_generation"]) < 1e-6
        assert report["recovery"]["maximum_parameter_difference"] < 1e-3
    else:
        original_path = Path(construction["parent_workspace"])
        assert sha256(original_path) == construction["parent_workspace_sha256"]
        original_args = argparse.Namespace(**vars(args)); original_args.input = original_path
        old_source, old, old_mc, old_data, _, _, old_nll = load_model(original_args)
        assert {d.GetName(): dataset_digest(d) for d in old.allData()} == construction["dataset_sha256"]
        ROOT.gInterpreter.Declare(CPP)
        u = np.load(args.input.parent / "reference_weights_normalized.npy", allow_pickle=False)
        np.testing.assert_array_equal(u, np.load(original_path.parent / "reference_weights_normalized.npy"))
        old_baseline = old_nll.getVal()
        cases = [("generation", dict(active))]
        for mu in (0., .01, .5, 1.5, 2.):
            values = dict(active, mu=mu)
            if mu in (.5, 1.5):
                for i, name in enumerate(pruning["kept_constrained"]):
                    values[name] += .12 * np.sin(i + mu)
                values["mu_qqZZ"] *= 1.02
            cases.append((f"restriction_mu_{mu}", values))
        for label, values in cases:
            for name, value in values.items():
                w.var(name).setVal(value); old.var(name).setVal(value)
            shapes, old_shapes = sample_shapes(w, data), sample_shapes(old, old_data)
            np.testing.assert_allclose(shapes, old_shapes, rtol=3e-12, atol=1e-12)
            integrals = shapes @ u
            np.testing.assert_allclose(integrals, 1., rtol=0, atol=3e-12)
            coefficients = np.array([w.function(f"nfFactor_{p}").getVal() for p in SR_PROCESSES])
            np.testing.assert_array_equal(coefficients, [old.function(f"nfFactor_{p}").getVal() for p in SR_PROCESSES])
            np.testing.assert_array_equal(cr_expectations(w), cr_expectations(old))
            intensity = coefficients @ shapes
            assert np.isfinite(intensity).all() and np.min(intensity) > 0
            residual = (nll.getVal()-generation_nll) - (old_nll.getVal()-old_baseline)
            assert abs(residual) < 1e-7
            direct_error = max(abs(w.function(f"reference_Z_{p}").getVal()/w.function(f"reference_Z_{p}").directValue()-1) for p in SR_PROCESSES)
            assert direct_error < 3e-12
            report["checks"].append({"case": label, "delta_nll_disagreement": residual,
                "minimum_intensity": float(np.min(intensity)),
                "maximum_process_integral_error": float(np.max(np.abs(integrals-1))),
                "maximum_normalizer_direct_error": direct_error})
            write_json(args.output, report)
        for name, value in active.items():
            w.var(name).setVal(value)
        report["scores"] = gradients(w, nll, list(active))
        report["maximum_absolute_score"] = max(abs(r["richardson"]) for r in report["scores"].values())
        assert report["maximum_absolute_score"] < 2e-5
        old_source.Close()
    report.update(passed=True, seconds=time.monotonic()-started)
    write_json(args.output, report)
    source.Close()
    print(json.dumps({k:v for k,v in report.items() if k not in ("scores", "fixed_parameters_at_generation")}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bestfit", type=Path)
    parser.add_argument("--workspace", default="combined")
    parser.add_argument("--model-config", default="ModelConfig_4l_offshell")
    parser.add_argument("--data", default="obsData_4l_offshell")
    parser.set_defaults(weighted_asimov=True)
    main(parser.parse_args())
