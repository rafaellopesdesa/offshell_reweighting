#!/usr/bin/env python3
"""Check interpolation, finite-reference closure, auxiliary data, and stationarity."""
import argparse
import gc
import json
from pathlib import Path
import time

import numpy as np
import ROOT

from build_shifted_asimov import CPP, sample_shapes, cr_expectations, set_generation, make_dataset
from root_asimov_inputs import SR_PROCESSES
from root_workspace_fit import load_model, sha256, write_json
from root_workspace_validation import validate_saved_fit


def gradients(workspace, nll, names):
    records = {}
    for name in names:
        var = workspace.var(name)
        original = var.getVal()
        estimates = []
        for h in (2e-4, 1e-4):
            var.setVal(original+h)
            up = nll.getVal()
            var.setVal(original-h)
            down = nll.getVal()
            estimates.append((up-down)/(2*h))
        var.setVal(original)
        records[name] = {"step_2e-4": estimates[0], "step_1e-4": estimates[1],
                         "richardson": (4*estimates[1]-estimates[0])/3}
    return records


def main(args):
    started = time.monotonic()
    source, w, mc, data, parameters, floating, nll = load_model(args)
    ROOT.gInterpreter.Declare(CPP)
    construction = json.loads((args.input.parent/"construction.json").read_text())
    generation = construction["generation_parameters"]
    u = np.load(args.input.parent/"reference_weights_normalized.npy", allow_pickle=False)
    columns = data.to_numpy(copy=True)
    assert construction["workspace_sha256"] == sha256(args.input)
    assert len(generation)==131 and generation["mu"]==1
    for name, value in generation.items():
        assert w.var(name).getVal()==value and not w.var(name).isConstant()
    for v in mc.GetGlobalObservables():
        assert v.isConstant() and v.getVal()==generation["alpha_"+v.GetName()[4:]]
    observed = json.loads(Path(construction["observed_bestfit"]).read_text())
    for row in observed["parameters"]:
        if row["name"] != "mu":
            assert generation[row["name"]]==row["value"]
    np.testing.assert_array_equal(u, np.load(Path(construction["source_workspace"]).parent/"reference_weights_normalized.npy"))
    np.testing.assert_array_equal(np.load(args.input.parent/"event_rows.npy"),
        np.load(Path(construction["source_workspace"]).parent/"event_rows.npy"))
    generation_nll = nll.getVal()
    report = {"passed": False, "input_sha256": sha256(args.input), "generation_nll": generation_nll,
              "all_130_nuisances_match_observed_fit": True, "mu_generation": 1.,
              "global_observables_fixed_at_injection": True, "reference_selection_and_weights_unchanged": True}
    if args.bestfit:
        best = json.loads(args.bestfit.read_text())
        # A strategy-0 start at the exact Asimov minimum can finish with the
        # initial positive-definite diagonal covariance (quality 1). The
        # likelihood test needs a converged minimum, not Hesse error estimates.
        # Recovery and NLL checks below remain mandatory.
        report["saved_artifacts"] = validate_saved_fit(args.bestfit, require_accurate_covariance=False)
        report["recovery"] = {"mu_hat": best["mu_hat"], "mu_difference": best["mu_hat"]-1,
            "nll_minus_generation": best["nll_min"]-generation_nll,
            "maximum_parameter_difference": max(abs(r["value"]-generation[r["name"]]) for r in best["parameters"]),
            "parameter_differences": {r["name"]: r["value"]-generation[r["name"]] for r in best["parameters"]}}
        assert abs(report["recovery"]["mu_difference"])<1e-4
        assert abs(report["recovery"]["nll_minus_generation"])<1e-6
        assert report["recovery"]["maximum_parameter_difference"]<1e-3
    else:
        # Check multiple generating points, including one-at-a-time and
        # simultaneous variations. Main injection checks all 131 scores.
        nominal = {n: (1. if n in ("mu","mu_qqZZ","mu_qqZZ_1","mu_qqZZ_2") else 0.) for n in generation}
        cases = [("observed_nuisances_mu1", generation),
                 ("single_lumi_mu05", dict(nominal, mu=.5, alpha_ATLAS_alpha_LUMI=.37)),
                 ("single_shape_mu15", dict(nominal, mu=1.5, alpha_ATLAS_QCD_ggZZk_Shape=-1.3)),
                 ("observed_nuisances_mu15", dict(generation, mu=1.5))]
        records = []
        for index, (label, values) in enumerate(cases):
            set_generation(w, values)
            shapes = sample_shapes(w, data)
            integrals = shapes @ u
            np.testing.assert_allclose(integrals, 1., rtol=0, atol=3e-12)
            direct_errors = {}
            for p in SR_PROCESSES:
                z = w.function(f"reference_Z_{p}")
                direct_errors[p] = z.getVal()/z.directValue()-1
            assert max(map(abs, direct_errors.values()))<3e-12
            coefficients = np.array([w.function(f"nfFactor_{p}").getVal() for p in SR_PROCESSES])
            intensity = coefficients @ shapes
            assert np.all(intensity>0)
            weights = u*intensity
            cr = cr_expectations(w)
            np.testing.assert_allclose(weights.sum(), coefficients.sum(), rtol=3e-12)
            if index==0:
                np.testing.assert_allclose(columns["total_weight"][3:], weights, rtol=2e-13, atol=1e-13)
                np.testing.assert_allclose(columns["total_weight"][:3], cr, rtol=2e-13)
                this_nll = nll
                names = list(generation)
            else:
                generated, _ = make_dataset(w, columns, data, weights, cr, name="validation_data")
                this_nll = mc.GetPdf().createNLL(generated, ROOT.RooFit.Extended(True),
                    ROOT.RooFit.Constrain(mc.GetNuisanceParameters()), ROOT.RooFit.GlobalObservables(mc.GetGlobalObservables()),
                    ROOT.RooFit.GlobalObservablesSource("model"), ROOT.RooFit.EvalBackend("cpu"), ROOT.RooFit.Offset("initial"))
                names = ["mu", "mu_qqZZ", "mu_qqZZ_1", "mu_qqZZ_2", "alpha_ATLAS_alpha_LUMI", "alpha_ATLAS_QCD_ggZZk_Shape"]
            scores = gradients(w, this_nll, names)
            largest = max(abs(r["richardson"]) for r in scores.values())
            post_cache_error = max(abs(w.function(f"reference_Z_{p}").getVal()/
                                       w.function(f"reference_Z_{p}").directValue()-1) for p in SR_PROCESSES)
            assert post_cache_error<3e-12
            record = {"case": label, "sr_yield": float(weights.sum()), "cr_counts": cr.tolist(),
                "minimum_intensity": float(intensity.min()), "maximum_process_integral_error": float(np.max(np.abs(integrals-1))),
                "maximum_cached_vs_direct_normalizer_error": max(map(abs, direct_errors.values())),
                "normalizer_error_after_gradient_updates": post_cache_error,
                "maximum_absolute_score": largest, "scores": scores}
            records.append(record)
            report["generating_point_checks"] = records
            write_json(args.output, report)
            print(json.dumps({k:v for k,v in record.items() if k!="scores"}), flush=True)
            if largest>2e-5:
                raise RuntimeError(f"Nonstationary generated likelihood: {label}, score {largest}")
            if index:
                del this_nll, generated
                gc.collect()
        set_generation(w, generation)
    report.update(passed=True, seconds=time.monotonic()-started)
    write_json(args.output, report)
    print(json.dumps({k:v for k,v in report.items() if k not in ("generating_point_checks", "recovery")}), flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bestfit", type=Path)
    parser.add_argument("--workspace", default="combined")
    parser.add_argument("--model-config", default="ModelConfig_4l_offshell")
    parser.add_argument("--data", default="obsData_4l_offshell")
    parser.set_defaults(weighted_asimov=True)
    main(parser.parse_args())
