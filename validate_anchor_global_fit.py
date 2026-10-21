#!/usr/bin/env python3
"""Validate a global-only anchor-normalized fit against explicit intensities."""
import argparse
import json
from pathlib import Path

import numpy as np
import ROOT

from root_workspace_fit import load_model, restore, sha256, write_json
from root_workspace_validation import CPP, intensity_check, validate_saved_fit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--bestfit", type=Path, required=True)
    parser.add_argument("--previous-workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.set_defaults(workspace="combined", model_config="ModelConfig_4l_offshell",
                        data="obsData_4l_offshell", weighted_asimov=True)
    args = parser.parse_args()
    construction = json.loads((args.input.parent / "construction.json").read_text())
    assert construction["systematic_normalization_mode"] == "same_sample_anchors_to_supplied_rates"
    old = json.loads((args.previous_workspace.parent / "construction.json").read_text())
    assert sha256(args.previous_workspace) == old["workspace_sha256"]
    unchanged = {}
    for name in ("event_rows.npy", "reference_weights_normalized.npy", "asimov_weights.npy", "jet_category.npy"):
        new_values = np.load(args.input.parent / name, mmap_mode="r")
        previous = np.load(args.previous_workspace.parent / name, mmap_mode="r")
        np.testing.assert_array_equal(new_values, previous)
        unchanged[name] = True
    with np.load(args.input.parent / "nominal_ratios.npz") as new, np.load(args.previous_workspace.parent / "nominal_ratios.npz") as previous:
        for name in new.files:
            np.testing.assert_array_equal(new[name], previous[name])
    for name in ("sr_component_yields", "cr_component_yields", "cr_asimov_counts", "generation_parameters"):
        assert construction[name] == old[name]
        unchanged[name] = True
    best = json.loads(args.bestfit.read_text())
    assert best["input_sha256"] == sha256(args.input) and best["complete"]
    source, w, mc, data, parameters, floating, nll = load_model(args)
    ROOT.gInterpreter.Declare(CPP)
    weights = np.load(args.input.parent / "reference_weights_normalized.npy", mmap_mode="r")
    reference = ROOT.std.vector("double")(weights.tolist())
    checks = []
    settings = [("nominal", None)]
    for path in sorted(args.bestfit.parent.glob("global_fit_attempt_[0-9]*.json")):
        attempt = json.loads(path.read_text())
        settings.append((path.stem, attempt["parameters"]))
    settings.append(("selected_minimum", best["parameters"]))
    for label, records in settings:
        if records is None:
            w.loadSnapshot("asimov_generation")
        else:
            restore(w, records)
        native = nll.getVal()
        (explicit, current, continuous, invalid), integral = intensity_check(w, data, [.01, 1.5], reference)
        difference = 0.0 if not checks else (native-checks[0]["native_nll"])-(explicit-checks[0]["explicit_nll"])
        assert current > 0 and invalid == 0 and abs(difference) < 1e-6
        checks.append({"name": label, "mu": w.var("mu").getVal(), "native_nll": native,
                       "explicit_nll": explicit, "delta_nll_disagreement": difference,
                       "minimum_actual_intensity": current, "invalid_intensities": invalid,
                       "fixed_nuisance_interval_minimum_diagnostic": continuous, **integral})
    assert abs(checks[-1]["native_nll"]-best["nll_min"]) < 1e-6
    assert best["input_sha256"] == sha256(args.input)
    artifacts = validate_saved_fit(args.bestfit, require_accurate_covariance=False)
    scores = json.loads((args.input.parent / "validation.json").read_text())["nominal_systematic_scores"]
    report = {"passed": True, "input_sha256": best["input_sha256"],
              "previous_workspace_unchanged_sha256": sha256(args.previous_workspace),
              "unchanged_inputs": unchanged, "checks": checks, "saved_artifacts": artifacts,
              "mu_hat": best["mu_hat"], "mu_hat_minus_one": best["mu_hat"]-1,
              "fit_status": best["fit_status"], "edm": best["edm"],
              "start_mu_estimates": [row["mu_hat"] for row in best["attempts"]],
              "start_nll_spread": max(row["nll"] for row in best["attempts"])-min(row["nll"] for row in best["attempts"]),
              "largest_nominal_systematic_score": max(abs(x) for x in scores.values()),
              "full_hesse_performed": False, "profile_scan_performed": False,
              "scope": "Two-start unconditional fit; no claim that finite starts certify global optimality. Covariance is the Migrad estimate."}
    write_json(args.output, report)
    print(json.dumps(report, indent=2), flush=True)
    source.Close()


if __name__ == "__main__":
    main()
