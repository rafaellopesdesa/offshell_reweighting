#!/usr/bin/env python3
"""Verify true/data-like models agree at shared parameters on the same sample."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import ROOT

from build_shifted_asimov import CPP, sample_shapes, cr_expectations, set_generation
from root_asimov_inputs import SR_PROCESSES
from root_workspace_fit import sha256, write_json
from same_sample_normalizer import load_normalizer


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Model comparisons require a Slurm allocation")
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    load_normalizer()
    ROOT.gInterpreter.Declare(CPP)
    paths = [args.true, args.data_like]
    hashes = [sha256(p) for p in paths]
    files = [ROOT.TFile.Open(str(p), "READ") for p in paths]
    models = [f.Get("combined") for f in files]
    reports = [json.loads(str(f.Get("asimov_metadata").GetString())) for f in files]
    data = [w.data("obsData_4l_offshell") for w in models]
    assert reports[0]["label"] == "true Asimov" and reports[1]["label"] == "data-like Asimov"
    assert reports[0]["reference_events"] == reports[1]["reference_events"]
    assert reports[0]["source_workspace_sha256"] == reports[1]["source_workspace_sha256"]
    for sidecar in ("event_rows.npy", "reference_weights_normalized.npy", "jet_category.npy"):
        np.testing.assert_array_equal(np.load(paths[0].parent/sidecar), np.load(paths[1].parent/sidecar))
    cases = [("nominal", reports[0]["generation_parameters"]),
             ("observed_nuisances", reports[1]["generation_parameters"])]
    simultaneous = dict(cases[1][1], mu=1.5)
    for index, name in enumerate(sorted(simultaneous)):
        if name.startswith("alpha_"):
            simultaneous[name] += .1*np.sin(index)
    cases.append(("simultaneous_shift", simultaneous))
    checks = []
    for label, values in cases:
        for w in models:
            set_generation(w, values)
        shapes = [sample_shapes(w, d) for w, d in zip(models, data)]
        np.testing.assert_allclose(*shapes, rtol=3e-12, atol=1e-12)
        rates = [np.array([w.function(f"nfFactor_{p}").getVal() for p in SR_PROCESSES]) for w in models]
        np.testing.assert_array_equal(*rates)
        np.testing.assert_array_equal(*[cr_expectations(w) for w in models])
        checks.append(dict(case=label, maximum_absolute_shape_difference=float(np.max(np.abs(shapes[0]-shapes[1]))),
                           physical_rates_identical=True, cr_expectations_identical=True))
    # True generation must preserve the earlier central Asimov observations.
    original = ROOT.TFile.Open(reports[0]["source_workspace"], "READ")
    old = original.Get("combined").data("obsData_4l_offshell").to_numpy(copy=True)
    new = data[0].to_numpy(copy=True)
    assert set(old) == set(new)
    for name in old:
        if str(name) == "total_weight":
            np.testing.assert_allclose(new[name], old[name], rtol=3e-12, atol=1e-13)
        else:
            np.testing.assert_array_equal(new[name], old[name])
    report = dict(passed=True, reference_events=reports[0]["reference_events"],
                  input_sha256=hashes, same_reference_and_shape_model=True, checks=checks,
                  nominal_dataset_maximum_weight_change=float(np.max(np.abs(new["total_weight"]-old["total_weight"]))),
                  nominal_observables_unchanged=True)
    original.Close()
    for f in files:
        f.Close()
    assert [sha256(p) for p in paths] == hashes
    write_json(args.output, report)
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--true", type=Path, required=True)
    parser.add_argument("--data-like", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args())
