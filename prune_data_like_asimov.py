#!/usr/bin/env python3
"""Restrict the saved data-like model to the first N ranked constrained NPs.

Excluded nuisance coordinates are fixed at their generating values. Their
shape/rate factors remain in the graph, so this is an exact restriction of
the existing likelihood, not a new approximation to the central prediction.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

import numpy as np
import ROOT

from root_workspace_fit import sha256, write_json
from same_sample_normalizer import load_normalizer


def dataset_digest(data):
    digest = hashlib.sha256()
    # PyROOT exposes C++ std::string keys here; convert before Python sorting
    # to avoid invoking cppyy's overloaded C++ comparison machinery.
    columns = {str(name): values for name, values in data.to_numpy(copy=False).items()}
    for name, values in sorted(columns.items()):
        array = np.ascontiguousarray(values)
        digest.update(name.encode() + b"\0")
        digest.update(array.dtype.str.encode() + str(array.shape).encode())
        digest.update(array.tobytes())
    return digest.hexdigest()


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Workspace processing requires a Slurm allocation")
    started = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=True)
    target = args.output / f"HZZ_4l_data_like_asimov_4096_pruned{args.keep}.root"
    if target.exists():
        raise FileExistsError(target)
    rows = list(csv.DictReader(args.ranking.open()))
    constrained = [r for r in rows if r["has_gaussian_auxiliary"] == "True"]
    if [int(r["auxiliary_rank"]) for r in constrained] != list(range(1, 128)):
        raise ValueError("Expected the published complete 127-NP shared ranking")
    if not 1 <= args.keep <= len(constrained):
        raise ValueError("Invalid retained NP count")
    keep = [r["parameter"] for r in constrained[:args.keep]]
    remove = [r["parameter"] for r in constrained[args.keep:]]
    source_sha = sha256(args.input)
    original = json.loads((args.input.parent / "construction.json").read_text())
    if original["workspace_sha256"] != source_sha or original["label"] != "data-like Asimov":
        raise ValueError("Input must be the preserved, validated data-like workspace")
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    load_normalizer()
    source = ROOT.TFile.Open(str(args.input), "READ")
    w = source.Get("combined")
    mc = w.obj("ModelConfig_4l_offshell")
    generation = original["generation_parameters"]
    if set(v.GetName() for v in mc.GetNuisanceParameters()) != set(keep + remove):
        raise ValueError("Ranking does not match the input constrained nuisances")
    w.loadSnapshot("asimov_generation")
    digests = {d.GetName(): dataset_digest(d) for d in w.allData()}
    for name, value in generation.items():
        var = w.var(name)
        if var.getVal() != value:
            raise ValueError("Generation snapshot disagrees with construction metadata")
        var.setConstant(name in remove)
    # Preserve the template's POI grouping: mu and the three free qqZZ scales.
    # The scan fixes only mu, leaving 54 alphas plus three qqZZ scales floating.
    if set(v.GetName() for v in mc.GetParametersOfInterest()) != {"mu", "mu_qqZZ", "mu_qqZZ_1", "mu_qqZZ_2"}:
        raise ValueError("Unexpected signal/background normalization inventory")
    mc.SetNuisanceParameters(w.argSet(",".join(keep)))
    mc.SetGlobalObservables(w.argSet(",".join("nom_" + n[6:] for n in keep)))
    w.defineSet("pruned_fixed_parameters", w.argSet(",".join(remove)))
    snapshot = w.argSet(",".join(generation))
    snapshot.add(w.argSet(",".join("nom_" + n[6:] for n in keep + remove)))
    w.saveSnapshot("asimov_generation", snapshot, True)
    w.SetTitle(f"Data-like Asimov: {args.keep} constrained nuisances plus three free qqZZ scales")
    floating = {v.GetName() for v in mc.GetPdf().getParameters(w.data("obsData_4l_offshell"))
                if v.InheritsFrom("RooRealVar") and not v.isConstant()}
    expected = set(keep) | {"mu", "mu_qqZZ", "mu_qqZZ_1", "mu_qqZZ_2"}
    if floating != expected:
        raise ValueError("Unexpected floating parameters after pruning")
    report = {**original, "label": f"data-like Asimov, pruned to {args.keep} constrained NPs",
        "parent_workspace": str(args.input.resolve()), "parent_workspace_sha256": source_sha,
        "ranking_file": str(args.ranking.resolve()), "ranking_sha256": sha256(args.ranking),
        "pruning": {"kept_constrained": keep, "fixed_constrained": remove,
            "fixed_values": {name: generation[name] for name in remove},
            "retained_qqzz": ["mu_qqZZ", "mu_qqZZ_1", "mu_qqZZ_2"],
            "constrained_count": args.keep, "profiled_nuisance_count": args.keep + 3,
            "global_floating_count": args.keep + 4,
            "definition": "Exact conditional restriction: excluded NPs fixed at generating values, removed from ModelConfig nuisance/global-observable sets. Their constant shape/rate factors are retained; the dataset is unchanged."},
        "active_generation_parameters": {name: generation[name] for name in generation if name in expected},
        "global_observables": {str(v.GetName()): v.getVal() for v in mc.GetGlobalObservables()},
        "fixed_auxiliary_observables": {"nom_" + n[6:]: generation[n] for n in remove},
        "dataset_sha256": digests, "dataset_unchanged": True,
        "slurm_job_id": os.environ["SLURM_JOB_ID"], "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_files": {Path(__file__).name: sha256(Path(__file__))}}
    for name in ("workspace_sha256", "seconds", "output_workspace"):
        report.pop(name, None)
    for name in ("reference_weights_normalized.npy", "event_rows.npy", "jet_category.npy", "nominal_ratios.npz", "asimov_weights.npy"):
        shutil.copyfile(args.input.parent / name, args.output / name)
    shutil.copyfile(args.ranking, args.output / "all_np_uncertainty_ranking.csv")
    with (args.output / "kept_nuisances.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(constrained[:args.keep])
    temporary = target.with_suffix(".tmp.root")
    output = ROOT.TFile(str(temporary), "RECREATE")
    w.Write("combined")
    ROOT.TObjString(json.dumps(report)).Write("asimov_metadata")
    output.Close()
    temporary.replace(target)
    report.update(output_workspace=str(target.resolve()), workspace_sha256=sha256(target), seconds=time.monotonic()-started)
    write_json(args.output / "construction.json", report)
    if sha256(args.input) != source_sha:
        raise RuntimeError("Source workspace changed")
    source.Close()
    print(json.dumps({k: report[k] for k in ("output_workspace", "workspace_sha256", "pruning", "seconds")}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--ranking", type=Path, required=True)
    parser.add_argument("--keep", type=int, default=54)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args())
