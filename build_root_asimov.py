#!/usr/bin/env python3
"""Build same-sample-normalized ROOT Asimov workspaces from staged reference rows.

The supplied observed workspace defines the graph, names, parameter bounds,
constraints, and interpolation codes. NumPy inputs define the new physics
yields and systematic variations. Neither input is modified.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import os
from pathlib import Path
import time

import numpy as np
import ROOT

from root_asimov_inputs import (PROCESSES, SR_PROCESSES, CR_ROWS, YIELD_INDEX,
                               nuisance_name, sha256, write_json)
import json


ASIMOV_PROCESSES = ("SBI", "EWSBI", "qqZZ", "ttV")


def load(directory, name):
    return np.load(directory / name, mmap_mode="r", allow_pickle=False)


def prepare_sample(directory, size, output):
    manifest = json.loads((directory / "stage_manifest.json").read_text())
    if not 0 < size <= manifest["maximum_subset"]:
        raise ValueError("Subset exceeds staged sample")
    rows = np.asarray(load(directory, "event_rows.npy")[:size])
    if len(np.unique(rows)) != size:
        raise ValueError("Repeated reference rows")
    weights = np.asarray(load(directory, "w_SM.npy")[:size], dtype=float)
    total = float(weights.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError("Nonpositive signed reference normalization")
    u = weights / total
    jets = np.asarray(load(directory, "njets.npy")[:size])
    if np.any(jets < 0) or np.any(jets != np.floor(jets)):
        raise ValueError("Invalid jet multiplicity")
    jets = np.minimum(jets, 2).astype(np.int32)
    raw_yields = load(directory, "yields_array_SR.npy")
    yields = {p: float(raw_yields[YIELD_INDEX[p]]) for p in PROCESSES}
    np.testing.assert_allclose(sum(float(raw_yields[YIELD_INDEX[f"qqZZ_{c}"]]) for c in range(3)),
                               yields["qqZZ"], rtol=1e-6)
    z, corrected = {}, {}
    for process in PROCESSES:
        ratio = np.asarray(load(directory, f"r_{process}vsRef.npy")[:size], dtype=float)
        if not np.isfinite(ratio).all() or np.any(ratio < 0):
            raise ValueError(f"Invalid nominal ratio: {process}")
        z[process] = float(np.sum(u * ratio))
        if not np.isfinite(z[process]) or z[process] <= 0:
            raise ValueError(f"Invalid normalization for {process}")
        corrected[process] = ratio / z[process]
        np.testing.assert_allclose(np.sum(u * corrected[process]), 1, rtol=0, atol=2e-14)
    fractions = [float(np.sum(u * corrected["qqZZ"] * (jets == c))) for c in range(3)]
    if min(fractions) <= 0:
        raise ValueError("A qqZZ jet fraction is nonpositive")
    np.testing.assert_allclose(sum(fractions), 1, rtol=0, atol=2e-14)
    components = {p: corrected[p] for p in SR_PROCESSES if not p.startswith("qqZZ_")}
    sr_yields = {p: yields[p] for p in components}
    for c, fraction in enumerate(fractions):
        p = f"qqZZ_{c}"
        components[p] = corrected["qqZZ"] * (jets == c) / fraction
        sr_yields[p] = yields["qqZZ"] * fraction
        np.testing.assert_allclose(np.sum(u * components[p]), 1, rtol=0, atol=2e-14)
    intensity = sum(yields[p] * corrected[p] for p in ASIMOV_PROCESSES)
    if not np.all(intensity > 0):
        raise ValueError("Nonpositive nominal event intensity; no clipping is allowed")
    asimov_weights = u * intensity
    expected_sr = sum(yields[p] for p in ASIMOV_PROCESSES)
    np.testing.assert_allclose(asimov_weights.sum(), expected_sr, rtol=2e-14)
    hist = np.asarray(load(directory, "hist_combined.npy"))
    cr_yields = {p: hist[i].astype(float) for i, p in enumerate(CR_ROWS)}
    raw_cr = load(directory, "yields_array_CR.npy")
    for p in PROCESSES:
        np.testing.assert_allclose(cr_yields[p].sum(), raw_cr[YIELD_INDEX[p]], rtol=1e-6)
    np.testing.assert_allclose(cr_yields["qqZZ"], raw_cr[[3, 4, 5]], rtol=1e-6)
    expected_cr = sum(cr_yields[p] for p in ASIMOV_PROCESSES)
    for name, array in (("event_rows", rows), ("reference_weights_normalized", u),
                        ("asimov_weights", asimov_weights), ("jet_category", jets)):
        np.save(output / f"{name}.npy", array)
    np.savez(output / "nominal_ratios.npz", **components)
    report = {
        "data_kind": "weighted_reference_asimov",
        "reference_events": size, "reference_population": manifest["population"],
        "seed": manifest["seed"], "selection": manifest["selection"],
        "subset_relationship": "Nested prefixes of one seeded permutation; normalization recomputed separately",
        "reference_subtypes": ["S", "EWSBI10"], "process_label": "subType",
        "signed_weight_convention": "Signed MC quadrature weights; no clipping or absolute-value replacement",
        "reference_selection_source": "Provided preselected reference arrays; no further process selection",
        "generation_parameters": {"mu": 1.0, "mu_qqZZ": 1.0, "mu_qqZZ_1": 1.0, "mu_qqZZ_2": 1.0,
                                  "all_alpha": 0.0, "all_global_observables": 0.0},
        "selected_reference_weight_sum": total, "negative_reference_weights": int(np.sum(weights < 0)),
        "negative_asimov_weights": int(np.sum(asimov_weights < 0)),
        "minimum_asimov_weight": float(asimov_weights.min()), "maximum_asimov_weight": float(asimov_weights.max()),
        "Z": z, "qqZZ_jet_fractions": fractions, "native_inclusive_yields": yields,
        "sr_component_yields": sr_yields,
        "cr_component_yields": {p: v.tolist() for p, v in cr_yields.items()},
        "sr_asimov_yield": expected_sr, "cr_asimov_counts": expected_cr.tolist(),
        "total_asimov_yield": float(expected_sr + expected_cr.sum()),
        "normalization_checks": {p: float(np.sum(u * r)) for p, r in components.items()},
        "normalization": "u=w/sum(w); Z_p=sum(u*r_p); rhat_p=r_p/Z_p; Asimov weights=u*sum_SM(lambda_p*rhat_p)",
        "systematics": "Provided relative up/down variations and full-sample rate factors unchanged; no per-subset systematic renormalization",
    }
    return components, sr_yields, cr_yields, asimov_weights, report


def import_arg(workspace, arg, recycle=True):
    options = [ROOT.RooFit.Silence()]
    if recycle:
        options.append(ROOT.RooFit.RecycleConflictNodes())
    if workspace.Import(arg, *options):
        raise RuntimeError(f"Could not import {arg.GetName()}")


def normalize_systematic_anchors(directory, size, components, output):
    """Match every SR anchor integral to its supplied full-sample rate factor."""
    u = load(output, "reference_weights_normalized.npy")
    corrections, records = {}, {}
    maximum_error = 0.0
    for process in SR_PROCESSES:
        base = "qqZZ" if process.startswith("qqZZ_") else process
        measure = u * components[process]
        corrections[process], records[process] = {}, {}
        for direction in ("down", "up"):
            variations = load(directory, f"arr_ratio_{direction}_{base}.npy")[:, :size]
            rates = np.asarray(load(directory, f"list_weight_{direction}_{process}.npy"), dtype=float)
            means = np.asarray(variations @ measure, dtype=float)
            if not np.all(np.isfinite(means)) or np.any(means <= 0) or np.any(rates <= 0):
                raise ValueError(f"Nonpositive anchor integral or rate: {process}, {direction}")
            scale = rates / means
            corrected_means = (variations * scale[:, None]) @ measure
            np.testing.assert_allclose(corrected_means, rates, rtol=3e-13, atol=1e-14)
            maximum_error = max(maximum_error, float(np.max(np.abs(corrected_means/rates - 1))))
            corrections[process][direction] = scale
            records[process][direction] = {"original_reference_means": means.tolist(),
                                           "target_rate_factors": rates.tolist(),
                                           "multipliers": scale.tolist()}
    report = {"formula": "v_corrected=v*R/sum(u*rhat*v); R=list_weight_{direction}_{process}",
              "scope": "Each subset, SR process, qqZZ jet bin, nuisance, and direction independently",
              "anchor_count": len(SR_PROCESSES)*127*2,
              "maximum_relative_closure_error": maximum_error,
              "rate_factors_unchanged": True, "cr_inputs_unchanged": True,
              "interpolation_codes_unchanged": True, "components": records}
    write_json(output / "anchor_normalization.json", report)
    return corrections, report


def create_model(template, directory, sr_yields, cr_yields):
    workspace = ROOT.RooWorkspace("combined", template.GetTitle())
    changes = []
    # Install the new nominal-yield functions before importing the template PDF.
    # RecycleConflictNodes then connects the unchanged graph to these functions.
    for process, value in sr_yields.items():
        name = f"nominal_yield_SR_{process}"
        old = template.function(name)
        lumi = template.var("ATLAS_LUMI")
        constant = ROOT.RooConstVar(f"asimov_yield_SR_{process}", "Input native yield", value)
        replacement = ROOT.RooFormulaVar(name, old.GetTitle(), ROOT.RooArgList(constant, lumi))
        import_arg(workspace, replacement)
        changes.append({"function": name, "template_nominal_yield": old.getVal(), "input_nominal_yield": value})
    for channel in range(3):
        for process in PROCESSES:
            name = f"n_CR{channel}_{process}"
            old = template.function(name)
            servers = list(old.servers())
            if len(servers) != 3 or servers[0].GetName() != "ATLAS_LUMI":
                raise ValueError(f"Unexpected template formula: {name}")
            value = float(cr_yields[process][channel])
            constant = ROOT.RooConstVar(f"asimov_yield_CR{channel}_{process}", "Input native yield", value)
            replacement = ROOT.RooFormulaVar(name, old.GetTitle(), ROOT.RooArgList(servers[0], servers[1], constant))
            import_arg(workspace, replacement)
            changes.append({"function": name, "template_nominal_yield": servers[2].getVal(), "input_nominal_yield": value})
    import_arg(workspace, template.pdf("simPdf"))
    syst = [str(x) for x in load(directory, "list_syst.npy")]
    mapped = [nuisance_name(x) for x in syst]
    configuration = template.obj("ModelConfig_4l_offshell")
    original_names = {v.GetName() for v in configuration.GetNuisanceParameters()}
    if len(mapped) != 127 or len(set(mapped)) != 127 or set(mapped) != original_names:
        raise ValueError(f"Systematic mapping is not bijective: {set(mapped) ^ original_names}")
    mapping = [{"array_index": i, "array_name": name, "workspace_name": mapped[i]} for i, name in enumerate(syst)]
    for process in SR_PROCESSES:
        interp = workspace.function(f"mod_{process}")
        for direction, setter in (("down", interp.setLow), ("up", interp.setHigh)):
            values = load(directory, f"list_weight_{direction}_{process}.npy")
            for name, value in zip(mapped, values):
                if value <= 0:
                    raise ValueError(f"Nonpositive exponential rate modifier: {process}, {name}")
                setter(workspace.var(name), float(value))
    for channel in range(3):
        for process in PROCESSES:
            interp = workspace.function(f"interpFactor_CR{channel}_{process}")
            for direction, setter in (("down", interp.setLow), ("up", interp.setHigh)):
                values = load(directory, f"arr_ratio_{direction}_{process}_CR.npy")[:, channel]
                for name, value in zip(mapped, values):
                    if value <= 0:
                        raise ValueError(f"Nonpositive CR rate modifier: {process}, {name}")
                    setter(workspace.var(name), float(value))
    mc = ROOT.RooStats.ModelConfig(configuration.GetName(), workspace)
    mc.SetPdf(workspace.pdf("simPdf"))
    for getter, setter in ((configuration.GetObservables, mc.SetObservables),
                           (configuration.GetParametersOfInterest, mc.SetParametersOfInterest),
                           (configuration.GetNuisanceParameters, mc.SetNuisanceParameters),
                           (configuration.GetGlobalObservables, mc.SetGlobalObservables)):
        setter(workspace.argSet(getter().contentsString()))
    if workspace.Import(mc):
        raise RuntimeError("Cannot import ModelConfig")
    for name in mapped:
        workspace.var(name).setVal(0)
    for name in ("mu", "mu_qqZZ", "mu_qqZZ_1", "mu_qqZZ_2"):
        workspace.var(name).setVal(1)
    for var in workspace.obj(mc.GetName()).GetGlobalObservables():
        var.setVal(0)
        var.setConstant(True)
    return workspace, mapping, changes


def create_data(workspace, template, directory, size, components, weights, cr_counts, mapping,
                anchor_corrections=None):
    original = template.data("obsData_4l_offshell")
    category = workspace.cat("cat")
    category_codes = {}
    for label in ("CR0", "CR1", "CR2", "SR"):
        category.setLabel(label)
        category_codes[label] = category.getIndex()
    nrows = size + 3
    columns = {"cat": np.array([category_codes[f"CR{c}"] for c in range(3)] + [category_codes["SR"]] * size, dtype=np.int32),
               "total_weight": np.concatenate([cr_counts, weights])}
    for channel in range(3):
        columns[f"obs_CR{channel}"] = np.full(nrows, 0.5)
    for process in SR_PROCESSES:
        columns[f"r_{process}"] = np.concatenate([np.zeros(3), components[process]])
        base = "qqZZ" if process.startswith("qqZZ_") else process
        for direction in ("down", "up"):
            variations = load(directory, f"arr_ratio_{direction}_{base}.npy")
            for item in mapping:
                name = f"var_{item['workspace_name'][6:]}_ratio_{process}_{direction}"
                values = variations[item["array_index"], :size]
                if anchor_corrections is not None:
                    # Preserve the intended full-MC rate while closing this subset's anchor.
                    values = values * anchor_corrections[process][direction][item["array_index"]]
                columns[name] = np.concatenate([np.ones(3), values])
    expected_names = {v.GetName() for v in original.get()} | {"total_weight"}
    if set(columns) != expected_names:
        raise ValueError(f"Dataset columns differ from template: {set(columns) ^ expected_names}")
    variables = ROOT.RooArgSet()
    for variable in original.get():
        variables.add(workspace.arg(variable.GetName()))
    weight_var = ROOT.RooRealVar("total_weight", "Asimov event weight", 0)
    variables.add(weight_var)
    data = ROOT.RooDataSet.from_numpy(columns, variables, name=original.GetName(),
                                     title=f"SM Asimov: {size} reweighted reference events", weight_name="total_weight")
    if data.numEntries() != nrows:
        raise RuntimeError("Dataset import discarded events")
    import_arg(workspace, data, recycle=False)
    # Preserve the template's individual CR datasets as well as the combined data.
    for channel in range(3):
        obs = workspace.var(f"obs_CR{channel}")
        cr_data = ROOT.RooDataSet.from_numpy(
            {obs.GetName(): np.array([0.5]), "total_weight": np.array([cr_counts[channel]])},
            ROOT.RooArgSet(obs, weight_var), name=f"data_CR{channel}", title="Asimov control-region count",
            weight_name="total_weight")
        import_arg(workspace, cr_data, recycle=False)
    # Save only parameter/global-observable values in the generation snapshot.
    mc = workspace.obj("ModelConfig_4l_offshell")
    snapshot = ROOT.RooArgSet(mc.GetParametersOfInterest())
    snapshot.add(mc.GetNuisanceParameters())
    snapshot.add(mc.GetGlobalObservables())
    workspace.saveSnapshot("asimov_generation", snapshot, True)
    return columns


def check_structure(workspace, template):
    for getter in ("allPdfs", "allFunctions", "allVars", "allCats"):
        old = {x.GetName(): x.ClassName() for x in getattr(template, getter)()}
        new = {x.GetName(): x.ClassName() for x in getattr(workspace, getter)()}
        if old != new:
            raise ValueError(f"Model structure differs in {getter}: {old.keys() ^ new.keys()}")
    old_mc, mc = template.obj("ModelConfig_4l_offshell"), workspace.obj("ModelConfig_4l_offshell")
    for getter in ("GetObservables", "GetParametersOfInterest", "GetNuisanceParameters", "GetGlobalObservables"):
        assert {x.GetName() for x in getattr(old_mc, getter)()} == {x.GetName() for x in getattr(mc, getter)()}
    for var in list(mc.GetParametersOfInterest()) + list(mc.GetNuisanceParameters()) + list(mc.GetGlobalObservables()):
        original = template.var(var.GetName())
        assert (var.getMin(), var.getMax(), var.isConstant()) == (original.getMin(), original.getMax(), original.isConstant())
    for process in SR_PROCESSES:
        for prefix in ("mod_", "piecewise_interp_"):
            before, after = template.function(prefix + process), workspace.function(prefix + process)
            assert list(before.interpolationCodes()) == list(after.interpolationCodes())
        assert not workspace.function("piecewise_interp_" + process).positiveDefinite()
    for channel in range(3):
        for process in PROCESSES:
            name = f"interpFactor_CR{channel}_{process}"
            assert list(template.function(name).interpolationCodes()) == list(workspace.function(name).interpolationCodes())
    return {"graph_names_and_classes_match": True, "modelconfig_sets_match": True,
            "parameter_bounds_and_constants_match": True, "interpolation_codes_preserved": True,
            "pdf_classes": dict(Counter(x.ClassName() for x in workspace.allPdfs()))}


def build(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Build the workspace inside a Slurm allocation")
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    args.output.mkdir(parents=True, exist_ok=True)
    target = args.output / f"HZZ_4l_asimov_{args.size}.root"
    if target.exists():
        raise FileExistsError(target)
    started = time.monotonic()
    stage = json.loads((args.staged / "stage_manifest.json").read_text())
    for name, item in stage["outputs"].items():
        if sha256(args.staged / name) != item["sha256"]:
            raise ValueError(f"Staged input checksum changed: {name}")
    components, yields, cr_yields, weights, report = prepare_sample(args.staged, args.size, args.output)
    anchor_corrections = None
    report["systematic_normalization_mode"] = "provided"
    if args.normalize_systematic_anchors:
        anchor_corrections, anchor_report = normalize_systematic_anchors(
            args.staged, args.size, components, args.output)
        report["systematic_normalization_mode"] = "same_sample_anchors_to_supplied_rates"
        report["systematics"] = "Up/down event factors scaled to supplied rate factors on this subset; interpolation and full-MC rate factors unchanged"
        report["anchor_normalization_sha256"] = sha256(args.output / "anchor_normalization.json")
        report["anchor_normalization_maximum_relative_error"] = anchor_report["maximum_relative_closure_error"]
    source = ROOT.TFile.Open(str(args.template), "READ")
    if not source or source.IsZombie():
        raise RuntimeError("Cannot open template")
    template = source.Get("combined")
    workspace, mapping, changes = create_model(template, args.staged, yields, cr_yields)
    columns = create_data(workspace, template, args.staged, args.size, components, weights,
                          np.array(report["cr_asimov_counts"]), mapping, anchor_corrections)
    report["structure_validation"] = check_structure(workspace, template)
    data = workspace.data("obsData_4l_offshell")
    exported = data.to_numpy(copy=False)
    for name, values in columns.items():
        np.testing.assert_array_equal(exported[name], values, err_msg=f"Dataset column mismatch: {name}")
    report.update({"created_utc": datetime.now(timezone.utc).isoformat(), "ROOT_version": ROOT.gROOT.GetVersion(),
                   "slurm_job_id": os.environ["SLURM_JOB_ID"], "workspace": "combined", "pdf": "simPdf",
                   "model_config": "ModelConfig_4l_offshell", "data": "obsData_4l_offshell",
                   "template_path": str(args.template.resolve()), "template_sha256": sha256(args.template),
                   "staging_manifest_sha256": sha256(args.staged / "stage_manifest.json"),
                   "dataset_rows": data.numEntries(), "dataset_columns_without_weight": data.get().getSize(),
                   "dataset_sum_weights": data.sumEntries(), "dataset_roundtrip_exact": True,
                   "systematic_mapping": mapping, "yield_replacements": changes})
    np.testing.assert_allclose(data.sumEntries(), report["total_asimov_yield"], rtol=2e-13)
    temporary = target.with_suffix(".tmp.root")
    output_file = ROOT.TFile(str(temporary), "RECREATE")
    workspace.Write("combined")
    ROOT.TObjString(json.dumps(report)).Write("asimov_metadata")
    output_file.Close()
    temporary.replace(target)
    report["output_workspace"] = str(target.resolve())
    report["workspace_sha256"] = sha256(target)
    report["seconds"] = time.monotonic() - started
    write_json(args.output / "construction.json", report)
    source.Close()
    print(json.dumps({k: report[k] for k in ("output_workspace", "dataset_rows", "dataset_sum_weights",
                                             "negative_asimov_weights", "Z", "qqZZ_jet_fractions", "seconds")}, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--staged", type=Path, required=True)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--normalize-systematic-anchors", action="store_true",
                        help="Normalize each SR up/down anchor to its supplied varied rate on this subset")
    build(parser.parse_args())


if __name__ == "__main__":
    main()
