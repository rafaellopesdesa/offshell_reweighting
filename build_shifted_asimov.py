#!/usr/bin/env python3
"""Generate a fully normalized, nuisance-shifted Asimov from a preserved sample."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import time

import numpy as np
import ROOT

from build_root_asimov import import_arg
from root_asimov_inputs import SR_PROCESSES
from root_workspace_fit import sha256, write_json
from same_sample_normalizer import load_normalizer


CPP = r'''
#include <RooWorkspace.h>
#include <RooAbsData.h>
#include <RooAbsReal.h>
#include <array>
#include <string>
#include <vector>
std::vector<std::vector<double>> offshell_sample_shapes(RooWorkspace &w, RooAbsData &data) {
  const std::array<std::string,10> names={"S","SBI","B","EWB","EWSBI","EWSBI10","ttV","qqZZ_0","qqZZ_1","qqZZ_2"};
  std::array<RooAbsReal*,10> shapes;
  for (int p=0;p<10;++p) shapes[p]=w.function(("prod_process_"+names[p]).c_str());
  std::vector<std::vector<double>> values(10);
  auto variables=w.allVars();
  for (int i=0;i<data.numEntries();++i) {
    const auto *row=data.get(i);
    if(std::string(row->getCatLabel("cat"))!="SR") continue;
    variables.assignValueOnly(*row);
    for (int p=0;p<10;++p) values[p].push_back(shapes[p]->getVal());
  }
  return values;
}
'''


def sample_shapes(workspace, data):
    return np.asarray([list(v) for v in ROOT.offshell_sample_shapes(workspace, data)])


def cr_expectations(workspace):
    return np.array([sum(v.getVal() for v in workspace.pdf(f"pdf_CR{c}").coefList()) for c in range(3)])


def set_generation(workspace, values):
    for name, value in values.items():
        var = workspace.var(name)
        if not var or not var.getMin() <= value <= var.getMax():
            raise ValueError(f"Invalid generating parameter {name}: {value}")
        var.setVal(value)
        var.setConstant(False)
    mc = workspace.obj("ModelConfig_4l_offshell")
    for var in mc.GetGlobalObservables():
        var.setVal(values["alpha_" + var.GetName()[4:]])
        var.setConstant(True)


def normalized_model(template, columns, u):
    workspace = ROOT.RooWorkspace("combined", "Same-reference normalized model at every nuisance point")
    template.cat("cat").setLabel("SR")
    mask = columns["cat"] == template.cat("cat").getIndex()
    if int(mask.sum()) != len(u):
        raise ValueError("Reference sample and dataset differ")
    for process in SR_PROCESSES:
        interp = template.function(f"piecewise_interp_{process}")
        if set(interp.interpolationCodes()) != {5} or interp.nominalHist().getVal() != 1 or interp.positiveDefinite():
            raise ValueError("Unsupported interpolation convention")
        low = np.array([columns[v.GetName()][mask] for v in interp.lowList()])
        high = np.array([columns[v.GetName()][mask] for v in interp.highList()])
        measure = u * columns[f"r_{process}"][mask]
        normalizer = ROOT.RooSameSampleNormalizer(f"reference_Z_{process}", "sum(u*rhat*V(alpha))",
            interp.paramList(), ROOT.std.vector("double")(measure),
            ROOT.std.vector("double")(low.ravel()), ROOT.std.vector("double")(high.ravel()))
        # h=rhat*V/R, so h*R/Zraw=rhat*V/sum(u*rhat*V).
        # Physical yield coefficients retain their original R(alpha).
        raw = template.function(f"prod_process_{process}").clone(f"raw_process_{process}")
        numerator = ROOT.RooProduct(f"normalized_numerator_{process}", "rhat*V", ROOT.RooArgList(raw, template.function(f"mod_{process}")))
        shape = ROOT.RooRatio(f"prod_process_{process}", "Same-reference normalized process shape", numerator, normalizer)
        import_arg(workspace, shape)
    import_arg(workspace, template.pdf("simPdf"))
    previous = template.obj("ModelConfig_4l_offshell")
    mc = ROOT.RooStats.ModelConfig(previous.GetName(), workspace)
    mc.SetPdf(workspace.pdf("simPdf"))
    for getter, setter in ((previous.GetObservables, mc.SetObservables),
                          (previous.GetParametersOfInterest, mc.SetParametersOfInterest),
                          (previous.GetNuisanceParameters, mc.SetNuisanceParameters),
                          (previous.GetGlobalObservables, mc.SetGlobalObservables)):
        setter(workspace.argSet(getter().contentsString()))
    if workspace.Import(mc):
        raise RuntimeError("Cannot import ModelConfig")
    return workspace


def make_dataset(workspace, columns, original, weights, cr_counts, name="obsData_4l_offshell"):
    columns = dict(columns)
    columns["total_weight"] = np.concatenate([cr_counts, weights])
    # The preserved source uses exactly CR0,CR1,CR2 then the selected SR rows.
    labels = [original.get(i).getCatLabel("cat") for i in range(4)]
    if labels != ["CR0", "CR1", "CR2", "SR"]:
        raise ValueError("Unexpected dataset row order")
    variables = ROOT.RooArgSet()
    for var in original.get():
        variables.add(workspace.arg(var.GetName()))
    weight_var = ROOT.RooRealVar("total_weight", "Fixed Asimov event weight", 0)
    variables.add(weight_var)
    data = ROOT.RooDataSet.from_numpy(columns, variables, name=name,
        title="Data-like Asimov at mu=1 with observed-fit nuisance parameters", weight_name="total_weight")
    if data.numEntries() != len(weights)+3:
        raise RuntimeError("Dataset import discarded events")
    return data, columns


def build(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run numerical workspace construction through Slurm")
    started = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=True)
    target = args.output / "HZZ_4l_data_like_asimov_4096.root"
    if target.exists():
        raise FileExistsError(target)
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    load_normalizer()
    ROOT.gInterpreter.Declare(CPP)
    source_sha = sha256(args.source)
    original_report = json.loads((args.source.parent/"construction.json").read_text())
    if original_report["workspace_sha256"] != source_sha or original_report["reference_events"] != 4096:
        raise ValueError("Input must be the preserved, checked 4096 workspace")
    u = np.load(args.source.parent/"reference_weights_normalized.npy", allow_pickle=False)
    if np.any(u < 0) or not np.isfinite(u).all() or abs(u.sum()-1)>1e-13:
        raise ValueError("This experiment requires nonnegative normalized reference weights")
    source = ROOT.TFile.Open(str(args.source), "READ")
    template = source.Get("combined")
    original = template.data("obsData_4l_offshell")
    columns = original.to_numpy(copy=True)
    workspace = normalized_model(template, columns, u)
    observed = json.loads(args.observed_bestfit.read_text())
    parameters = {r["name"]: r["value"] for r in observed["parameters"]}
    parameters["mu"] = args.mu
    expected = {v.GetName() for v in workspace.obj("ModelConfig_4l_offshell").GetNuisanceParameters()} | {"mu", "mu_qqZZ", "mu_qqZZ_1", "mu_qqZZ_2"}
    if parameters.keys() != expected:
        raise ValueError("Observed fit does not provide exactly the 131 model parameters")
    set_generation(workspace, parameters)
    shapes = sample_shapes(workspace, original)
    integrals = shapes @ u
    np.testing.assert_allclose(integrals, 1, rtol=0, atol=3e-12)
    coefficients = np.array([workspace.function(f"nfFactor_{p}").getVal() for p in SR_PROCESSES])
    intensity = coefficients @ shapes
    if np.any(intensity <= 0) or not np.isfinite(intensity).all():
        raise ValueError("Nonpositive generated intensity; no clipping is allowed")
    weights = u * intensity
    cr_counts = cr_expectations(workspace)
    np.testing.assert_allclose(weights.sum(), coefficients.sum(), rtol=3e-12)
    data, generated_columns = make_dataset(workspace, columns, original, weights, cr_counts)
    import_arg(workspace, data, recycle=False)
    weight_var = ROOT.RooRealVar("total_weight", "Asimov event weight", 0)
    for c in range(3):
        obs = workspace.var(f"obs_CR{c}")
        cr_data = ROOT.RooDataSet.from_numpy({obs.GetName(): np.array([.5]), "total_weight": cr_counts[c:c+1]},
            ROOT.RooArgSet(obs, weight_var), name=f"data_CR{c}", title="Shifted Asimov control count", weight_name="total_weight")
        import_arg(workspace, cr_data, recycle=False)
    mc = workspace.obj("ModelConfig_4l_offshell")
    snapshot = ROOT.RooArgSet(mc.GetParametersOfInterest())
    snapshot.add(mc.GetNuisanceParameters())
    snapshot.add(mc.GetGlobalObservables())
    workspace.saveSnapshot("asimov_generation", snapshot, True)
    for name, values in workspace.data(data.GetName()).to_numpy(copy=False).items():
        np.testing.assert_array_equal(values, generated_columns[name])
    # Keep the custom class implementation with the persistent workspace too.
    embedded = bool(workspace.importClassCode("RooSameSampleNormalizer", True))
    report = {
        "data_kind": "weighted_reference_asimov", "label": "data-like Asimov",
        "normalization_implementation": "RooSameSampleNormalizer_v1",
        "normalization": "rho=rhat*V/sum_reference(u*rhat*V), recomputed at every trial alpha; physical rates unchanged",
        "reference_events": 4096, "negative_reference_weights": 0,
        "source_workspace": str(args.source.resolve()), "source_workspace_sha256": source_sha,
        "observed_bestfit": str(args.observed_bestfit.resolve()), "observed_bestfit_sha256": sha256(args.observed_bestfit),
        "generation_parameters": parameters, "observed_mu_not_injected": observed["mu_hat"],
        "qqZZ_parameterization": "cumulative: beta0, beta0*beta1, beta0*beta1*beta2",
        "global_observables": {v.GetName(): v.getVal() for v in mc.GetGlobalObservables()},
        "sr_asimov_yield": float(weights.sum()), "cr_asimov_counts": cr_counts.tolist(),
        "total_asimov_yield": float(weights.sum()+cr_counts.sum()),
        "minimum_intensity": float(intensity.min()), "minimum_asimov_weight": float(weights.min()),
        "process_integrals_at_generation": dict(zip(SR_PROCESSES, integrals.tolist())),
        "dataset_roundtrip_exact": True, "dataset_rows": data.numEntries(),
        "class_code_embedded": embedded, "root_version": ROOT.gROOT.GetVersion(),
        "slurm_job_id": os.environ["SLURM_JOB_ID"], "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_files": {name: sha256(Path(__file__).parent/name) for name in
                         ("build_shifted_asimov.py", "same_sample_normalizer.py", "RooSameSampleNormalizer.h", "RooSameSampleNormalizer.cxx")},
    }
    for name in ("reference_weights_normalized.npy", "event_rows.npy", "jet_category.npy", "nominal_ratios.npz"):
        shutil.copyfile(args.source.parent/name, args.output/name)
    np.save(args.output/"asimov_weights.npy", weights)
    temporary = target.with_suffix(".tmp.root")
    output = ROOT.TFile(str(temporary), "RECREATE")
    workspace.Write("combined")
    ROOT.TObjString(json.dumps(report)).Write("asimov_metadata")
    output.Close()
    temporary.replace(target)
    report.update(output_workspace=str(target.resolve()), workspace_sha256=sha256(target), seconds=time.monotonic()-started)
    write_json(args.output/"construction.json", report)
    if sha256(args.source) != source_sha:
        raise RuntimeError("Preserved source changed during construction")
    source.Close()
    print(json.dumps({k:report[k] for k in ("output_workspace", "sr_asimov_yield", "cr_asimov_counts", "class_code_embedded", "seconds")}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--observed-bestfit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mu", type=float, default=1.)
    build(parser.parse_args())
