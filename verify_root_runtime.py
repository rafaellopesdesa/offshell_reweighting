#!/usr/bin/env python3
"""Check the installed ROOT/RooFit runtime using a small synthetic fit."""

import argparse
import json
import math
from pathlib import Path

import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ROOT.gROOT.SetBatch(True)
    assert ROOT.gROOT.GetVersionInt() >= 63404, ROOT.gROOT.GetVersion()
    for library in ("libRooFit", "libRooStats", "libHistFactory", "libMinuit2"):
        assert ROOT.gSystem.Load(library) >= 0, library

    ROOT.RooRandom.randomGenerator().SetSeed(12345)
    x = ROOT.RooRealVar("x", "x", -10.0, 10.0)
    mean = ROOT.RooRealVar("mean", "mean", 0.4, -3.0, 3.0)
    sigma = ROOT.RooRealVar("sigma", "sigma", 1.2, 0.1, 4.0)
    model = ROOT.RooGaussian("model", "model", x, mean, sigma)
    observables = ROOT.RooArgSet(x)
    data = model.generate(observables, 1000)
    data.SetName("data")
    mean.setVal(0.0)
    sigma.setVal(2.0)
    result = model.fitTo(
        data,
        ROOT.RooFit.Save(True),
        ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
        ROOT.RooFit.PrintLevel(-1),
    )
    assert result.status() == 0, result.status()
    assert result.covQual() == 3, result.covQual()
    assert abs(mean.getVal() - 0.4) < 0.15, mean.getVal()
    assert abs(sigma.getVal() - 1.2) < 0.15, sigma.getVal()

    workspace = ROOT.RooWorkspace("w", "ROOT installation verification")
    import_object = getattr(workspace, "import")
    assert not import_object(model)
    assert not import_object(data)
    configuration = ROOT.RooStats.ModelConfig("ModelConfig", workspace)
    configuration.SetPdf(workspace.pdf("model"))
    configuration.SetObservables(ROOT.RooArgSet(workspace.var("x")))
    configuration.SetParametersOfInterest(ROOT.RooArgSet(workspace.var("mean")))
    configuration.SetNuisanceParameters(ROOT.RooArgSet(workspace.var("sigma")))
    assert not import_object(configuration)
    assert not import_object(result, "fit_result")
    root_path = args.output_dir / "root_smoke_workspace.root"
    workspace.writeToFile(str(root_path))
    root_file = ROOT.TFile.Open(str(root_path), "READ")
    assert root_file and not root_file.IsZombie()
    restored = root_file.Get("w")
    assert restored and restored.data("data").numEntries() == 1000
    assert restored.obj("ModelConfig")
    assert restored.obj("fit_result").status() == 0
    density = restored.pdf("model").getVal(ROOT.RooArgSet(restored.var("x")))
    assert math.isfinite(density) and density > 0, density
    assert ROOT.RooStats.HistFactory.Measurement("measurement", "measurement")
    report = {
        "root_version": ROOT.gROOT.GetVersion(),
        "fit_status": result.status(),
        "covariance_quality": result.covQual(),
        "fitted_mean": mean.getVal(),
        "fitted_sigma": sigma.getVal(),
        "events": data.numEntries(),
        "workspace_round_trip": True,
        "roofit_roostats_histfactory_minuit2": True,
    }
    root_file.Close()
    (args.output_dir / "root_smoke_result.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
