#!/usr/bin/env python3
"""Cross-check the native RooFit NLL against its extended-intensity expression."""

import argparse
import json
import math
from pathlib import Path

import ROOT

from root_workspace_fit import load_model, restore, write_json


CPP = r'''
#include <RooWorkspace.h>
#include <RooAbsData.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// The coefficients and per-event functions are taken from the workspace.
// This assembles the extended likelihood explicitly, independently of
// RooAddPdf/RooSimultaneous normalization and constraint extraction.
std::vector<double> offshell_check_intensity(RooWorkspace &w, RooAbsData &data,
                                           double muLow, double muHigh) {
    const std::array<std::string,10> processes = {
        "S", "SBI", "B", "EWB", "EWSBI", "EWSBI10", "qqZZ_0", "qqZZ_1", "qqZZ_2", "ttV"};
    std::array<RooAbsReal*,10> shapes, coefficients;
    for (int p=0; p<10; ++p) {
        shapes[p] = w.function(("prod_process_"+processes[p]).c_str());
        coefficients[p] = w.function(("nfFactor_"+processes[p]).c_str());
    }
    auto crYield = [&](int channel) {
        auto *pdf = static_cast<RooAddPdf*>(w.pdf(("pdf_CR"+std::to_string(channel)).c_str()));
        double total = 0.;
        for (auto *arg : pdf->coefList()) total += static_cast<RooAbsReal*>(arg)->getVal();
        return total;
    };
    const double originalMu = w.var("mu")->getVal();
    std::array<double,10> current;
    std::array<double,3> currentCR;
    double srYield = 0.;
    for (int p=0;p<10;++p) {current[p]=coefficients[p]->getVal(); srYield+=current[p];}
    for (int c=0;c<3;++c) currentCR[c]=crYield(c);
    // At fixed nuisance parameters every coefficient is quadratic in sqrt(mu).
    std::array<std::array<double,10>,3> basis;
    std::array<std::array<double,3>,3> crBasis;
    for (int t=0;t<3;++t) {
        w.var("mu")->setVal(t*t);
        for (int p=0;p<10;++p) basis[t][p]=coefficients[p]->getVal();
        for (int c=0;c<3;++c) crBasis[t][c]=crYield(c);
    }
    w.var("mu")->setVal(originalMu);
    const double tLow=std::sqrt(muLow), tHigh=std::sqrt(muHigh);
    auto minimum = [&](std::array<double,3> values) {
        const double a=(values[2]-2*values[1]+values[0])/2;
        const double b=values[1]-values[0]-a;
        auto value = [&](double t) {return (a*t+b)*t+values[0];};
        double result=std::min(value(tLow),value(tHigh));
        if (a>0) {const double t=-b/(2*a); if(t>tLow && t<tHigh) result=std::min(result,value(t));}
        return result;
    };
    double continuousMin=std::numeric_limits<double>::infinity();
    double currentMin=std::numeric_limits<double>::infinity();
    double nll=srYield, bad=0.;
    for(int c=0;c<3;++c) {
        if(currentCR[c]<=0) ++bad;
        nll+=currentCR[c];
        continuousMin=std::min(continuousMin,minimum({crBasis[0][c],crBasis[1][c],crBasis[2][c]}));
    }
    RooArgSet variables=w.allVars();
    for(int i=0;i<data.numEntries();++i) {
        const RooArgSet *row=data.get(i);
        const std::string channel=row->getCatLabel("cat");
        const double weight=data.weight();
        if(channel!="SR") {
            const int c=channel.back()-'0';
            nll-=weight*std::log(currentCR.at(c));
            continue;
        }
        variables.assignValueOnly(*row);
        double intensity=0.;
        std::array<double,3> values={0.,0.,0.};
        for(int p=0;p<10;++p) {
            const double shape=shapes[p]->getVal();
            intensity+=current[p]*shape;
            for(int t=0;t<3;++t) values[t]+=basis[t][p]*shape;
        }
        currentMin=std::min(currentMin,intensity);
        continuousMin=std::min(continuousMin,minimum(values));
        if(intensity<=0 || !std::isfinite(intensity)) ++bad;
        nll-=weight*std::log(intensity);
    }
    for(auto *arg : variables) {
        const std::string name=arg->GetName();
        if(name.rfind("alpha_",0)!=0) continue;
        auto *alpha=static_cast<RooRealVar*>(arg);
        auto *nominal=w.var(("nom_"+name.substr(6)).c_str());
        const double pull=alpha->getVal()-nominal->getVal();
        nll+=.5*pull*pull;
    }
    return {nll,currentMin,continuousMin,bad};
}
'''


def validate_saved_fit(bestfit_path):
    best = json.loads(bestfit_path.read_text())
    directory = bestfit_path.parent
    output = ROOT.TFile.Open(str(directory / "bestfit.root"), "READ")
    if not output or output.IsZombie():
        raise RuntimeError("Cannot open saved fit result")
    result = output.Get("fit_result")
    if not result or result.status() != 0 or result.covQual() != 3:
        raise RuntimeError("Saved RooFitResult has invalid convergence diagnostics")
    saved = {v.GetName(): v.getVal() for v in result.floatParsFinal()}
    expected = {row["name"]: row["value"] for row in best["parameters"]}
    if saved.keys() != expected.keys() or any(abs(saved[name]-value) > 1e-12 for name, value in expected.items()):
        raise RuntimeError("ROOT and JSON fit parameters disagree")
    covariance = result.covarianceMatrix()
    decomposition = ROOT.TMatrixDSymEigen(covariance)
    smallest_eigenvalue = min(decomposition.GetEigenValues()[i] for i in range(covariance.GetNrows()))
    if smallest_eigenvalue <= 0:
        raise RuntimeError("Covariance matrix is not positive definite")
    postfit_file = ROOT.TFile.Open(str(directory / "postfit_workspace.root"), "READ")
    workspace = postfit_file.Get(best["workspace"])
    if not workspace.getSnapshot("bestfit"):
        raise RuntimeError("Missing saved bestfit workspace snapshot")
    workspace.var("mu").setVal(1.0)
    workspace.loadSnapshot("bestfit")
    if any(abs(workspace.var(name).getVal()-value) > 1e-12 for name, value in expected.items()):
        raise RuntimeError("Workspace snapshot does not restore the best-fit parameters")
    report = {"parameter_count": len(expected), "root_json_agreement": True,
              "snapshot_restores_all_parameters": True,
              "covariance_minimum_eigenvalue": smallest_eigenvalue}
    postfit_file.Close()
    output.Close()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bestfit", type=Path)
    parser.add_argument("--workspace", default="combined")
    parser.add_argument("--model-config", default="ModelConfig_4l_offshell")
    parser.add_argument("--data", default="obsData_4l_offshell")
    parser.add_argument("--weighted-asimov", action="store_true")
    parser.add_argument("--range", type=float, nargs=2, default=[0.5, 1.5])
    parser.add_argument("--scan-points", type=Path,
                        help="Also validate every saved profile point independently")
    args = parser.parse_args()
    source, w, mc, data, parameters, floating, nll = load_model(args)
    ROOT.gInterpreter.Declare(CPP)
    if args.bestfit:
        restore(w, json.loads(args.bestfit.read_text())["parameters"])
    baseline_mu = w.var("mu").getVal()
    checks = [("baseline", baseline_mu, None, 0.0),
              ("mu_low", args.range[0], None, 0.0), ("mu_high", args.range[1], None, 0.0),
              ("lumi_shift", baseline_mu, "alpha_ATLAS_alpha_LUMI", 0.2),
              ("shape_shift", baseline_mu, "alpha_ATLAS_QCD_ggZZk_Shape", -0.2),
              ("qqzz_shift", baseline_mu, "mu_qqZZ_1", 0.03)]
    records = []
    for name, mu, parameter, shift in checks:
        w.var("mu").setVal(mu)
        if parameter:
            original = w.var(parameter).getVal()
            w.var(parameter).setVal(original + shift)
        native = nll.getVal()
        explicit, min_current, min_interval, invalid = ROOT.offshell_check_intensity(w, data, *args.range)
        row = {"name": name, "mu": mu, "native_nll": native,
               "explicit_nll": explicit, "minimum_current_intensity": min_current,
               "minimum_continuous_interval_intensity": min_interval,
               "invalid_intensities": invalid}
        if records:
            row["delta_nll_disagreement"] = (native-records[0]["native_nll"])-(explicit-records[0]["explicit_nll"])
        else:
            row["delta_nll_disagreement"] = 0.0
        records.append(row)
        print(json.dumps(row), flush=True)
        if parameter:
            w.var(parameter).setVal(original)
    valid = all(r["invalid_intensities"] == 0 and
                r["minimum_continuous_interval_intensity"] > 0 and
                math.isfinite(r["delta_nll_disagreement"]) and
                abs(r["delta_nll_disagreement"]) < 1e-6 for r in records)
    report = {"passed": valid, "checks": records, "mu_interval": args.range}
    if args.bestfit:
        report["saved_artifacts"] = validate_saved_fit(args.bestfit)
    if args.scan_points:
        baseline = records[0]
        scan_checks = []
        for point in json.loads(args.scan_points.read_text()):
            for name, value in point["parameters"].items():
                w.var(name).setVal(value)
            if abs(w.var("mu").getVal() - point["mu"]) > 1e-12:
                raise RuntimeError("Saved scan coordinate and parameter disagree")
            native = nll.getVal()
            explicit, current, interval, invalid = ROOT.offshell_check_intensity(w, data, *args.range)
            residual = (native-baseline["native_nll"])-(explicit-baseline["explicit_nll"])
            scan_checks.append({"mu": point["mu"], "delta_nll_disagreement": residual,
                                "saved_nll_difference": native-point["nll"],
                                "minimum_current_intensity": current,
                                "minimum_continuous_interval_intensity": interval,
                                "invalid_intensities": invalid})
            valid = valid and invalid == 0 and current > 0 and math.isfinite(residual) and abs(residual) < 1e-6 and abs(native-point["nll"]) < 1e-6
        report["scan_checks"] = scan_checks
        report["passed"] = valid
    write_json(args.output, report)
    if not valid:
        raise RuntimeError("Extended-intensity cross-check failed; inspect the saved report")
    source.Close()


if __name__ == "__main__":
    main()
