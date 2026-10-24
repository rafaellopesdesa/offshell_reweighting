#!/usr/bin/env python3
"""Primary-data score covariance at a fitted point, retaining all profiled NPs.

Weighted Asimov rows are quadrature masses of a Poisson intensity, not
experimental observations with unit variance or SumW2 errors. Signed
quadrature requires explicit opt-in and is reported as an integration estimate.
"""
import argparse
import json
import os
from pathlib import Path
import time

import numpy as np
import ROOT

from root_workspace_fit import load_model, restore, sha256, write_json


CPP = r'''
#include <RooWorkspace.h>
#include <RooAbsData.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <array>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
std::vector<double> offshell_log_intensities(RooWorkspace &w, RooAbsData &data) {
  const std::array<std::string,10> processes={"S","SBI","B","EWB","EWSBI","EWSBI10","qqZZ_0","qqZZ_1","qqZZ_2","ttV"};
  std::array<RooAbsReal*,10> shapes;
  std::array<double,10> coefficients;
  for(int p=0;p<10;++p) {
    shapes[p]=w.function(("prod_process_"+processes[p]).c_str());
    coefficients[p]=w.function(("nfFactor_"+processes[p]).c_str())->getVal();
  }
  std::array<double,3> cr={0.,0.,0.};
  for(int c=0;c<3;++c) {
    auto *pdf=static_cast<RooAddPdf*>(w.pdf(("pdf_CR"+std::to_string(c)).c_str()));
    for(auto *arg:pdf->coefList()) cr[c]+=static_cast<RooAbsReal*>(arg)->getVal();
  }
  RooArgSet variables=w.allVars();
  std::vector<double> values;
  for(int i=0;i<data.numEntries();++i) {
    const auto *row=data.get(i);
    const std::string channel=row->getCatLabel("cat");
    double intensity=0.;
    if(channel=="SR") {
      variables.assignValueOnly(*row);
      for(int p=0;p<10;++p) intensity+=coefficients[p]*shapes[p]->getVal();
    } else intensity=cr.at(channel.back()-'0');
    if(!(intensity>0.) || !std::isfinite(intensity)) throw std::runtime_error("Invalid intensity in score derivative");
    values.push_back(std::log(intensity));
  }
  return values;
}
'''


def score_integral(weights, response):
    """Quadrature of a squared score, preserving signed integration masses."""
    weights=np.asarray(weights,dtype=float);response=np.asarray(response,dtype=float)
    if weights.shape!=response.shape or not np.all(np.isfinite(weights)) or not np.all(np.isfinite(response)):
        raise ValueError("Invalid score quadrature")
    terms=weights*response**2
    return float(terms.sum()), {"positive":float(terms[terms>=0].sum()),
                               "negative":float(terms[terms<0].sum())}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--bestfit",type=Path,required=True)
    p.add_argument("--covariance",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--workspace",default="combined")
    p.add_argument("--model-config",default="ModelConfig_4l_offshell")
    p.add_argument("--data",default="obsData_4l_offshell")
    p.add_argument("--weighted-asimov",action="store_true")
    p.add_argument("--allow-signed-quadrature",action="store_true",
                   help="Preserve signed integration masses; do not interpret negative rows as Poisson means")
    args=p.parse_args()
    source,w,mc,data,parameters,floating,nll=load_model(args)
    best=json.loads(args.bestfit.read_text())
    assert best["input_sha256"]==sha256(args.input)
    restore(w,best["parameters"])
    covariance=json.loads(args.covariance.read_text())
    assert covariance["passed"] and covariance["input_sha256"]==best["input_sha256"]
    point={r["name"]:r["value"] for r in covariance["parameters"]}
    assert set(point)=={v.GetName() for v in floating}
    for name,value in point.items(): w.var(name).setVal(value)
    # Inspect the actual constraint graph; only independent unit Gaussians
    # permit assigning the covariance column itself as a unit-GO response.
    auxiliaries={}
    for pdf in mc.GetPdf().getComponents():
        if not pdf.InheritsFrom("RooGaussian"): continue
        names={pdf.getX().GetName(),pdf.getMean().GetName()}
        alpha=[n for n in names if n.startswith("alpha_")]
        if not alpha: continue
        assert len(alpha)==1
        a=alpha[0];g="nom_"+a[6:]
        assert names=={a,g} and abs(pdf.getSigma().getVal()-1)<1e-14
        assert pdf.getSigma().isConstant() and w.var(g).isConstant()
        assert a not in auxiliaries
        auxiliaries[a]={"global":g,"precision":1.,"value":w.var(g).getVal()}
    assert set(auxiliaries)=={n for n in point if n.startswith("alpha_")}
    ROOT.gInterpreter.Declare(CPP)
    channels=[];weights=[]
    for i in range(data.numEntries()):
        row=data.get(i);channels.append(row.getCatLabel("cat"));weights.append(data.weight())
    channels=np.array(channels);weights=np.array(weights)
    metadata=json.loads(str(source.Get("asimov_metadata").GetString()))
    assert np.all(np.isfinite(weights)) and len(weights)==metadata["reference_events"]+3
    assert np.count_nonzero(channels=="SR")==metadata["reference_events"]
    for channel in ("CR0","CR1","CR2"):
        assert np.count_nonzero(channels==channel)==1 and np.all(weights[channels==channel]>0)
    negative=weights<0
    if np.any(negative) and not args.allow_signed_quadrature:
        raise RuntimeError("Signed quadrature requires --allow-signed-quadrature")
    assert float(weights[channels=="SR"].sum())>0
    names=covariance["parameter_order"]
    V=np.array(covariance["covariance"]);mu=names.index("mu")
    # Only one POI is requested. The covariance column is the probe: its
    # directional log-intensity derivative is exactly the complete J_mu,i.
    direction=V[:,mu]
    started=time.monotonic()
    derivatives=[]
    for step in (1e-4,5e-5):
        for name,shift in zip(names,direction):w.var(name).setVal(point[name]+step*shift)
        plus=np.array(ROOT.offshell_log_intensities(w,data))
        for name,shift in zip(names,direction):w.var(name).setVal(point[name]-step*shift)
        minus=np.array(ROOT.offshell_log_intensities(w,data))
        derivatives.append((plus-minus)/(2*step))
    for name,value in point.items():w.var(name).setVal(value)
    response=(4*derivatives[1]-derivatives[0])/3
    components={};differences={};contributions={}
    for channel in ("SR","CR0","CR1","CR2"):
        mask=channels==channel
        variance,contributions[channel]=score_integral(weights[mask],response[mask])
        coarse,_=score_integral(weights[mask],derivatives[0][mask])
        components[channel]=variance
        differences[channel]=abs(variance-coarse)/max(abs(variance),1e-30)
    nll_at_point=nll.getVal()
    result={"input_sha256":best["input_sha256"],"bestfit_sha256":sha256(args.bestfit),
        "covariance_sha256":sha256(args.covariance),"parameters":point,"parameter_order":names,
        "auxiliaries":auxiliaries,"primary_variance":components,
        "primary_variance_convention":("Signed quadrature approximation to the population Poisson score integral; negative rows are integration coefficients, not negative Poisson variances. No absolute-value replacement, clipping or SumW2 correction."
            if np.any(negative) else "Independent Poisson observation fluctuations; the fixed positive Asimov row weights supply the generating mean measure. No SumW2 correction."),
        "quadrature":{"reference_events":metadata["reference_events"],"negative_rows":int(negative.sum()),
            "negative_weight_sum":float(weights[negative].sum()),"positive_weight_sum":float(weights[~negative].sum()),
            "signed_variance_contributions":contributions},
        "difference_method":"Directional log-intensity derivative along V[:,mu], with h=1e-4 and h/2; Richardson extrapolation",
        "step_relative_variance_differences":differences,"nll":nll_at_point,
        "seconds":time.monotonic()-started,"slurm_job_id":os.environ.get("SLURM_JOB_ID"),
        "passed":bool(max(differences.values())<1e-5 and
            all(np.isfinite(v) and v>=0 for v in components.values()) and
            abs(nll_at_point-covariance["nll"])<1e-6)}
    args.output.mkdir(parents=True,exist_ok=True)
    write_json(args.output/"observation_response.json",result)
    source.Close()
    if not result["passed"]:raise RuntimeError("Score derivatives failed step-size stability")


if __name__=="__main__":main()
