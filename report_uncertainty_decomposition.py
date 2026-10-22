#!/usr/bin/env python3
"""Source-covariance propagation and an impact-style auxiliary ranking.

Implements arXiv:2608.15312, Eqs. (22), (30), (44), (51), (81).
Observed-Hessian closure is measured, never enforced by rescaling sources.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from export_np_ranking import write_ranking


def group(name):
    short=name.removeprefix("alpha_ATLAS_")
    if short.startswith("JET_"):return "Jets"
    if short.startswith(("EG_","EL_")):return "Electrons"
    if short.startswith("MUON_"):return "Muons"
    if short=="alpha_LUMI":return "Luminosity"
    if short=="PRW_DATASF":return "Pileup"
    if "PDF" in short:return "Parton distributions"
    if short.startswith(("PS_","H4l_Shower_")):return "Parton shower"
    if short.startswith(("HOEW","HOQCD","QCD_")):return "Higher-order theory"
    raise ValueError("Unassigned nuisance: "+name)


def csv_write(path,rows):
    with path.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(rows)


def analyze(directory):
    h=json.loads((directory/"covariance.json").read_text())
    b=json.loads((directory/"observation_response.json").read_text())
    check=json.loads((directory/"curvature_validation.json").read_text())
    assert check["passed"]
    covariance_hash=hashlib.sha256((directory/"covariance.json").read_bytes()).hexdigest()
    assert b["covariance_sha256"]==check["covariance_sha256"]==covariance_hash
    assert h["passed"] and b["passed"] and h["input_sha256"]==b["input_sha256"]
    assert h["parameter_order"]==b["parameter_order"]
    names=h["parameter_order"];V=np.array(h["covariance"]);mu=names.index("mu")
    values={r["name"]:r["value"] for r in h["parameters"]}
    assert max(abs(values[n]-b["parameters"][n]) for n in names)<1e-12
    assert np.linalg.eigvalsh(V)[0]>0
    row=V[mu];total=float(V[mu,mu]);sources=[];impacts=[]
    for i,name in enumerate(names):
        if name=="mu":continue
        sigma=float(np.sqrt(V[i,i]));constrained=name in b["auxiliaries"]
        response=float(row[i]) if constrained else None
        nominal=b["auxiliaries"][name]["value"] if constrained else None
        impacts.append({"name":name,"constrained":constrained,"group":group(name) if constrained else "Unconstrained qqZZ",
            "value":values[name],"anchor":nominal,"pull_from_anchor":values[name]-nominal if constrained else None,
            "posterior_sigma":sigma,"unit_auxiliary_response":response,
            "auxiliary_variance":response**2 if constrained else None,
            "conditional_postfit_np_impact":float(row[i]/sigma),
            "correlation_with_mu":float(row[i]/np.sqrt(total*V[i,i]))})
        if constrained:sources.append({"name":name,"group":group(name),"variance":response**2,"sigma":abs(response)})
    auxiliary=sum(s["variance"] for s in sources)
    primary=b["primary_variance"]
    assert min(primary.values())>=0
    statistical=sum(primary.values());propagated=statistical+auxiliary
    profile_curvature=next(r["richardson_curvature"] for r in check["checks"] if r["direction"]=="profile_mu_tangent")
    a=[i for i,n in enumerate(names) if n in b["auxiliaries"]]
    q=[i for i,n in enumerate(names) if n.startswith("mu_qqZZ")]
    penalty=lambda ids:float(row[ids]@np.linalg.solve(V[np.ix_(ids,ids)],row[ids]))
    grouped=[]
    for label in dict.fromkeys(s["group"] for s in sources):
        variance=sum(s["variance"] for s in sources if s["group"]==label)
        grouped.append({"source":label,"kind":"auxiliary","variance":variance,"sigma":float(np.sqrt(variance)),"fraction_of_hessian_variance":variance/total})
    grouped.sort(key=lambda r:r["variance"],reverse=True)
    for channel,variance in primary.items():
        grouped.append({"source":channel+" statistics","kind":"primary","variance":variance,"sigma":float(np.sqrt(variance)),"fraction_of_hessian_variance":variance/total})
    grouped.extend([{"source":label,"kind":"subtotal","variance":variance,"sigma":float(np.sqrt(variance)),"fraction_of_hessian_variance":variance/total}
                    for label,variance in [("Auxiliary subtotal",auxiliary),("Primary-data subtotal",statistical),("Propagated source total",propagated),("Observed-Hessian total",total)]])
    out={"mu_hat":values["mu"],"hessian_variance":total,"hessian_sigma":float(np.sqrt(total)),
        "auxiliary_variance":auxiliary,"auxiliary_sigma":float(np.sqrt(auxiliary)),
        "primary_variance":statistical,"primary_sigma":float(np.sqrt(statistical)),
        "primary_components":primary,"propagated_variance":propagated,"propagated_sigma":float(np.sqrt(propagated)),
        "closure_absolute":propagated-total,"closure_relative":propagated/total-1,
        "native_covariance_probe_curvature":profile_curvature,
        "source_to_native_curvature_probe_relative_difference":propagated/(total*profile_curvature)-1,
        "curvature_validation":check,
        "groups":grouped,"impacts":impacts,
        "conditional_diagnostics":{"qqzz_block_profiling_variance_penalty":penalty(q),
            "alpha_block_profiling_variance_penalty":penalty(a),
            "alpha_fixed_sigma":float(np.sqrt(total-penalty(a))),
            "note":"Freeze/float penalties are conditional diagnostics, not additive physical-source uncertainties."},
        "covariance_quality":h["covariance_quality"],"hesse_status":h["hesse_status"],
        "hesse_seconds":h["seconds"],"maximum_step_relative_difference":max(b["step_relative_variance_differences"].values()),
        "input_sha256":h["input_sha256"],
        "source_sha256":{name:hashlib.sha256((directory/name).read_bytes()).hexdigest() for name in ("covariance.json","observation_response.json","curvature_validation.json")}}
    return out


def plot(samples,output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"svg.fonttype":"none","svg.hashsalt":"offshell-uncertainty","axes.spines.top":False,"axes.spines.right":False})
    records={key:{r["name"]:r for r in s["impacts"] if r["constrained"]} for key,s in samples.items()}
    names=sorted(records["true"],key=lambda n:max(abs(records[k][n]["unit_auxiliary_response"]) for k in records),reverse=True)[:12]
    fig,(pull,impact)=plt.subplots(1,2,figsize=(12,6.2),sharey=True,gridspec_kw={"width_ratios":[1,1.6],"wspace":.10})
    for key,label,color,offset in [("true","True Asimov","#126c96",-.15),("data_like","Data-like Asimov","#c95639",.15)]:
        rows=[records[key][n] for n in names];y=np.arange(len(names))+offset
        pull.errorbar([r["pull_from_anchor"] for r in rows],y,xerr=[r["posterior_sigma"] for r in rows],fmt="o",ms=3,color=color,capsize=2,label=label)
        positive=np.array([r["unit_auxiliary_response"] for r in rows])
        impact.barh(y,positive,height=.25,color=color,alpha=.85)
        impact.barh(y,-positive,height=.25,facecolor="none",edgecolor=color,linewidth=1.1)
    for axis in (pull,impact):
        axis.axvline(0,color="#66727b",lw=.8);axis.grid(axis="x",alpha=.16)
    pull.set_yticks(np.arange(len(names)),[n.removeprefix("alpha_ATLAS_") for n in names],fontsize=9)
    pull.invert_yaxis();pull.set_xlim(-1.3,1.3)
    pull.set_xlabel(r"$\hat\alpha-a\;\pm\;\sigma_{\alpha,\mathrm{post}}$")
    pull.legend(loc="lower left",bbox_to_anchor=(0,1.01),frameon=False,fontsize=9)
    impact.set_xlabel(r"Local $\Delta\hat\mu$ from a $\pm1\sigma$ auxiliary shift")
    impact.set_title("Filled: +1σ anchor shift · outline: −1σ",fontsize=10,pad=15)
    fig.suptitle("4096 events · largest auxiliary-source responses",fontsize=14,x=.63,y=.985)
    fig.subplots_adjust(left=.32,right=.98,bottom=.12,top=.88)
    for ext in ("svg","pdf","png"):
        fig.savefig(output/f"auxiliary_impacts.{ext}",dpi=180,metadata={"Date":None} if ext=="svg" else None)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    samples={name:analyze(args.results/name) for name in ("true","data_like")}
    result={"complete":True,"created_utc":datetime.now(timezone.utc).isoformat(),
        "reference":"https://arxiv.org/abs/2608.15312","information":"Observed full Hesse; all 131 parameters floating",
        "approximation":"Local linear response; not asymmetric profile intervals or finite one-sigma refits",
        "source_model":"Primary Poisson counting measure plus 127 independent unit Gaussian auxiliary observations",
        "closure_policy":"Report propagated and Hessian totals separately; do not force closure or rescale source components",
        "samples":samples}
    for name,s in samples.items():
        csv_write(args.output/f"{name}_uncertainties.csv",s["groups"])
        csv_write(args.output/f"{name}_impacts.csv",s["impacts"])
        covariance=json.loads((args.results/name/"covariance.json").read_text())
        with (args.output/f"{name}_covariance.csv").open("w",newline="") as f:
            writer=csv.writer(f,lineterminator="\n");writer.writerow(["parameter",*covariance["parameter_order"]])
            for parameter,row in zip(covariance["parameter_order"],covariance["covariance"]):writer.writerow([parameter,*row])
    (args.output/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    write_ranking(args.output/"summary.json", args.output)
    plot(samples,args.output)
    print(json.dumps({name:{k:v for k,v in s.items() if k not in ("groups","impacts","source_sha256")} for name,s in samples.items()},indent=2))


if __name__=="__main__":main()
