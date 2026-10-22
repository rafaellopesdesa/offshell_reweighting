#!/usr/bin/env python3
"""Export validated true/data-like 4096 profiles and exact HEPData comparisons."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

from compare_hepdata_scan import series_points, csv_series_points


def write_csv(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer=csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(args):
    output=args.output
    output.mkdir(parents=True, exist_ok=True)
    hashes={}
    def read(path):
        raw=path.read_bytes()
        hashes[str(path.relative_to(args.results))]=hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    base=args.results/"asimov_data_like_4096"
    table=read(args.results/"hepdata_comparison/table_2.json")
    expected=series_points(table,"Exp NSBI")
    assert expected==csv_series_points(args.results/"hepdata_comparison/table_2.csv","Exp NSBI")
    construction=read(base/"data_like/workspace/construction.json")
    check=read(base/"data_like/workspace/validation.json")
    recovery=read(base/"data_like/fit/global_validation.json")
    initial_fit=read(base/"data_like/fit_initial/bestfit.json")
    assert check["passed"] and recovery["passed"]
    result={"complete":True,"created_utc":datetime.now(timezone.utc).isoformat(),
        "reference_events":4096,"reference_doi":table["doi"],"reference_series":"Exp NSBI",
        "comparison_method":"Exact shared coordinates; continuous fitted baselines; no re-zeroing, rescaling or reference interpolation",
        "model_distinction":"True Asimov retains nominal/anchor normalization. Data-like shifts observed-fit nuisances AND normalizes the fully interpolated shape at every trial point; physical rates stay separate.",
        "generation":{k:construction[k] for k in ("generation_parameters","sr_asimov_yield","cr_asimov_counts","observed_bestfit_sha256","source_workspace_sha256","workspace_sha256")},
        "construction_validation":{ "maximum_absolute_score":max(c["maximum_absolute_score"] for c in check["generating_point_checks"]),
            "maximum_process_integral_error":max(c["maximum_process_integral_error"] for c in check["generating_point_checks"]),
            "maximum_normalizer_disagreement":max(max(c["maximum_cached_vs_direct_normalizer_error"],c["normalizer_error_after_gradient_updates"]) for c in check["generating_point_checks"]),
            "generating_point_cases":[c["case"] for c in check["generating_point_checks"]]},
        "recovery":recovery["recovery"],"samples":{}}
    grids={}
    for label,fit,final in [("true", args.results/"asimov_anchor_normalized/n_4096/fit",base/"true_extension/final"),
                            ("data_like",base/"data_like/fit",base/"data_like/final")]:
        best=read(fit/"bestfit.json");summary=read(final/"scan_summary.json")
        validation=read(final/"likelihood_validation.json");rows=read(final/"scan_points.json")
        assert best["complete"] and validation["passed"] and summary["all_converged"]
        assert summary["input_sha256"]==best["input_sha256"] and summary["baseline_nll"]==best["nll_min"]
        grid=[0.,.01]+[round(i*.05,2) for i in range(1,41)]
        assert [r["mu"] for r in rows]==grid and len(validation["scan_checks"])==42
        aggregate=[]
        for r in rows:
            assert r["migrad_status"]==0 and r["edm"]<=1e-4
            q=2*(r["nll"]-best["nll_min"])
            assert abs(q-r["q"])<1e-9 and q>=-2e-4
            reference=expected.get(r["mu"])
            aggregate.append({"mu":r["mu"],"q":q,"hepdata_q":reference,"delta_q":q-reference if reference is not None else None,
                              "nll":r["nll"],"edm":r["edm"],"migrad_status":r["migrad_status"],"seconds":r["seconds"]})
        grids[label]=aggregate
        metrics={}
        for key,low,high in [("full",0,2),("previous_range",.05,1.5),("positive_published",.05,2)]:
            common=[r for r in aggregate if r["delta_q"] is not None and low<=r["mu"]<=high]
            extreme=max(common,key=lambda r:abs(r["delta_q"]))
            metrics[key]={"range":[low,high],"points":len(common),"max_abs_delta_q":abs(extreme["delta_q"]),
                          "mu_at_maximum":extreme["mu"],"rms_delta_q":math.sqrt(sum(r["delta_q"]**2 for r in common)/len(common))}
        diagnostics=validation["scan_checks"]
        result["samples"][label]={"mu_hat":best["mu_hat"],"nll_min":best["nll_min"],"edm":best["edm"],
            "covariance_quality":best["covariance_quality"],"points":42,"q_zero":aggregate[0]["q"],"q_two":aggregate[-1]["q"],
            "max_profile_edm":summary["max_edm"],"max_overlap_delta_q":max(r["absolute_q_difference"] for r in summary["overlap_checks"]),
            "maximum_native_explicit_disagreement":max(abs(r["delta_nll_disagreement"]) for r in diagnostics),
            "maximum_sr_integral_mismatch":max(abs(r["sr_integral_relative_difference"]) for r in diagnostics),
            "minimum_fitted_intensity":min(r["minimum_current_intensity"] for r in diagnostics),
            "comparison":metrics,"attempts":best["attempts"]}
        if label=="data_like":
            result["samples"][label]["attempts"]=[
                {**r,"role":"initial diagnostic; covariance rejected" if i==0 else "initial displaced start"}
                for i,r in enumerate(initial_fit["attempts"])] + [
                {**r,"role":"saved profile baseline"} for r in best["attempts"]]
            result["samples"][label]["initial_fit_note"]="Initial stationary-start covariance was not positive definite. Its files are preserved in fit_initial; the displaced start was repeated to save a valid fit result."
        write_csv(output/f"{label}_scan.csv",aggregate,list(aggregate[0]))
        write_csv(output/f"{label}_validation.csv",diagnostics,list(diagnostics[0]))
        params=[dict(r,injected=(construction["generation_parameters"][r["name"]] if label=="data_like" else (0. if r["name"].startswith("alpha_") else 1.))) for r in best["parameters"]]
        write_csv(output/f"{label}_parameters.csv",params,["name","injected","value","error","min","max"])
    common=[{"mu":t["mu"],"true_q":t["q"],"data_like_q":d["q"],"hepdata_q":t["hepdata_q"],"data_like_minus_true":d["q"]-t["q"]}
            for t,d in zip(grids["true"],grids["data_like"])]
    write_csv(output/"comparison.csv",common,list(common[0]))
    write_csv(output/"generation_parameters.csv",[{"name":k,"value":v} for k,v in construction["generation_parameters"].items()],["name","value"])
    result["source_sha256"]=hashes
    (output/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    make_plot(grids,expected,result,output)
    print(json.dumps(result["samples"],indent=2))


def make_plot(grids,expected,summary,output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,
                         "axes.spines.right":False,"svg.fonttype":"none","svg.hashsalt":"offshell-data-like"})
    fig,(ax,res)=plt.subplots(2,1,figsize=(10.5,6.3),sharex=True,gridspec_kw={"height_ratios":[3.2,1],"hspace":.08})
    ref=sorted((x,y) for x,y in expected.items() if 0<=x<=2)
    ax.plot(*zip(*ref),color="#253746",lw=1.7,label="HEPData expected NSBI (4ℓ, stat+syst)")
    for key,color,label in [("true","#126c96","True Asimov · original anchor normalization"),
                            ("data_like","#c95639","Data-like Asimov · full nuisance normalization")]:
        rows=grids[key]
        ax.plot([r["mu"] for r in rows],[r["q"] for r in rows],"o-",ms=3.4,lw=1.5,color=color,label=label)
        ax.plot(summary["samples"][key]["mu_hat"],0,"*",ms=11,color=color)
        shared=[r for r in rows if r["delta_q"] is not None]
        res.plot([r["mu"] for r in shared],[r["delta_q"] for r in shared],"o-",ms=3.2,lw=1.2,color=color)
    ax.set_ylabel(r"$q(\mu)=-2\ln\lambda(\mu)$")
    ax.set_ylim(bottom=-.04)
    ax.legend(frameon=False,fontsize=9.5,loc="upper center")
    ax.set_title("4,096 fixed reference events · generating μ = 1",fontsize=13,pad=12)
    res.axhline(0,color="#82909c",lw=1)
    res.set_ylabel("Local − HEPData",fontsize=10)
    res.set_xlabel(r"Off-shell signal strength $\mu$")
    res.set_xlim(0,2)
    for axis in (ax,res):axis.grid(alpha=.16)
    fig.subplots_adjust(left=.12,right=.98,bottom=.12,top=.92)
    for extension in ("svg","pdf","png"):
        fig.savefig(output/f"true_vs_data_like.{extension}",dpi=180,
                    metadata={"Date":None} if extension=="svg" else None)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
