#!/usr/bin/env python3
"""Export validated covariance source tables and ranking into the presentation."""
import argparse
import json
from pathlib import Path
import shutil


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results",type=Path,required=True)
    p.add_argument("--analysis-commit",required=True)
    p.add_argument("--output",type=Path,default=Path(__file__).resolve().parents[1])
    args=p.parse_args()
    result=json.loads((args.results/"report/summary.json").read_text())
    resources=json.loads((args.results/"resources.json").read_text())
    assert result["complete"] and resources["complete"]
    assert all(s["curvature_validation"]["passed"] for s in result["samples"].values())
    result.update(resources=resources,analysis_commit=args.analysis_commit)
    encoded=json.dumps(result,indent=2,allow_nan=False)
    (args.output/"data/decomposition_snapshot.json").write_text(encoded+"\n")
    (args.output/"data/decomposition_snapshot.js").write_text("window.DECOMPOSITION_RESULTS = "+encoded+";\n")
    for sample in ("true","data_like"):
        for kind in ("uncertainties","impacts","covariance"):
            name=f"{sample}_{kind}.csv";shutil.copyfile(args.results/"report"/name,args.output/"data"/name)
    shutil.copyfile(args.results/"report/auxiliary_impacts.svg",args.output/"assets/auxiliary_impacts.svg")
    for name in ("all_np_uncertainty_ranking.csv", "all_np_uncertainty_ranking_README.md"):
        shutil.copyfile(args.results/"report"/name,args.output/"data"/name)
    shutil.copyfile(args.results/"resources.csv",args.output/"data/decomposition_resources.csv")
    print("Exported complete covariance decomposition")


if __name__=="__main__":main()
