#!/usr/bin/env python3
"""Append the validated 4096 data-like study to the frozen presentation."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results",type=Path,required=True)
    p.add_argument("--analysis-commit",required=True)
    p.add_argument("--output",type=Path,default=Path(__file__).resolve().parents[1])
    args=p.parse_args()
    comparison=args.results/"comparison"
    summary=json.loads((comparison/"summary.json").read_text())
    resources=json.loads((args.results/"resources.json").read_text())
    timing=json.loads((args.results/"timing/summary.json").read_text())
    assert summary["complete"] and resources["complete"]
    assert all(s["points"]==42 for s in summary["samples"].values())
    assert timing["complete"]
    snapshot={**summary,"resources":resources,"timing":timing,"analysis_commit":args.analysis_commit,
              "created_utc":datetime.now(timezone.utc).isoformat(),
              "export_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output/"data/shifted_snapshot.json").write_text(json.dumps(snapshot,indent=2,allow_nan=False)+"\n")
    (args.output/"data/shifted_snapshot.js").write_text("window.SHIFTED_RESULTS = "+json.dumps(snapshot,allow_nan=False)+";\n")
    for name in ("true_scan.csv","data_like_scan.csv","true_parameters.csv","data_like_parameters.csv",
                 "true_validation.csv","data_like_validation.csv","comparison.csv","generation_parameters.csv"):
        shutil.copyfile(comparison/name,args.output/"data"/("shifted_"+name))
    shutil.copyfile(args.results/"resources.csv",args.output/"data/shifted_resources.csv")
    shutil.copyfile(args.results/"timing/summary.json",args.output/"data/normalizer_timing.json")
    shutil.copyfile(comparison/"true_vs_data_like.svg",args.output/"assets/true_vs_data_like.svg")
    print(json.dumps({"complete":True,"samples":list(summary["samples"]),"points_per_scan":42}))


if __name__=="__main__":main()
