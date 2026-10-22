#!/usr/bin/env python3
"""Account for the covariance study, including failed setup and superseded diagnostics."""
import argparse
import csv
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import subprocess

from collect_anchor_resources import seconds,memory_gib


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--results",type=Path,required=True)
    args=p.parse_args();manifest=json.loads((args.results/"jobs.json").read_text())
    specs=[{**row,"phase":phase} for phase in ("initial_jobs","jobs","superseded_response_jobs","probe_jobs") for row in manifest[phase]]
    fields=["JobID","State","ExitCode","ElapsedRaw","TotalCPU","MaxRSS","Start","End","NodeList","ReqMem"]
    cmd=["sacct","-n","-P","--units=K","-j",",".join(r["job_id"] for r in specs),"--format="+",".join(fields)]
    raw=subprocess.check_output(cmd,text=True,env={**os.environ,"TZ":"UTC"})
    rows=[dict(zip(fields,line.split("|"))) for line in raw.splitlines() if line.strip()]
    lookup={r["JobID"]:r for r in rows};jobs=[]
    for spec in specs:
        row=lookup[spec["job_id"]]
        if spec["phase"] in ("jobs","probe_jobs"):
            assert row["State"]=="COMPLETED" and row["ExitCode"]=="0:0",row
        elif spec["phase"]=="initial_jobs":assert row["State"]=="FAILED",row
        else:assert row["State"].startswith("CANCELLED"),row
        peaks=[memory_gib(r["MaxRSS"]) for r in rows if r["JobID"].startswith(spec["job_id"]+".") and r["MaxRSS"]]
        assert peaks and seconds(row["TotalCPU"])>0
        jobs.append({"model":spec["model"],"phase":spec["phase"],"job_id":spec["job_id"],
            "state":row["State"],"elapsed_seconds":int(row["ElapsedRaw"]),"cpu_seconds":seconds(row["TotalCPU"]),
            "peak_rss_gib":max(peaks),"start":row["Start"],"end":row["End"],"node":row["NodeList"],"requested_memory":row["ReqMem"]})
    totals=[]
    for model in ("true","data_like"):
        group=[r for r in jobs if r["model"]==model]
        totals.append({"model":model,"wall_seconds":(max(datetime.fromisoformat(r["end"]) for r in group)-min(datetime.fromisoformat(r["start"]) for r in group)).total_seconds(),
            "cpu_seconds":sum(r["cpu_seconds"] for r in group),"largest_job_peak_rss_gib":max(r["peak_rss_gib"] for r in group)})
    result={"complete":True,"created_utc":datetime.now(timezone.utc).isoformat(),"jobs":jobs,"totals":totals,
        "scope":"Separate covariance study: setup failures, full Hesse with Migrad initialization, superseded full-Jacobian diagnostics, covariance-column probes and directional curvature checks are all included.",
        "definitions":"Wall time spans first job start to last end per model; CPU sums parent-job usage; peak RSS is the largest individual job step, not concurrent total memory. One CPU and 4 GiB requested per job."}
    (args.results/"resources.json").write_text(json.dumps(result,indent=2)+"\n")
    (args.results/"slurm_accounting.txt").write_text("|".join(fields)+"\n"+raw)
    with (args.results/"resources.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(totals[0]),lineterminator="\n");writer.writeheader();writer.writerows(totals)
    print(json.dumps(totals,indent=2))


if __name__=="__main__":main()
