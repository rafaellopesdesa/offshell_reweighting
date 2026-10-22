#!/usr/bin/env python3
"""Collect measured accounting for the new 4096 work, retaining failed checks."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess

from collect_anchor_resources import seconds, memory_gib
from report_shifted_asimov import write_csv


def main(args):
    manifest=json.loads((args.results/"jobs.json").read_text())
    extension=json.loads((args.results/"true_extension/jobs.json").read_text())
    recovery=json.loads((args.results/"data_like/fit/global_validation.json").read_text())
    assert recovery["passed"]
    specs=[]
    for j in extension["jobs"]:
        specs.append({"job_id":j["job_id"],"sample":"true","phase":"upper_extension","kind":j["kind"]})
    for j in manifest["jobs"]:
        if j["kind"] in ("superseded_pending_finalizer", "initial_scan", "initial_data_like_finalize"):continue
        sample="true" if j["kind"]=="true_finalize" else "data_like"
        phase="global_workflow" if j["kind"] in ("build_global_validation","global_recovery_validation","initial_global_recovery_validation","global_refit") else "scan" if j["kind"]=="scan" else "finalize"
        specs.append({"job_id":j["job_id"],"sample":sample,"phase":phase,"kind":j["kind"]})
    fields=["JobID","JobName","State","ElapsedRaw","TotalCPU","CPUTimeRAW","MaxRSS","ReqMem","AllocCPUS","Start","End","NodeList","ExitCode"]
    cmd=["sacct","-n","-P","--units=K","-j",",".join(j["job_id"] for j in specs),"--format="+",".join(fields)]
    raw=subprocess.check_output(cmd,text=True,env={**os.environ,"TZ":"UTC"})
    records=[dict(zip(fields,line.split("|"))) for line in raw.splitlines() if line.strip()]
    lookup={r["JobID"]:r for r in records};jobs=[]
    for spec in specs:
        row=lookup[spec["job_id"]]
        exception=None
        if row["State"]!="COMPLETED" or row["ExitCode"]!="0:0":
            if spec["kind"]=="build_global_validation" and row["State"]=="FAILED":
                log=(args.results/f'data_like/global_{spec["job_id"]}.log').read_text()
                if "Saved RooFitResult has invalid convergence diagnostics" not in log:
                    raise RuntimeError("Unexpected global workflow failure")
                exception="Fits completed; selected stationary-start covariance failed validation. The displaced start was repeated separately and its saved result validated."
            elif spec["kind"]=="initial_global_recovery_validation" and row["State"]=="FAILED":
                log=(args.results/f'global_recovery_{spec["job_id"]}.log').read_text()
                if "Covariance matrix is not positive definite" not in log:
                    raise RuntimeError("Unexpected recovery-check failure")
                exception="The checker rejected a non-positive-definite covariance; a separate displaced-start fit supplies the validated scan baseline."
            else:raise RuntimeError(f"Unfinished or failed required job: {row}")
        peaks=[memory_gib(r["MaxRSS"]) for r in records if r["JobID"].startswith(spec["job_id"]+".") and r["MaxRSS"]]
        cpu=seconds(row["TotalCPU"])
        if not peaks or cpu<=0:raise RuntimeError("Accounting is not final yet")
        jobs.append({**spec,"state":row["State"],"exit_code":row["ExitCode"],"exception":exception,
            "elapsed_seconds":int(row["ElapsedRaw"]),"total_cpu_seconds":cpu,"peak_rss_gib":max(peaks),
            "allocated_cpus":int(row["AllocCPUS"]),"requested_memory":row["ReqMem"],
            "start_utc":row["Start"],"end_utc":row["End"],"node":row["NodeList"]})
    phases=[]
    for sample,phase in [("true","upper_extension"),("true","finalize"),("data_like","global_workflow"),("data_like","scan"),("data_like","finalize")]:
        group=[r for r in jobs if (r["sample"],r["phase"])==(sample,phase)]
        phases.append({"sample":sample,"phase":phase,"job_count":len(group),
            "wall_seconds":(max(datetime.fromisoformat(r["end_utc"]) for r in group)-min(datetime.fromisoformat(r["start_utc"]) for r in group)).total_seconds(),
            "total_cpu_seconds":sum(r["total_cpu_seconds"] for r in group),
            "largest_job_peak_rss_gib":max(r["peak_rss_gib"] for r in group)})
    accounting="|".join(fields)+"\n"+raw
    report={"complete":True,"created_utc":datetime.now(timezone.utc).isoformat(),"jobs":jobs,"phases":phases,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "accounting_sha256":hashlib.sha256(accounting.encode()).hexdigest(),
        "definitions":{"wall_seconds":"First start to last finish within each phase, including intervening scheduling gaps",
            "total_cpu_seconds":"Summed Slurm TotalCPU, without double-counting steps",
            "largest_job_peak_rss_gib":"Largest single job's step MaxRSS; not summed concurrent memory",
            "scope":"New work only. Earlier true-Asimov construction, original 31-point scan and zero attempt are reported separately in the preserved experiment. Prototype development, installation and staging are excluded.",
            "configuration":"One CPU/job; true extension/finalizer 2 GiB; data-like build/global 6 GiB, recovery check/refit 3 GiB, five parallel profile chunks and finalizer 4 GiB. Different nodes; no controlled hardware comparison.",
            "covariance":"No explicit full Hesse. All fitting time, both covariance-check failures and the displaced-start repeat are retained in global-workflow costs."}}
    (args.results/"slurm_accounting.txt").write_text(accounting)
    (args.results/"resources.json").write_text(json.dumps(report,indent=2)+"\n")
    write_csv(args.results/"resources.csv",phases,list(phases[0]))
    print(json.dumps(phases,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results",type=Path,required=True)
    main(p.parse_args())
