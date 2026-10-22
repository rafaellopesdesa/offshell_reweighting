#!/usr/bin/env python3
"""Collect final Slurm time and memory accounting for the corrected experiment."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess


def seconds(value):
    days, clock = value.split("-", 1) if "-" in value else ("0", value)
    result = 0.0
    for part in clock.split(":"):
        result = 60*result+float(part)
    return 86400*int(days)+result


def memory_gib(value):
    if not value:
        return None
    unit = value[-1].upper()
    return float(value[:-1])*{"K": 2**10, "M": 2**20, "G": 2**30, "T": 2**40}[unit]/2**30


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path, required=True)
    args = p.parse_args()
    global_jobs = json.loads((args.results/"jobs.json").read_text())["jobs"]
    scans = json.loads((args.results/"scan_jobs.json").read_text())["jobs"]
    specs = [{"job_id": j["job_id"], "size": int(size), "phase": "build_global_validate"}
             for size,j in global_jobs.items()]
    specs += [{"job_id": j["job_id"], "size": j["size"],
               "phase": "scan" if j["kind"].startswith("scan_") else j["kind"]} for j in scans]
    fields = ["JobID", "JobName", "State", "ElapsedRaw", "TotalCPU", "CPUTimeRAW",
              "MaxRSS", "ReqMem", "AllocCPUS", "Start", "End", "NodeList", "ExitCode"]
    command = ["sacct", "-n", "-P", "--units=K", "-j", ",".join(j["job_id"] for j in specs),
               "--format="+",".join(fields)]
    raw = subprocess.check_output(command, text=True, env={**os.environ, "TZ": "UTC"})
    records = [dict(zip(fields, line.split("|"))) for line in raw.splitlines() if line.strip()]
    lookup = {r["JobID"]:r for r in records}
    jobs = []
    for spec in specs:
        row = lookup[spec["job_id"]]
        if row["State"] != "COMPLETED" or row["ExitCode"] != "0:0":
            raise RuntimeError(f"Wait for successful completion of {spec['job_id']}: {row['State']}")
        peaks = [memory_gib(r["MaxRSS"]) for r in records
                 if r["JobID"].startswith(spec["job_id"]+".") and r["MaxRSS"]]
        if not peaks:
            raise RuntimeError(f"Final memory accounting is not available yet for {spec['job_id']}")
        cpu_seconds = seconds(row["TotalCPU"])
        if cpu_seconds <= 0:
            raise RuntimeError(f"Final CPU accounting is not available yet for {spec['job_id']}")
        jobs.append({**spec, "elapsed_seconds": int(row["ElapsedRaw"]),
                     "total_cpu_seconds": cpu_seconds,
                     "allocated_cpu_seconds": int(row["CPUTimeRAW"]),
                     "peak_rss_gib": max(peaks), "requested_memory": row["ReqMem"],
                     "allocated_cpus": int(row["AllocCPUS"]), "start_utc": row["Start"],
                     "end_utc": row["End"], "node": row["NodeList"], "state": row["State"]})
    phases = []
    for size in (4096,8192):
        for phase in ("build_global_validate", "scan", "finalize"):
            group = [r for r in jobs if r["size"] == size and r["phase"] == phase]
            starts = [datetime.fromisoformat(r["start_utc"]) for r in group]
            ends = [datetime.fromisoformat(r["end_utc"]) for r in group]
            phases.append({"reference_events": size, "phase": phase, "job_count": len(group),
                           "wall_seconds": (max(ends)-min(starts)).total_seconds(),
                           "total_cpu_seconds": sum(r["total_cpu_seconds"] for r in group),
                           "sum_job_elapsed_seconds": sum(r["elapsed_seconds"] for r in group),
                           "largest_job_peak_rss_gib": max(r["peak_rss_gib"] for r in group)})
    accounting_text = "|".join(fields)+"\n"+raw
    report = {"created_utc": datetime.now(timezone.utc).isoformat(), "complete": True,
              "collector_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "accounting_command": command, "accounting_sha256": hashlib.sha256(accounting_text.encode()).hexdigest(),
              "jobs": jobs, "phases": phases,
              "definitions": {"wall_seconds": "First job start to last job end in the phase; excludes waiting before the first start",
                              "total_cpu_seconds": "Sum of Slurm TotalCPU for parent jobs, without double-counting job steps",
                              "largest_job_peak_rss_gib": "Largest Slurm-reported MaxRSS among the individual jobs' steps; not summed concurrent memory",
                              "configuration": "One CPU and 4 GiB requested per job; four scan chunks per sample run in parallel on different CPU nodes",
                              "global_stage": "Includes construction from already staged inputs, construction checks, two unconditional starts, serialization, and independent global validation; no explicit full Hesse",
                              "excluded_shared_setup": "ROOT installation and earlier input staging are excluded; the shared final HEPData comparison job is excluded from per-sample stage totals",
                              "finalize_stage": "Merge scan chunks, make figures, and independently validate the global minimum and all 31 profile points"},
              "documentation": "https://slurm.schedmd.com/sacct.html"}
    (args.results/"slurm_accounting.txt").write_text(accounting_text)
    (args.results/"resources.json").write_text(json.dumps(report, indent=2)+"\n")
    with (args.results/"resources.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(phases[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(phases)
    print(json.dumps(phases, indent=2))


if __name__ == "__main__":
    main()
