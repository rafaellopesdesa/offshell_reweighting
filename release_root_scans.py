#!/usr/bin/env python3
"""Release pending Slurm scans from fit checkpoints without interrupting fits.

Reads fit_jobs.json, scan_jobs.json, and finalize_jobs.json from a run directory.
Finalization continues to depend on the completed fit, so provisional scan
baselines are always replaced by the final continuous minimum.
"""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import time


def run(*command):
    return subprocess.check_output(command, text=True).strip()


def save(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def states(ids):
    output = run("squeue", "--noheader", "--user", os.environ["USER"], "--format=%i|%T")
    return {job: state for job, state in (line.split("|", 1) for line in output.splitlines() if line)
            if job in ids}


def release(directory, size, fit_job):
    checkpoint = directory/f"n_{size}"/"bestfit_checkpoint.json"
    final = directory/f"n_{size}"/"bestfit.json"
    if final.exists():
        return "Final fit already available; existing dependencies retained"
    if states([fit_job]).get(fit_job) != "RUNNING":
        return None
    if not checkpoint.exists():
        return None
    best = json.loads(checkpoint.read_text())
    if len(best["parameters"]) != 131 or any(a["migrad_status"] for a in best["attempts"]):
        raise RuntimeError("Invalid global-fit checkpoint")
    scans = json.loads((directory/"scan_jobs.json").read_text())
    old = [job for job in scans if job["size"] == int(size)]
    ids = [job["job_id"] for job in old]
    if fit_job in ids:
        raise RuntimeError("A fit job must never appear in the scan job list")
    if set(states(ids).values()) != {"PENDING"} or len(states(ids)) != len(ids):
        return "Existing scans already active; dependencies retained"
    run("scontrol", "hold", ",".join(ids))
    replacements = []
    original_finals = json.loads((directory/"finalize_jobs.json").read_text())
    previous_dependency = next(arg.split("=", 1)[1] for arg in original_finals[size]["command"]
                               if arg.startswith("--dependency="))
    try:
        if set(states(ids).values()) != {"PENDING"}:
            run("scontrol", "release", ",".join(ids))
            return "A scan started during the check; existing work retained"
        for job in old:
            command = [str(checkpoint.resolve()) if arg == str(final.resolve()) else arg
                       for arg in job["command"] if not arg.startswith("--dependency=")]
            job_id = run(*command).split(";", 1)[0]
            replacements.append({**job, "job_id": job_id, "command": command,
                                 "checkpoint_baseline": str(checkpoint.resolve())})
        finals = json.loads((directory/"finalize_jobs.json").read_text())
        dependency = "afterok:" + ":".join([fit_job]+[job["job_id"] for job in replacements])
        run("scontrol", "update", "JobId="+finals[size]["job_id"], "Dependency="+dependency)
        finals[size]["command"] = ["--dependency="+dependency if arg.startswith("--dependency=") else arg
                                    for arg in finals[size]["command"]]
        save(directory/"finalize_jobs.json", finals)
        save(directory/f"scan_jobs_before_checkpoint_{size}.json", old)
        save(directory/"scan_jobs.json", [job for job in scans if job["size"] != int(size)]+replacements)
        run("scancel", *ids)  # These are held, pending scans, never the fit.
    except Exception:
        if replacements:
            run("scancel", *[job["job_id"] for job in replacements])
        run("scontrol", "update", "JobId="+original_finals[size]["job_id"],
            "Dependency="+previous_dependency)
        save(directory/"finalize_jobs.json", original_finals)
        save(directory/"scan_jobs.json", scans)
        run("scontrol", "release", ",".join(ids))
        raise
    return {"replaced_pending_scans": ids, "released_scans": [job["job_id"] for job in replacements],
            "fit_job_untouched": fit_job, "finalization_dependency": dependency,
            "checkpoint_mu_hat": best["mu_hat"], "checkpoint_nll": best["nll_min"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run this lightweight controller in a Slurm allocation")
    fits = json.loads((args.run/"fit_jobs.json").read_text())
    completed = {}
    while len(completed) < len(fits):
        for size, job in fits.items():
            if size in completed:
                continue
            current = states([job["job_id"]])
            if job["job_id"] not in current:
                result = "Fit left the queue; existing scan dependencies retained"
            else:
                result = release(args.run, size, job["job_id"])
            if result is not None:
                completed[size] = {"utc": datetime.now(timezone.utc).isoformat(), "result": result}
                save(args.run/"scan_release_history.json", completed)
                print(size, json.dumps(completed[size]), flush=True)
        if len(completed) < len(fits):
            time.sleep(30)


if __name__ == "__main__":
    main()
