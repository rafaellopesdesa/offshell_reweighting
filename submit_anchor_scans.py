#!/usr/bin/env python3
"""Submit new anchor-normalized profiles only after both global jobs finish."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess


def run(*args):
    return subprocess.check_output(args, text=True).strip()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--setup", type=Path, required=True)
    args = p.parse_args()
    repo = Path(__file__).resolve().parent
    base = args.results.resolve()
    manifest = base/"scan_jobs.json"
    if manifest.exists():
        raise RuntimeError("A scan submission manifest already exists; inspect it before resubmitting")
    if run("git", "-C", str(repo), "status", "--porcelain"):
        raise RuntimeError("Commit the analysis source before submitting reproducible scans")
    globals_ = json.loads((base/"jobs.json").read_text())
    for size in (4096, 8192):
        job = globals_["jobs"][str(size)]["job_id"]
        states = run("sacct", "-X", "-n", "-P", "-j", job, "--format=JobID,State,ExitCode")
        if f"{job}|COMPLETED|0:0" not in states.splitlines():
            raise RuntimeError(f"Global job {job} has not completed successfully; leave it untouched")
        fit = base/f"n_{size}"/"fit"
        validation = json.loads((fit/"global_validation.json").read_text())
        best = json.loads((fit/"bestfit.json").read_text())
        if not validation["passed"] or not best["complete"] or best["fit_status"]:
            raise RuntimeError(f"Global fit validation missing or failed for {size}")
    report = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "commit": run("git", "-C", str(repo), "rev-parse", "HEAD"),
              "source_sha256": {name: hashlib.sha256((repo/name).read_bytes()).hexdigest()
                                for name in ("root_workspace_fit.py", "root_workspace_validation.py",
                                             "compare_asimov_hepdata.py", "submit_anchor_scans.py",
                                             "run_root_scan.slurm", "run_root_finalize.slurm",
                                             "run_anchor_comparison.slurm")},
              "global_jobs_completed_before_submission": True,
              "grid": [.01]+[i/20 for i in range(1, 31)], "jobs": []}

    def submit(kind, size, script, arguments, dependencies=(), hours="02:00:00"):
        output = base/f"n_{size}"/"fit" if size else base/"comparison"
        output.mkdir(exist_ok=True)
        command = ["sbatch", "--parsable", "--partition=cpu", "--account=pi_rclsa_umass_edu",
                   "--qos=normal", "--nodes=1", "--ntasks=1", "--cpus-per-task=1", "--mem=4G",
                   "--time="+hours, f"--job-name=anchor_{size}_{kind}",
                   f"--output={output}/{kind}_%j.log", "--kill-on-invalid-dep=yes"]
        if dependencies:
            command += ["--dependency=afterok:"+":".join(dependencies)]
        command += [str(repo/script), *map(str, arguments)]
        job = run(*command).split(";")[0]
        report["jobs"].append({"kind": kind, "size": size, "job_id": job, "command": command})
        manifest.write_text(json.dumps(report, indent=2)+"\n")
        print(f"{size} {kind}: {job}", flush=True)
        return job

    grid = report["grid"]
    chunks = [grid[:9], grid[8:17], grid[16:25], grid[24:]]
    finals = []
    for size in (4096, 8192):
        fit = base/f"n_{size}"/"fit"
        workspace = base/f"n_{size}"/"workspace"/f"HZZ_4l_asimov_{size}.root"
        scans = [submit(f"scan_{i}", size, "run_root_scan.slurm",
                        [args.setup, repo, workspace, fit/f"scan_{i}", fit/"bestfit.json", *points],
                        hours="03:00:00" if size == 8192 else "02:00:00")
                 for i, points in enumerate(chunks)]
        finals.append(submit("finalize", size, "run_root_finalize.slurm",
                             [args.setup, repo, workspace, fit, *grid], scans, "00:20:00"))
    submit("compare", None, "run_anchor_comparison.slurm",
           [args.setup, repo, base, args.reference.resolve()], finals, "00:15:00")


if __name__ == "__main__":
    main()
