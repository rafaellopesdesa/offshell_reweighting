#!/usr/bin/env python3
"""Run and export the six preserved workspace scans on a common 61-point grid."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def write_csv(path, fields, rows):
    tmp = path.with_suffix(".csv.tmp")
    with tmp.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def execute(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Use Slurm for fits and validation")
    run = args.run.resolve()
    config = read(run / "configuration.json")
    source = run / "source"
    for name, expected in config["source_sha256"].items():
        if digest(source / name) != expected:
            raise RuntimeError(f"Frozen source changed: {name}")
    if args.stage == "export":
        export_all(run, config)
        return
    spec = config["workspaces"][args.model_index]
    input_path = Path(spec["source_path"])
    if digest(input_path) != spec["sha256"]:
        raise RuntimeError("Workspace checksum changed")
    directory = run / spec["label"]
    best = directory / "fit/bestfit.json"
    base = ["--input", str(input_path), "--weighted-asimov"]

    def call(script, *arguments):
        subprocess.run([sys.executable, "-u", str(source / script), *map(str, arguments)], check=True)

    if args.stage == "global":
        call("root_workspace_fit.py", "fit", *base, "--output", directory / "fit",
             "--starts", ".6", "1.4", "--fit-strategy", "0", "--fit-eps", ".0001")
        call("root_workspace_validation.py", *base, "--bestfit", best,
             "--range", "0", "3", "--allow-approximate-covariance",
             "--output", directory / "fit/likelihood_validation.json")
    elif args.stage == "scan":
        if not read(directory / "fit/likelihood_validation.json")["passed"]:
            raise RuntimeError("Global likelihood validation is incomplete")
        chunk = int(os.environ["SLURM_ARRAY_TASK_ID"])
        grid = config["chunks"][chunk]
        call("root_workspace_fit.py", "scan", *base, "--bestfit", best,
             "--output", directory / f"scan/chunk_{chunk}", "--range", "0", "3",
             "--allow-zero", "--grid", *grid)
    elif args.stage == "finalize":
        parts = [directory / f"scan/chunk_{i}/scan_points.json" for i in range(len(config["chunks"]))]
        final = directory / "final"
        call("root_workspace_fit.py", "finalize", *base, "--bestfit", best,
             "--output", final, "--range", "0", "3", "--points", "61",
             "--allow-zero", "--scan-parts", *parts)
        call("root_workspace_validation.py", *base, "--bestfit", best,
             "--range", "0", "3", "--allow-approximate-covariance",
             "--scan-points", final / "scan_points.json",
             "--output", final / "likelihood_validation.json")
        export_one(run, spec)
    else:
        raise ValueError(args.stage)
    if digest(input_path) != spec["sha256"]:
        raise RuntimeError("Input changed during calculation")


def export_one(run, spec):
    directory = run / spec["label"]
    best = read(directory / "fit/bestfit.json")
    points = read(directory / "final/scan_points.json")
    check = read(directory / "final/likelihood_validation.json")
    assert check["passed"] and len(points) == 61
    assert [round(r["mu"], 10) for r in points] == [i / 20 for i in range(61)]
    checks = {round(r["mu"], 10): r for r in check["scan_checks"]}
    assert best["input_sha256"] == spec["sha256"] and len(checks) == 61
    records = []
    for point in points:
        if point["migrad_status"] or not math.isfinite(point["edm"]) or point["edm"] > 1e-4:
            raise RuntimeError("Unconverged profile cannot be exported")
        tested = checks[round(point["mu"], 10)]
        records.append({"mu": point["mu"], "q": point["q"], "profile_nll": point["nll"],
            "global_nll": best["nll_min"], "mu_hat": best["mu_hat"],
            "status": point["migrad_status"], "edm": point["edm"],
            "wall_seconds": point["seconds"], "cpu_seconds": point["cpu_seconds"],
            "process_peak_rss_gib": point["process_peak_rss_gib"],
            "likelihood_evaluations": point["likelihood_evaluations"],
            "minimum_intensity": tested["minimum_current_intensity"],
            "independent_delta_nll_error": tested["delta_nll_disagreement"],
            "sr_integral_relative_error": tested["sr_integral_relative_difference"],
            "profiled_parameters": spec["constrained_nps"] + 3,
            "input_sha256": spec["sha256"]})
    destination = run / "csv"
    write_csv(destination / f"{spec['label']}.csv", list(records[0]), records)
    parameter_rows = [{"mu": p["mu"], "parameter": name, "value": value}
                      for p in points for name, value in sorted(p["parameters"].items())]
    write_csv(destination / f"{spec['label']}_profiled_parameters.csv",
              ["mu", "parameter", "value"], parameter_rows)
    shutil.copyfile(directory / "fit/bestfit_parameters.csv",
                    destination / f"{spec['label']}_bestfit_parameters.csv")
    save(directory / "export_complete.json", {"passed": True, "points": 61,
         "mu_hat": best["mu_hat"], "input_sha256": spec["sha256"],
         "created_utc": datetime.now(timezone.utc).isoformat()})


def export_all(run, config):
    rows, baselines = [], []
    for spec in config["workspaces"]:
        assert read(run / spec["label"] / "export_complete.json")["passed"]
        assert digest(Path(spec["source_path"])) == spec["sha256"]
        with (run / "csv" / f"{spec['label']}.csv").open() as f:
            for row in csv.DictReader(f):
                rows.append({"workspace": spec["label"], **row})
        best = read(run / spec["label"] / "fit/bestfit.json")
        baselines.append({"workspace": spec["label"], "mu_hat": best["mu_hat"],
            "nll_min": best["nll_min"], "status": best["fit_status"], "edm": best["edm"],
            "covariance_quality": best["covariance_quality"], "input_sha256": spec["sha256"]})
    assert len(rows) == 366
    write_csv(run / "csv/all_scans.csv", list(rows[0]), rows)
    write_csv(run / "csv/global_fits.csv", list(baselines[0]), baselines)
    archive = run / "likelihood_scans_csv.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
        for path in sorted((run / "csv").glob("*.csv")):
            z.write(path, arcname=path.name)
        z.write(run / "README.md", arcname="README.md")
    save(run / "completion.json", {"complete": True, "workspaces": 6, "points_per_workspace": 61,
         "total_points": 366, "all_independently_validated": True,
         "original_input_hashes_unchanged": True, "archive_sha256": digest(archive),
         "created_utc": datetime.now(timezone.utc).isoformat()})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("global", "scan", "finalize", "export"))
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--model-index", type=int, default=0)
    execute(p.parse_args())
