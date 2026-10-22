#!/usr/bin/env python3
"""Export validated aggregate anchor-normalized scan results for the repository."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    sources = {}

    def read(relative):
        raw = (args.results/relative).read_bytes()
        sources[str(relative)] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    def csv_write(name, rows, fields):
        with (args.output/name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    jobs = read("scan_jobs.json")
    globals_ = read("summary.json")
    comparison = read("comparison/comparison_summary.json")
    resources = read("resources.json")
    report = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "export_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "global_fit_commit": globals_["source_commit"], "scan_commit": jobs["commit"],
              "reference_doi": comparison["reference_doi"], "reference_series": comparison["reference_series"],
              "comparison_method": comparison["comparison_method"],
              "sample_difference_4096_minus_8192": comparison["sample_difference_4096_minus_8192"],
              "resource_phases": resources["phases"], "resource_definitions": resources["definitions"],
              "evaluated_grid": [.01]+[i/20 for i in range(1,31)],
              "model_caveat": comparison["model_caveat"], "samples": {}}
    args.output.mkdir(parents=True, exist_ok=True)
    for size in (4096,8192):
        prefix = Path(f"n_{size}/fit")
        best = read(prefix/"bestfit.json")
        scan = read(prefix/"scan_points.json")
        summary = read(prefix/"scan_summary.json")
        validation = read(prefix/"likelihood_validation.json")
        global_validation = read(prefix/"global_validation.json")
        construction = read(prefix.parent/"workspace/construction.json")
        if not validation["passed"] or not global_validation["passed"] or not summary["all_converged"]:
            raise ValueError("Require completed independent validation")
        if best["input_sha256"] != construction["workspace_sha256"] or best["input_sha256"] != summary["input_sha256"]:
            raise ValueError("Workspace checksums disagree")
        if [r["mu"] for r in scan] != report["evaluated_grid"] or len(validation["scan_checks"]) != 31:
            raise ValueError("Incomplete profile grid")
        for row in scan:
            if row["migrad_status"] or row["edm"] > 1e-4 or row["q"] < -2e-4 or abs(row["q"]-2*(row["nll"]-best["nll_min"])) > 1e-8:
                raise ValueError("Unconverged profile or inconsistent baseline")
        diagnostics = next(r for r in globals_["rows"] if r["reference_events"] == size)
        sample = {**diagnostics, "workspace_sha256": best["input_sha256"],
                  "nll_min": best["nll_min"], "start_mu_estimates": global_validation["start_mu_estimates"],
                  "covariance_quality": best["covariance_quality"],
                  "full_hesse_performed": global_validation["full_hesse_performed"],
                  "profile_points": len(scan), "maximum_profile_edm": max(r["edm"] for r in scan),
                  "maximum_overlap_delta_q": max(r["absolute_q_difference"] for r in summary["overlap_checks"]),
                  "maximum_profile_delta_nll_disagreement": max(abs(r["delta_nll_disagreement"]) for r in validation["scan_checks"]),
                  "minimum_profile_intensity": min(r["minimum_current_intensity"] for r in validation["scan_checks"]),
                  "maximum_profile_absolute_sr_integral_mismatch": max(abs(r["sr_integral_relative_difference"]) for r in validation["scan_checks"]),
                  "comparison": comparison["results"][str(size)]["comparison"],
                  "q_at_0p01": scan[0]["q"]}
        report["samples"][str(size)] = sample
        csv_write(f"scan_{size}.csv", scan, ["mu", "q", "nll", "migrad_status", "edm", "seconds"])
        csv_write(f"parameters_{size}.csv", best["parameters"], ["name", "value", "error", "min", "max"])
        csv_write(f"validation_{size}.csv", validation["scan_checks"], list(validation["scan_checks"][0]))
    csv_write("global_fit_comparison.csv", globals_["rows"], list(globals_["rows"][0]))
    csv_write("resources.csv", resources["phases"], list(resources["phases"][0]))
    accounting = (args.results/"slurm_accounting.txt").read_bytes()
    sources["slurm_accounting.txt"] = hashlib.sha256(accounting).hexdigest()
    if sources["slurm_accounting.txt"] != resources["accounting_sha256"]:
        raise ValueError("Slurm accounting checksum disagrees with the resource report")
    report["sources_sha256"] = sources
    (args.output/"summary.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps(report["samples"], indent=2))


if __name__ == "__main__":
    main()
