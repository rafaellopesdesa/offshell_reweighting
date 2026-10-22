#!/usr/bin/env python3
"""Collect the two completed anchor-normalized global fits into JSON and CSV."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    rows, sources = [], {}
    for size in (4096, 8192):
        sample = args.results / f"n_{size}"
        inputs = {"fit": sample / "fit/bestfit.json",
                  "global_validation": sample / "fit/global_validation.json",
                  "construction_validation": sample / "workspace/validation.json",
                  "construction": sample / "workspace/construction.json"}
        records = {name: json.loads(path.read_text()) for name, path in inputs.items()}
        fit, validation = records["fit"], records["global_validation"]
        construction_check = records["construction_validation"]
        assert fit["complete"] and validation["passed"] and construction_check["passed"]
        assert fit["input_sha256"] == validation["input_sha256"] == records["construction"]["workspace_sha256"]
        assert fit["fit_status"] == 0 and fit["edm"] < 1e-6
        alphas = [p["value"] for p in fit["parameters"] if p["name"].startswith("alpha_")]
        starts = validation["start_mu_estimates"]
        rows.append({"reference_events": size, "mu_hat": fit["mu_hat"],
                     "mu_hat_minus_one": fit["mu_hat"]-1,
                     "percent_above_one": 100*(fit["mu_hat"]-1),
                     "fit_status": fit["fit_status"], "edm": fit["edm"],
                     "maximum_absolute_alpha": max(map(abs, alphas)),
                     "start_mu_spread": max(starts)-min(starts),
                     "start_nll_spread": validation["start_nll_spread"],
                     "maximum_global_delta_nll_disagreement": max(abs(x["delta_nll_disagreement"]) for x in validation["checks"]),
                     "maximum_anchor_relative_closure_error": construction_check["maximum_anchor_relative_closure_error"],
                     "sr_integral_relative_difference_at_minimum": validation["checks"][-1]["sr_integral_relative_difference"],
                     "largest_nominal_systematic_score": validation["largest_nominal_systematic_score"]})
        sources[str(size)] = {name: {"path": str(path.relative_to(args.results)), "sha256": digest(path)}
                              for name, path in inputs.items()}
    report = {"complete": True, "rows": rows, "sources": sources,
              "source_commit": json.loads((args.results / "jobs.json").read_text())["commit"],
              "report_source_sha256": digest(Path(__file__)),
              "scope": "Unconditional fits from two starts with all 131 parameters free; no scans or explicit full Hesse. Anchor closure does not enforce closure between anchors."}
    (args.results / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    with (args.results / "global_fit_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
