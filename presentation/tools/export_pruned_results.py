#!/usr/bin/env python3
"""Copy the validated pruned scan and attributable per-mu resource results."""
import argparse
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--analysis-commit", required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    summary = json.loads((args.report / "summary.json").read_text())
    assert summary["complete"] and summary["resources"]["complete"]
    assert summary["points"] == 42 and summary["profiled_nuisances"] == 57
    assert all(summary["checks"]["protected_file_hashes_match"].values())
    summary["analysis_commit"] = args.analysis_commit
    encoded = json.dumps(summary, indent=2, allow_nan=False)
    (args.output / "data/pruned54_snapshot.json").write_text(encoded + "\n")
    (args.output / "data/pruned54_snapshot.js").write_text("window.PRUNED_RESULTS = " + encoded + ";\n")
    for name in ("profile_resources", "resources", "job_resources", "scan_validation", "bestfit_parameters", "kept_nuisances", "fixed_nuisances"):
        shutil.copyfile(args.report / (name + ".csv"), args.output / "data" / ("pruned54_" + name + ".csv"))
    for name in ("pruned_vs_hepdata", "profile_cost_vs_mu"):
        for extension in ("svg", "pdf"):
            shutil.copyfile(args.report / (name + "." + extension), args.output / "assets" / ("pruned54_" + name + "." + extension))
    print("Exported the complete 54-NP pruned study")


if __name__ == "__main__":
    main()
