#!/usr/bin/env python3
"""Publish six validated CSV/HEPData comparisons in a consecutive slide block."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1])
    args = p.parse_args()
    data = json.loads((args.report / "summary.json").read_text())
    completion = json.loads((args.run / "completion.json").read_text())
    assert data["complete"] and len(data["models"]) == 6
    assert completion["complete"] and completion["total_points"] == 366
    assert [m["reference_events"] for m in data["models"]] == [4096] * 3 + [8192] * 3
    assert [m["constrained_nps"] for m in data["models"]] == [127, 127, 54] * 2
    for model in data["models"]:
        assert model["points"] == 61 and len(model["rows"]) == 61
        label = model["label"]
        assert digest(args.report / f"{label}.csv") == model["scan_csv_sha256"]
        for suffix in (".csv", "_comparison.csv"):
            shutil.copyfile(args.report / (label + suffix), args.output / "data" / ("family_" + label + suffix))
        for extension in ("svg", "pdf"):
            shutil.copyfile(args.report / f"{label}_vs_hepdata.{extension}",
                            args.output / "assets" / f"family_{label}_vs_hepdata.{extension}")
    for name in ("comparison_summary.csv", "hepdata_expected_0_3.csv"):
        shutil.copyfile(args.report / name, args.output / "data" / ("family_" + name))
    archive = args.run / "likelihood_scans_csv.zip"
    assert digest(archive) == completion["archive_sha256"]
    shutil.copyfile(archive, args.output / "data/family_scans_csv.zip")
    encoded = json.dumps(data, indent=2, allow_nan=False)
    (args.output / "data/family_snapshot.json").write_text(encoded + "\n")
    (args.output / "data/family_snapshot.js").write_text("window.FAMILY_RESULTS = " + encoded + ";\n")
    print("Exported all six complete 61-point comparisons and the CSV bundle")


if __name__ == "__main__":
    main()
