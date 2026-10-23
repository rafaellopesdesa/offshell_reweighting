#!/usr/bin/env python3
"""Compare each current 61-point scan with unchanged HEPData Table 2 Exp NSBI."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil

from compare_hepdata_scan import series_points, csv_series_points


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot(model, reference, output, limits):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
        "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none", "svg.hashsalt": "offshell-family-61"})
    rows = model["rows"]
    color = model["color"]
    fig, (ax, residual) = plt.subplots(2, 1, figsize=(10.5, 6.3), sharex=True,
                                      gridspec_kw={"height_ratios": [3.1, 1], "hspace": .08})
    ax.plot(list(reference), list(reference.values()), "--", color="#253746", lw=1.8,
            label="HEPData expected NSBI (4ℓ, stat+syst)")
    ax.plot([r["mu"] for r in rows], [r["q"] for r in rows], "o-", color=color,
            lw=1.8, ms=3.6, label=model["title"])
    ax.plot(model["mu_hat"], 0., "*", color=color, ms=12, zorder=5)
    ax.set_ylabel(r"$q(\mu)=-2\ln\Lambda(\mu)$")
    ax.set_ylim(-.15, limits["q_max"])
    ax.set_title(f"{model['reference_events']:,} reference events · generating μ = 1", fontsize=13)
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    residual.plot([r["mu"] for r in rows], [r["delta_q"] for r in rows], "o-",
                  color=color, lw=1.5, ms=3.3)
    residual.axhline(0, color="#82909c", lw=1)
    residual.set_ylabel("Scan − HEPData", fontsize=10)
    residual.set_xlabel(r"Off-shell signal strength $\mu$")
    residual.set_ylim(-limits["residual_abs"], limits["residual_abs"])
    residual.set_xlim(0, 3)
    residual.set_xticks([i / 2 for i in range(7)])
    for axis in (ax, residual):
        axis.grid(alpha=.17)
    fig.subplots_adjust(left=.12, right=.985, bottom=.14, top=.925)
    fig.text(.12, .016, "61 exact shared coordinates · HEPData 156805 v1, Table 2 · star: continuous fitted minimum",
             fontsize=8.5, color="#536777")
    for extension in ("svg", "pdf", "png"):
        fig.savefig(output / f"{model['label']}_vs_hepdata.{extension}", dpi=180,
                    metadata={"Date": None} if extension == "svg" else None)
    plt.close(fig)


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run report generation inside a Slurm allocation")
    args.output.mkdir(parents=True, exist_ok=True)
    configuration = json.loads((args.run / "configuration.json").read_text())
    table = json.loads(args.hepdata.read_text())
    expected = series_points(table, "Exp NSBI")
    assert expected == csv_series_points(args.hepdata.with_suffix(".csv"), "Exp NSBI")
    grid = [i / 20 for i in range(61)]
    assert configuration["grid"] == grid and sorted(expected) == grid
    models, missing = [], []
    for spec in configuration["workspaces"]:
        label = spec["label"]
        directory = args.run / label
        path = args.run / "csv" / f"{label}.csv"
        check_path = directory / "final/likelihood_validation.json"
        if not path.exists() or not check_path.exists():
            missing.append(label)
            continue
        check = json.loads(check_path.read_text())
        best = json.loads((directory / "fit/bestfit.json").read_text())
        merge = json.loads((directory / "final/scan_summary.json").read_text())
        assert check["passed"] and best["complete"] and merge["all_converged"]
        assert best["input_sha256"] == merge["input_sha256"] == spec["sha256"]
        assert len(check["scan_checks"]) == 61
        with path.open() as f:
            scan = list(csv.DictReader(f))
        assert [float(r["mu"]) for r in scan] == grid
        rows = []
        for r in scan:
            mu, q = float(r["mu"]), float(r["q"])
            assert r["input_sha256"] == spec["sha256"]
            assert int(r["status"]) == 0 and math.isfinite(float(r["edm"])) and float(r["edm"]) <= 1e-4
            assert math.isfinite(q) and q >= -2e-4
            assert abs(q - 2 * (float(r["profile_nll"]) - best["nll_min"])) < 1e-9
            assert abs(float(r["mu_hat"]) - best["mu_hat"]) < 1e-12
            rows.append({"mu": mu, "q": q, "hepdata_q": expected[mu], "delta_q": q - expected[mu]})
        title = "True Asimov" if label.startswith("true_") else "Data-like Asimov"
        if spec["constrained_nps"] == 54:
            title += " · pruned to 54 NPs"
        color = "#126c96" if label.startswith("true_") else "#7d4698" if spec["constrained_nps"] == 54 else "#c95639"
        extreme = max(rows, key=lambda r: abs(r["delta_q"]))
        model = {"label": label, "title": title, "color": color,
            "reference_events": spec["reference_events"], "constrained_nps": spec["constrained_nps"],
            "free_qqzz": 3, "points": 61, "mu_hat": best["mu_hat"],
            "q_zero": rows[0]["q"], "q_three": rows[-1]["q"],
            "rms_delta_q": math.sqrt(sum(r["delta_q"] ** 2 for r in rows) / 61),
            "max_abs_delta_q": abs(extreme["delta_q"]), "mu_at_max_abs_delta_q": extreme["mu"],
            "max_profile_edm": max(float(r["edm"]) for r in scan),
            "max_independent_nll_error": max(abs(float(r["independent_delta_nll_error"])) for r in scan),
            "minimum_fitted_intensity": min(float(r["minimum_intensity"]) for r in scan),
            "scan_csv_sha256": digest(path), "workspace_sha256": spec["sha256"],
            "validation_sha256": digest(check_path), "global_fit_sha256": digest(directory / "fit/bestfit.json"),
            "rows": rows}
        models.append(model)
        shutil.copyfile(path, args.output / f"{label}.csv")
        write_csv(args.output / f"{label}_comparison.csv", rows)
    if missing and not args.preview:
        raise RuntimeError("Wait for validated scans: " + ", ".join(missing))
    assert models
    limits = {"q_max": math.ceil(2 * 1.1 * max(max(expected.values()), *(r["q"] for m in models for r in m["rows"]))) / 2,
              "residual_abs": math.ceil(10 * 1.1 * max(m["max_abs_delta_q"] for m in models)) / 10}
    for model in models:
        plot(model, expected, args.output, limits)
    write_csv(args.output / "hepdata_expected_0_3.csv", [{"mu": mu, "q": expected[mu]} for mu in grid])
    write_csv(args.output / "comparison_summary.csv", [{k: v for k, v in m.items() if k not in ("rows", "color")} for m in models])
    summary = {"complete": not missing, "missing": missing, "models": models,
        "created_utc": datetime.now(timezone.utc).isoformat(), "grid": grid, "limits": limits,
        "reference_doi": table["doi"], "reference_series": "Exp NSBI", "reference_points": 61,
        "reference_json_sha256": digest(args.hepdata), "reference_csv_sha256": digest(args.hepdata.with_suffix(".csv")),
        "comparison_definition": "q_scan(mu) - q_HEPData(mu) at 61 exact shared coordinates; no interpolation, re-zeroing, rescaling, or clipping",
        "normalization_implementation": "RooSameSampleNormalizer_v1 for all six workspaces",
        "scan_source_commit": json.loads((args.run / "publication.json").read_text())["commit"]}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"complete": summary["complete"], "models": len(models), "missing": missing, "limits": limits}, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--hepdata", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--preview", action="store_true", help="Render available validated models for private layout checks only")
    main(p.parse_args())
