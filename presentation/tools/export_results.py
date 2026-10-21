#!/usr/bin/env python3
"""Export a frozen, aggregate-only presentation snapshot from saved fit results.

Does not read ROOT workspaces or event arrays, run fits, or change inputs.
Requires matplotlib; run inside an allocation when generating the figures.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from compare_hepdata_scan import series_points

COLORS = {4096: "#126c96", 8192: "#c95639"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#b4bdc5", "axes.labelcolor": "#223447",
                     "text.color": "#223447", "xtick.color": "#526475",
                     "ytick.color": "#526475", "svg.fonttype": "none",
                     "svg.hashsalt": "offshell-results", "figure.facecolor": "white"})


def save_csv(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_scan(rows, reference, output, label, color, bounds, completed):
    fig, (ax, residual) = plt.subplots(2, 1, figsize=(10.5, 6.3), sharex=True,
                                      gridspec_kw={"height_ratios": [3.2, 1], "hspace": .08})
    points = sorted((x, y) for x, y in reference.items() if bounds[0] <= x <= bounds[1])
    ax.plot(*zip(*points), color="#263746", lw=1.9, label="HEPData · " + label)
    if rows:
        # Disconnected markers for unfinished scans: never imply that a gap was fitted.
        ax.plot([r["mu"] for r in rows], [r["q"] for r in rows],
                "o-" if completed else "o", color=color, ms=4.5, lw=1.5,
                label="Local · complete" if completed else "Local · converged points only")
    ax.set_ylabel(r"$q(\mu)=-2\ln\lambda(\mu)$")
    ax.set_ylim(bottom=-.025)
    ax.legend(loc="best", frameon=False, fontsize=10)
    ax.grid(alpha=.15)
    common = [r for r in rows if round(r["mu"], 10) in reference]
    residual.axhline(0, color="#82909c", lw=1)
    if bounds[0] >= .5:
        residual.axhspan(-.0005, .0005, color="#d9e2e8", alpha=.65)
    residual.plot([r["mu"] for r in common],
                  [r["q"]-reference[round(r["mu"], 10)] for r in common],
                  "o-" if completed else "o", color=color, ms=4)
    residual.set_ylabel("Local − HEPData", fontsize=10)
    residual.set_xlabel(r"Off-shell signal strength $\mu$")
    residual.set_xlim(*bounds)
    residual.grid(alpha=.15)
    residual.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    fig.subplots_adjust(left=.12, right=.98, bottom=.12, top=.96)
    fig.savefig(output, metadata={"Date": None})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=REPO/"presentation")
    args = parser.parse_args()
    out = args.output
    for sub in ("data", "assets"):
        (out/sub).mkdir(parents=True, exist_ok=True)
    sources = {}

    def read(relative):
        path = args.results/relative
        raw = path.read_bytes()
        sources[str(relative)] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    hep = read("hepdata_comparison/table_2.json")
    expected, observed = series_points(hep, "Exp NSBI"), series_points(hep, "Obs NSBI")
    data_fit = read("data_workspace/bestfit.json")
    data_rows = read("data_workspace/scan_points.json")
    observed_comparison = read("hepdata_comparison/comparison_summary.json")
    construction = read("asimov_workspaces/summary.json")
    quadrature = read("asimov_fits/quadrature_diagnostic.json")
    source = read("asimov_fits/source_manifest.json")
    snapshot = {"created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "analysis_commit": source["commit"], "root_version": data_fit["root_version"],
                "reference_doi": hep["doi"],
                "observed": {k: data_fit[k] for k in ("mu_hat", "edm", "covariance_quality")},
                "samples": {}, "quadrature_first_start": quadrature}
    snapshot["observed"].update({k: observed_comparison[k] for k in
                                 ("max_abs_delta_q", "rms_delta_q", "comparison_points")})
    plot_scan(data_rows, observed, out/"assets/observed.svg", "observed NSBI", "#126c96", [.5, 1.5], True)
    save_csv(out/"data/observed_scan.csv", data_rows, ["mu", "q", "nll", "edm", "migrad_status"])
    save_csv(out/"data/observed_parameters.csv", data_fit["parameters"], ["name", "value", "error", "min", "max"])
    save_csv(out/"data/hepdata.csv", [{"mu": mu, "observed": observed[mu], "expected": expected[mu]}
                                   for mu in sorted(expected)], ["mu", "observed", "expected"])
    for size in (4096, 8192):
        prefix = Path(f"asimov_fits/n_{size}")
        best_name = "bestfit.json" if (args.results/prefix/"bestfit.json").exists() else "bestfit_checkpoint.json"
        best = read(prefix/best_name)
        by_mu, overlap = {}, []
        final_scan = args.results/prefix/"scan_points.json"
        paths = [final_scan] if final_scan.exists() else sorted((args.results/prefix).glob("scan_*/scan_points.json"))
        for path in paths:
            for row in read(path.relative_to(args.results)):
                if row["migrad_status"] or row["edm"] > 1e-4 or not math.isfinite(row["nll"]):
                    continue
                mu = round(row["mu"], 10)
                if mu in by_mu:
                    overlap.append(2*abs(row["nll"]-by_mu[mu]["nll"]))
                if mu not in by_mu or row["nll"] < by_mu[mu]["nll"]:
                    by_mu[mu] = row
        rows = []
        for mu, row in sorted(by_mu.items()):
            q = 2*(row["nll"]-best["nll_min"])
            if q < -2e-4:
                raise ValueError("Profile improves the global minimum; resolve before publishing")
            rows.append({"mu": mu, "q": q, "nll": row["nll"], "edm": row["edm"],
                         "migrad_status": row["migrad_status"], "reference_q": expected.get(mu),
                         "delta_q": q-expected[mu] if mu in expected else None})
        validation_file = args.results/prefix/"likelihood_validation.json"
        validated = bool(validation_file.exists() and read(prefix/validation_file.name)["passed"])
        complete = best["complete"] and len(rows) == 31 and validated and max(overlap, default=0) <= 1e-4
        record = next(r for r in construction["workspaces"] if r["reference_events"] == size)
        params = best["parameters"]
        attempts = best["attempts"]
        sample = {"mu_hat": best["mu_hat"], "nll_min": best["nll_min"],
                  "complete": complete, "fit_complete": best["complete"], "validated": validated,
                  "points": len(rows), "rows": rows,
                  "covariance_quality": best.get("covariance_quality"),
                  "edm": best.get("edm", min(a["edm"] for a in attempts)),
                  "start_mu_estimates": [a["mu_hat"] for a in attempts],
                  "start_nll_spread": max(a["nll"] for a in attempts)-min(a["nll"] for a in attempts),
                  "parameters": [{k: p[k] for k in ("name", "value", "min", "max")} for p in params],
                  "construction": {k: record[k] for k in ("dataset_rows", "sum_weights", "negative_weights", "max_delta_nll_difference")}}
        shared = [r for r in rows if r["delta_q"] is not None]
        sample["max_abs_delta_q_available"] = max((abs(r["delta_q"]) for r in shared), default=None)
        sample["shared_points"] = len(shared)
        snapshot["samples"][str(size)] = sample
        plot_scan(rows, expected, out/f"assets/asimov_{size}.svg", "expected NSBI", COLORS[size], [0, 1.5], complete)
        save_csv(out/f"data/asimov_{size}_scan.csv", rows,
                 ["mu", "q", "reference_q", "delta_q", "nll", "edm", "migrad_status"])
        save_csv(out/f"data/asimov_{size}_parameters.csv", params, ["name", "value", "error", "min", "max"])
    snapshot["sources_sha256"] = sources
    snapshot["status"] = "complete" if all(s["complete"] for s in snapshot["samples"].values()) else "preliminary"
    (out/"data/snapshot.json").write_text(json.dumps(snapshot, indent=2, allow_nan=False)+"\n")
    (out/"data/snapshot.js").write_text("window.RESULTS = "+json.dumps(snapshot, ensure_ascii=True, allow_nan=False)+";\n")
    print(json.dumps({"created_utc": snapshot["created_utc"], "status": snapshot["status"],
                      "points": {n:s["points"] for n,s in snapshot["samples"].items()}}, indent=2))


if __name__ == "__main__":
    main()
