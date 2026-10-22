#!/usr/bin/env python3
"""Report the completed 54-NP data-like scan and measured per-coordinate cost."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess

from collect_anchor_resources import seconds, memory_gib
from compare_hepdata_scan import series_points, csv_series_points


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def accounting(run, manifest, grid):
    selected = [j for j in manifest["jobs"] if j["kind"] in
                ("build_validation", "global_fit", "scan_array", "finalize", "superseded_build_validation")]
    fields = ["JobID", "JobIDRaw", "JobName", "State", "ElapsedRaw", "TotalCPU", "MaxRSS",
              "ReqMem", "AllocCPUS", "Start", "End", "NodeList", "ExitCode"]
    command = ["sacct", "-n", "-P", "--array", "--units=K", "-j",
               ",".join(j["job_id"] for j in selected), "--format=" + ",".join(fields)]
    raw = subprocess.check_output(command, text=True, env={**os.environ, "TZ": "UTC"})
    (run / "slurm_accounting.txt").write_text("|".join(fields) + "\n" + raw)
    records = [dict(zip(fields, line.split("|"))) for line in raw.splitlines() if line.strip()]
    lookup = {r["JobID"]: r for r in records}
    jobs = []
    for spec in selected:
        tasks = [(f"{spec['job_id']}_{i}", mu) for i, mu in enumerate(grid)] if spec["kind"] == "scan_array" else [(spec["job_id"], None)]
        for job_id, mu in tasks:
            row = lookup[job_id]
            superseded = spec["kind"] == "superseded_build_validation"
            if not superseded and (row["State"] != "COMPLETED" or row["ExitCode"] != "0:0"):
                raise RuntimeError(f"Required job not successfully complete: {job_id} {row['State']}")
            if superseded and not row["State"].startswith("CANCELLED"):
                raise RuntimeError("Unexpected superseded-check status")
            peaks = [memory_gib(r["MaxRSS"]) for r in records if r["JobID"].startswith(job_id + ".") and r["MaxRSS"]]
            if not peaks or seconds(row["TotalCPU"]) <= 0:
                raise RuntimeError("Slurm accounting is not yet final: " + job_id)
            jobs.append({"job_id": job_id, "phase": spec["kind"], "mu": mu, "state": row["State"],
                "elapsed_seconds": int(row["ElapsedRaw"]), "cpu_seconds": seconds(row["TotalCPU"]),
                "peak_rss_gib": max(peaks), "requested_memory": row["ReqMem"], "cpus": int(row["AllocCPUS"]),
                "start_utc": row["Start"], "end_utc": row["End"], "node": row["NodeList"]})
    phases = []
    for phase in ("superseded_build_validation", "build_validation", "global_fit", "scan_array", "finalize"):
        group = [j for j in jobs if j["phase"] == phase]
        if not group:
            continue
        phases.append({"phase": phase, "jobs": len(group),
            "wall_seconds": (max(datetime.fromisoformat(j["end_utc"]) for j in group)
                             - min(datetime.fromisoformat(j["start_utc"]) for j in group)).total_seconds(),
            "cpu_seconds": sum(j["cpu_seconds"] for j in group),
            "largest_job_peak_rss_gib": max(j["peak_rss_gib"] for j in group)})
    return {"complete": True, "jobs": jobs, "phases": phases,
        "accounting_sha256": digest(run / "slurm_accounting.txt"),
        "definitions": {
            "fit_seconds": "Time inside Migrad and result capture for this coordinate, including strategy retry if needed; excludes load and serialization",
            "fit_cpu_seconds": "Process CPU time over the same measured minimization interval",
            "process_peak_rss_gib": "Linux getrusage process high-water RSS through completion of this fit, including workspace load; one process per mu",
            "slurm_peak_rss_gib": "Largest step MaxRSS for that one-coordinate job, sampled by Slurm",
            "phase_wall_seconds": "First job start to last job finish, including gaps within the phase",
            "phase_cpu_seconds": "Summed Slurm parent TotalCPU, without double-counting steps",
            "phase_peak_rss": "Largest individual job RSS, not summed concurrent memory",
            "scope": "Includes the canceled checksum-helper prototype separately. Pending canceled jobs never ran. Aggregate report generation is excluded.",
            "comparability": "Independent global-fit starts at each mu; older unpruned scans used warm-start chunks and different nodes. Runtime differences are not a controlled pruning benchmark."}}


def plots(rows, expected, summary, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none", "svg.hashsalt": "offshell-pruned54"})
    color = "#7d4698"
    fig, (ax, residual) = plt.subplots(2, 1, figsize=(10.5, 6.3), sharex=True,
                                     gridspec_kw={"height_ratios": [3.2, 1], "hspace": .08})
    ref = sorted((x, y) for x, y in expected.items() if 0 <= x <= 2)
    ax.plot(*zip(*ref), color="#253746", lw=1.8, label="HEPData expected NSBI (4ℓ, stat+syst)")
    ax.plot([r["mu"] for r in rows], [r["q"] for r in rows], "o-", color=color, ms=3.8,
            label="Data-like Asimov · 54 constrained NPs + 3 free qqZZ scales")
    ax.plot(summary["mu_hat"], 0, "*", color=color, ms=12)
    ax.set_ylabel(r"$q(\mu)=-2\ln\lambda(\mu)$")
    ax.set_ylim(bottom=-.04)
    ax.set_title("4,096 fixed reference events · pruned data-like model · generating μ = 1", fontsize=12.5)
    ax.legend(frameon=False, fontsize=9.5, loc="upper center")
    shared = [r for r in rows if r["hepdata_q"] is not None]
    residual.plot([r["mu"] for r in shared], [r["delta_q"] for r in shared], "o-", color=color, ms=3.4)
    residual.axhline(0, color="#82909c", lw=1)
    residual.set_ylabel("Pruned − HEPData", fontsize=9.5)
    residual.set_xlabel(r"Off-shell signal strength $\mu$")
    residual.set_xlim(0, 2)
    for axis in (ax, residual): axis.grid(alpha=.16)
    fig.subplots_adjust(left=.13, right=.98, bottom=.12, top=.92)
    save_figure(fig, output / "pruned_vs_hepdata")
    plt.close(fig)
    fig, (timing, memory) = plt.subplots(2, 1, figsize=(10.5, 6.3), sharex=True,
                                      gridspec_kw={"hspace": .15})
    mu = [r["mu"] for r in rows]
    timing.plot(mu, [r["fit_seconds"] for r in rows], "o-", ms=3.5, color=color, label="Minimizer wall time")
    timing.plot(mu, [r["fit_cpu_seconds"] for r in rows], "--", color="#d18b34", label="Minimizer CPU time")
    timing.set_ylabel("Seconds")
    timing.set_ylim(bottom=0)
    timing.legend(frameon=False, fontsize=10)
    timing.set_title("Pruned 4096 scan · one CPU and one independent process per μ", fontsize=13)
    memory.plot(mu, [r["process_peak_rss_gib"] for r in rows], "o-", color=color, ms=3.5)
    memory.set_ylabel("Process peak RSS [GiB]\n(workspace load + fit)", fontsize=10)
    memory.set_xlabel(r"Off-shell signal strength $\mu$")
    memory.set_xlim(0, 2); memory.set_ylim(bottom=0)
    for axis in (timing, memory): axis.grid(alpha=.16)
    fig.subplots_adjust(left=.13, right=.98, bottom=.12, top=.92)
    save_figure(fig, output / "profile_cost_vs_mu")
    plt.close(fig)


def save_figure(figure, stem):
    for extension in ("svg", "pdf", "png"):
        figure.savefig(stem.with_suffix("." + extension), dpi=180,
                       metadata={"Date": None} if extension == "svg" else None)


def main(args):
    run = args.run; output = args.output
    output.mkdir(parents=True, exist_ok=True)
    inputs = {}
    def read(path):
        inputs[str(path)] = digest(path)
        return json.loads(path.read_text())
    manifest = read(run / "jobs.json")
    construction = read(run / "workspace/construction.json")
    build_check = read(run / "workspace/validation.json")
    global_check = read(run / "fit/global_validation.json")
    best = read(run / "fit/bestfit.json")
    scan = read(run / "final/scan_points.json")
    scan_summary = read(run / "final/scan_summary.json")
    validation = read(run / "final/likelihood_validation.json")
    grid = read(run / "grid.json")
    table = read(args.results / "hepdata_comparison/table_2.json")
    expected = series_points(table, "Exp NSBI")
    assert expected == csv_series_points(args.results / "hepdata_comparison/table_2.csv", "Exp NSBI")
    assert build_check["passed"] and global_check["passed"] and validation["passed"] and best["complete"]
    assert scan_summary["all_converged"] and scan_summary["profiled"] == 57
    assert [r["mu"] for r in scan] == grid and len(grid) == 42
    assert len(validation["scan_checks"]) == 42 and len(best["parameters"]) == 58
    assert best["input_sha256"] == construction["workspace_sha256"] == scan_summary["input_sha256"]
    resources = accounting(run, manifest, grid)
    per_mu = {j["mu"]: j for j in resources["jobs"] if j["phase"] == "scan_array"}
    unpruned = {r["mu"]: r["q"] for r in read(args.results / "asimov_data_like_4096/data_like/final/scan_points.json")}
    rows = []
    for r in scan:
        assert r["migrad_status"] == 0 and r["edm"] <= 1e-4
        assert abs(r["q"] - 2*(r["nll"]-best["nll_min"])) < 1e-9
        assert r["q"] >= -2e-4
        # An exact restricted model cannot improve on the full profile,
        # modulo the already measured numerical minimization tolerance.
        assert r["q"] - unpruned[r["mu"]] >= -2e-4
        j = per_mu[r["mu"]]; ref = expected.get(r["mu"])
        rows.append({"mu": r["mu"], "q": r["q"], "hepdata_q": ref,
            "delta_q": r["q"]-ref if ref is not None else None,
            "pruned_minus_unpruned_q": r["q"]-unpruned[r["mu"]],
            "nll": r["nll"], "migrad_status": r["migrad_status"], "edm": r["edm"],
            "fit_seconds": r["seconds"], "fit_cpu_seconds": r["cpu_seconds"],
            "likelihood_evaluations": r["likelihood_evaluations"], "process_peak_rss_gib": r["process_peak_rss_gib"],
            "slurm_elapsed_seconds": j["elapsed_seconds"], "slurm_cpu_seconds": j["cpu_seconds"],
            "slurm_peak_rss_gib": j["peak_rss_gib"], "job_id": j["job_id"], "node": j["node"],
            "requested_memory": j["requested_memory"], "cpus": j["cpus"]})
    shared = [r for r in rows if r["delta_q"] is not None]
    extreme = max(shared, key=lambda r: abs(r["delta_q"]))
    protected = {path: digest(Path(path)) == value for path, value in manifest["protected_files_sha256"].items()}
    assert all(protected.values())
    source = run / manifest["active_source"]
    assert all(digest(source / name) == value for name, value in manifest["source_v2_sha256"].items())
    summary = {"complete": True, "created_utc": datetime.now(timezone.utc).isoformat(),
        "reference_events": 4096, "constrained_nuisances": 54, "profiled_nuisances": 57,
        "mu_hat": best["mu_hat"], "nll_min": best["nll_min"], "edm": best["edm"],
        "covariance_quality": best["covariance_quality"], "global_attempts": best["attempts"],
        "generation_mu": 1., "pruning": construction["pruning"], "workspace_sha256": construction["workspace_sha256"],
        "parent_workspace_sha256": construction["parent_workspace_sha256"], "ranking_sha256": construction["ranking_sha256"],
        "points": 42, "q_zero": rows[0]["q"], "q_two": rows[-1]["q"],
        "reference_doi": table["doi"], "reference_series": "Exp NSBI",
        "comparison_method": "41 exact shared coordinates in [0,2]; no interpolation, re-zeroing or rescaling; mu=.01 excluded from residual statistics",
        "comparison": {"shared_points": len(shared), "rms_delta_q": math.sqrt(sum(r["delta_q"]**2 for r in shared)/len(shared)),
            "maximum_absolute_delta_q": abs(extreme["delta_q"]), "mu_at_maximum": extreme["mu"]},
        "checks": {"generation_maximum_score": build_check["maximum_absolute_score"],
            "maximum_restriction_delta_nll_difference": max(abs(r["delta_nll_disagreement"]) for r in build_check["checks"]),
            "maximum_profile_delta_nll_difference": max(abs(r["delta_nll_disagreement"]) for r in validation["scan_checks"]),
            "maximum_sr_integral_difference": max(abs(r["sr_integral_relative_difference"]) for r in validation["scan_checks"]),
            "minimum_fitted_intensity": min(r["minimum_current_intensity"] for r in validation["scan_checks"]),
            "maximum_profile_edm": max(r["edm"] for r in scan),
            "minimum_pruned_minus_unpruned_q": min(r["pruned_minus_unpruned_q"] for r in rows),
            "maximum_pruned_minus_unpruned_q": max(r["pruned_minus_unpruned_q"] for r in rows),
            "protected_file_hashes_match": protected, "frozen_source_hashes_match": True},
        "recovery": global_check["recovery"],
        "profile_cost": {key: {"minimum": min(r[key] for r in rows), "median": statistics.median(r[key] for r in rows),
                              "maximum": max(r[key] for r in rows)} for key in ("fit_seconds", "fit_cpu_seconds", "process_peak_rss_gib")},
        "resources": resources, "input_sha256": inputs}
    write_csv(output / "profile_resources.csv", rows)
    write_csv(output / "resources.csv", resources["phases"])
    write_csv(output / "job_resources.csv", resources["jobs"])
    write_csv(output / "scan_validation.csv", validation["scan_checks"])
    write_csv(output / "bestfit_parameters.csv", best["parameters"])
    write_csv(output / "fixed_nuisances.csv", [{"name": n, "value": v} for n, v in construction["pruning"]["fixed_values"].items()])
    shutil.copyfile(run / "workspace/kept_nuisances.csv", output / "kept_nuisances.csv")
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    plots(rows, expected, summary, output)
    print(json.dumps({k: summary[k] for k in ("mu_hat", "q_zero", "q_two", "comparison", "profile_cost")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args())
