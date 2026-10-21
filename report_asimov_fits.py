#!/usr/bin/env python3
"""Write the completed Asimov comparison into the project's status document."""

import argparse
import json
from pathlib import Path


def completed_report(directory):
    source = json.loads((directory / "source_manifest.json").read_text())
    integrity = json.loads((directory / "input_integrity_check.json").read_text())
    summary = json.loads((directory / "run_summary.json").read_text())
    if source["pipeline_state"] != "complete" or not integrity["all_unchanged"]:
        raise ValueError("Require completed validation and unchanged inputs")
    lines = ["## Completed Asimov fits and expected comparison", "",
             "Numerical pipeline completed at " + source["completed_utc"] + ".",
             "Earlier progress entries below are historical. Visual inspection of the",
             "saved plots remains a separate manual check.", "",
             "| Reference events | Fitted mu | EDM | Maximum absolute Delta q | RMS Delta q |",
             "|---:|---:|---:|---:|---:|"]
    for size in (4096, 8192):
        row = summary[str(size)]
        if row["covariance_quality"] != 3 or row["scan_points"] != 31:
            raise ValueError("Incomplete fit or scan")
        metric = row["comparison"]
        lines.append(f"| {size} | {row['mu_hat']:.9f} | {row['edm']:.4g} | "
                     f"{metric['max_abs_delta_q']:.6g} | {metric['rms_delta_q']:.6g} |")
    lines += ["", "The reference is [HEPData Table 2, Exp NSBI](https://doi.org/10.17182/hepdata.156805.v1/t2).",
              "Both local scans profile 127 systematic parameters and three qqZZ factors.",
              "The plot spans [0,1.5]; evaluated points are 0.01 and 0.05, 0.10, ..., 1.50.",
              "Residuals use the 30 shared coordinates >=0.05, without shifts or rescaling.",
              "The published zero point is not equated with the local 0.01 proxy.", ""]
    for size in (4096, 8192):
        path = directory / f"n_{size}"
        best = json.loads((path / "bestfit.json").read_text())
        validation = json.loads((path / "likelihood_validation.json").read_text())
        if not best["complete"] or best["fit_status"] or not validation["passed"]:
            raise ValueError("Invalid final fit or validation")
        near_bounds = [p["name"] for p in best["parameters"]
                       if any(b is not None and abs(p["value"] - b) < 1e-4
                              for b in (p["min"], p["max"]))]
        integral = validation["checks"][0]["sr_integral_relative_difference"]
        residual = summary[str(size)]["maximum_absolute_delta_nll_disagreement"]
        lines += [f"For {size} events, the largest independent Delta-NLL disagreement is {residual:.4g}.",
                  f"At the global minimum the finite-reference SR integral differs from the",
                  f"model's full-MC expected SR yield by {100 * integral:+.6f}%.",
                  "Parameters within 1e-4 of a bound: " + (", ".join(near_bounds) or "none") + ".", ""]
    lines += ["Nominal same-sample normalization is independent for each subset; supplied",
              "systematic variations and full-MC yields remain unchanged. The quadrature",
              "diagnostic does not modify the model. These two nested samples do not establish",
              "a sampling distribution or isolate differences in the supplied physics inputs.", "",
              "All original and generated input ROOT files passed the unchanged-input check.",
              f"Source commit: [{source['commit'][:7]}]({source['repository']}/commit/{source['commit']}).", "",
              f"- Comparison plot: [{(directory / 'hepdata_comparison/expected_scan_comparison.png').name}]"
              f"({directory / 'hepdata_comparison/expected_scan_comparison.png'})",
              f"- Numerical comparison: [{(directory / 'hepdata_comparison/comparison_summary.json').name}]"
              f"({directory / 'hepdata_comparison/comparison_summary.json'})"]
    for size in (4096, 8192):
        path = directory / f"n_{size}" / "bestfit.json"
        lines.append(f"- {size} fitted mu, alpha, and qqZZ factors: [bestfit.json]({path})")
        path = directory / "hepdata_comparison" / f"expected_scan_comparison_{size}.png"
        lines.append(f"- {size} comparison plot: [{path.name}]({path})")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits", type=Path, required=True)
    parser.add_argument("--status-doc", type=Path, required=True)
    args = parser.parse_args()
    report = completed_report(args.fits.resolve())
    start, end = "<!-- ASIMOV_FIT_REPORT_BEGIN -->", "<!-- ASIMOV_FIT_REPORT_END -->"
    text = args.status_doc.read_text()
    block = start + "\n" + report + end + "\n"
    if start in text:
        before, remainder = text.split(start, 1)
        _, after = remainder.split(end, 1)
        text = before + block + after.lstrip("\n")
    else:
        title, remainder = text.split("\n", 1)
        text = title + "\n\n" + block + "\n" + remainder.lstrip("\n")
    temporary = args.status_doc.with_suffix(".md.tmp")
    temporary.write_text(text)
    temporary.replace(args.status_doc)
    print(report)


if __name__ == "__main__":
    main()
