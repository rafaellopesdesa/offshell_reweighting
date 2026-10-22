#!/usr/bin/env python3
"""Export every nuisance in the auxiliary-response order used by the slides."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def write_ranking(summary_path, output):
    summary = json.loads(summary_path.read_text())
    if not summary["complete"]:
        raise ValueError("A completed decomposition is required")
    samples = summary["samples"]
    keys = ("true", "data_like")
    records = {key: {r["name"]: r for r in samples[key]["impacts"]} for key in keys}
    if records["true"].keys() != records["data_like"].keys():
        raise ValueError("The two models have different nuisance inventories")
    for name in records["true"]:
        if records["true"][name]["constrained"] != records["data_like"][name]["constrained"]:
            raise ValueError("Different constraint status for " + name)
    constrained = [name for name, row in records["true"].items() if row["constrained"]]
    score = lambda name: max(abs(records[key][name]["unit_auxiliary_response"]) for key in keys)
    names = sorted(constrained, key=score, reverse=True)
    ranks = {
        key: {name: i for i, name in enumerate(sorted(constrained,
              key=lambda n: abs(records[key][n]["unit_auxiliary_response"]), reverse=True), 1)}
        for key in keys
    }
    unconstrained = sorted(set(records["true"]) - set(constrained))
    rows = []
    for i, name in enumerate(names + unconstrained, 1):
        constrained_np = name in ranks["true"]
        row = {"auxiliary_rank": i if constrained_np else None,
               "parameter": name, "group": records["true"][name]["group"],
               "has_gaussian_auxiliary": constrained_np,
               "ranking_sigma_mu": score(name) if constrained_np else None}
        for key in keys:
            r = records[key][name]
            variance = r["auxiliary_variance"]
            row.update({
                key + "_auxiliary_rank": ranks[key].get(name),
                key + "_sigma_mu": abs(r["unit_auxiliary_response"]) if constrained_np else None,
                key + "_auxiliary_variance": variance,
                key + "_percent_hessian_variance": 100 * variance / samples[key]["hessian_variance"] if constrained_np else None,
                key + "_percent_auxiliary_variance": 100 * variance / samples[key]["auxiliary_variance"] if constrained_np else None,
            })
            for field in ("unit_auxiliary_response", "value", "anchor", "pull_from_anchor",
                          "posterior_sigma", "conditional_postfit_np_impact", "correlation_with_mu"):
                row[key + "_" + field] = r[field]
        rows.append(row)
    output.mkdir(parents=True, exist_ok=True)
    filename = "all_np_uncertainty_ranking.csv"
    with (output / filename).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    notes = f"""# Complete nuisance-parameter uncertainty ranking

`{filename}` contains all {len(rows)} floating nuisance parameters for the
true and data-like 4096-event Asimov fits, side by side. It uses the saved
full observed-Hesse decomposition; no new fits or rescaling are performed.

## Ordering and columns

- `auxiliary_rank`: descending importance among the {len(names)} constrained
  nuisances, using the same order as the slide's top-12 impact plot.
- `ranking_sigma_mu`: max(abs(V_true[mu, alpha]), abs(V_data_like[mu, alpha])).
  The separate `true_auxiliary_rank` and `data_like_auxiliary_rank` columns
  give each model's own ordering.
- `*_unit_auxiliary_response`: signed local shift of mu_hat from a +1-sigma
  shift of that independent unit Gaussian auxiliary observation. A -1-sigma
  shift has the opposite linear response.
- `*_sigma_mu`: absolute auxiliary response; `*_auxiliary_variance`: its
  square. Variances add; uncertainties add in quadrature.
- `*_percent_hessian_variance`: 100 times that variance divided by V[mu,mu].
  `*_percent_auxiliary_variance` uses the auxiliary subtotal instead.
- `*_value`, `*_anchor`, `*_pull_from_anchor`, `*_posterior_sigma`: fitted
  NP value, generating auxiliary anchor, their difference, and fitted error.
  For the unit Gaussian auxiliaries, the difference is the pull in prior-sigma units.
- `*_conditional_postfit_np_impact`: V[mu,NP]/sqrt(V[NP,NP]), the conventional
  conditional one-post-fit-sigma NP impact. Its square is not an additive
  physical-source variance. `*_correlation_with_mu` is the fitted correlation.

The final {len(unconstrained)} rows are the unconstrained qqZZ normalizations.
Their auxiliary ranks, source variances and auxiliary-anchor fields are blank
(not zero): they have no Gaussian auxiliary source. Their fitted errors,
correlations and conditional impacts are retained. They remain profiled in
all calculations; their uncertainty enters through the SR/CR observations.

## Totals and interpretation

This file ranks nuisance auxiliary sources; SR/CR statistical contributions
are in the existing `true_uncertainties.csv` and `data_like_uncertainties.csv`.
The values are local linear responses, not finite-shift refits or asymmetric
profile intervals. The preserved true model normalizes nominal/up/down anchors;
the data-like model normalizes every nuisance point and shifts the generating NPs.

| Variance of mu | True Asimov | Data-like Asimov |
|---|---:|---:|
| All Gaussian auxiliaries | {samples['true']['auxiliary_variance']:.16g} | {samples['data_like']['auxiliary_variance']:.16g} |
| Observed Hesse | {samples['true']['hessian_variance']:.16g} | {samples['data_like']['hessian_variance']:.16g} |
| All primary + auxiliary sources | {samples['true']['propagated_variance']:.16g} | {samples['data_like']['propagated_variance']:.16g} |

The full source sum differs from the Hesse variance by
{100*samples['true']['closure_relative']:.4f}% (true) and
{100*samples['data_like']['closure_relative']:.4f}% (data-like). No source is rescaled.

Input: completed decomposition `summary.json`.
SHA256: `{hashlib.sha256(summary_path.read_bytes()).hexdigest()}`.
"""
    (output / "all_np_uncertainty_ranking_README.md").write_text(notes)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = write_ranking(args.summary, args.output)
    print(f"Exported {len(rows)} nuisance parameters to {args.output}")


if __name__ == "__main__":
    main()
