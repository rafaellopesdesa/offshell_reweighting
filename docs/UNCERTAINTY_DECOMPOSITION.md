# Covariance-based uncertainty decomposition

The original 2026-10-22 study compared the then-preserved anchor-normalized
true 4096 Asimov and the fully normalized data-like 4096 Asimov. Their profile-scan baselines and input workspaces remain unchanged.
The numerical table and ranking are exported under
[`reports/data_like_4096/decomposition/`](../reports/data_like_4096/decomposition/).

## What the table measures

We follow [arXiv:2608.15312](https://arxiv.org/abs/2608.15312), particularly
Eqs. (22), (30), (44), (51), and (81). With covariance `V` and independent
unit Gaussian auxiliary observations, a unit shift of auxiliary observation
`a_k` induces the local response `delta mu = V[mu, alpha_k]`. Its variance
contribution is `V[mu, alpha_k]^2`. Group variances add; standard deviations
add in quadrature.

The conventional conditional impact is instead
`V[mu, alpha_k] / sqrt(V[alpha_k, alpha_k])` for a one-post-fit-sigma
nuisance displacement. These named-coordinate impacts are not an additive
source decomposition. An unconstrained qqZZ parameter has no corresponding
Gaussian auxiliary source. It remains profiled in every response calculation.

These are local linear-response uncertainties. They are neither finite
one-sigma refits nor asymmetric profile-likelihood intervals. The ranking
uses an impact-plot layout but displays auxiliary-source responses.

## Covariance and primary observations in this implementation

`compute_fit_covariance.py` restores the saved global fit and performs the
Migrad initialization required by RooMinimizer before a full strategy-1
Hesse calculation. Results go in a separate directory. It records the small
parameter/NLL change, status, positive definiteness, parameter ordering,
ROOT version, source hashes, calls and time. No SumW2 correction is applied.
The original MIGRAD estimate is not used for the decomposition.

`compute_observation_response.py` evaluates the native workspace's SR
intensities and CR expectations at precisely that refined point. The probe direction is the complete covariance column `V[:, mu]`, including
all three cumulative qqZZ parameters. Four log-intensity evaluations at
`point +/- h V[:, mu]`, with `h=1e-4` and `h/2`, give the observation response
by Richardson extrapolation. This avoids constructing 131 event-gradient
columns while evaluating the same scalar variance. Step-size stability is
checked separately for every region.

For each region `s`, the code constructs

```
B_s = sum_{i in s} a_i grad(log I_i) grad(log I_i)^T
variance_mu[s] = V[mu, :] B_s V[:, mu]
```

Here `a_i` is the fixed positive Asimov weight: a quadrature mass of the
generating Poisson intensity. The 4096 reference rows do not each represent
one experimental event. The CR row weights are the corresponding generating
Poisson means. Squared MC weights would describe a different experiment;
we do not use them as the observation variances. Finite-reference MC
uncertainty is not an additional fitted source in these workspaces.

The script inspects the actual model graph to require one independent,
fixed-width unit Gaussian for every alpha and a fixed matching global
observable. It fails if this assumption is not satisfied.

## Closure and checks

The report uses the **observed Hessian**, evaluated near the saved optimum.
It reports both the propagated source sum and `V[mu, mu]`, and their relative
difference. It does not rescale the components to make the totals agree.
Exact closure requires the source score covariance to equal the information
used for `V`. In that original study, the preserved true model only normalized its anchors;
the data-like model normalized every nuisance point. The current six-workspace
family now uses RooSameSampleNormalizer throughout. Observed curvature,
generation/fit differences and this model distinction can affect closure.

`validate_covariance_curvature.py` independently checks four native-NLL
directions against the saved covariance: the profiled-mu tangent, the
largest covariance eigenvector and two reproducible whitened directions.
The finite differences use two scales and Richardson extrapolation.
`test_uncertainty_decomposition.py` checks an analytically soluble Gaussian
experiment with both a constrained nuisance and an unconstrained background;
it distinguishes source propagation from freezing the nuisance.

All numerical work uses Slurm. Run the covariance calculation first, the
observation-response calculation at its saved refinement point, then the
directional check and `report_uncertainty_decomposition.py`. Preserve the
per-model covariance JSON, ROOT result, response JSON and validation JSON
alongside the Slurm logs. The public aggregates include all named nuisances,
explicit grouping, conditional-impact diagnostics and measured resource use.


## Completed 4096 study

The [results report](DATA_LIKE_ASIMOV_RESULTS.md) gives the numerical table.
The source/Hesse variance mismatch is 2.5265% for the preserved true model
and 0.1722% for the data-like model. Dividing the propagated source sum by
the independently evaluated quadratic form along `V[:,mu]` instead of its
Hesse prediction leaves 2.3030% and 0.00551%. This is a closure diagnostic;
it does not replace the saved covariance or rescale any source component.
Both totals remain visible in the tables and slides.


## Complete nuisance ranking download

The [combined CSV](../reports/data_like_4096/decomposition/all_np_uncertainty_ranking.csv)
contains both models and all 130 NPs. The 127 Gaussian auxiliary sources are
sorted by the larger absolute response across the models, exactly as in the
slide plot; each model's individual rank is also included. The three
unconstrained qqZZ parameters follow with blank auxiliary-source fields
and their fitted errors and conditional impacts retained.
[Column definitions and totals](../reports/data_like_4096/decomposition/all_np_uncertainty_ranking_README.md)
accompany the table. Regenerate it from saved results without ROOT or refitting:

```bash
python export_np_ranking.py \
  --summary reports/data_like_4096/decomposition/summary.json \
  --output reports/data_like_4096/decomposition
```

## Structured-summary update: unpruned data-like 4096 and 8192

The 2026-10-24 report uses the current unpruned data-like workspaces at both
reference sizes. Generating mu remains one. It reuses the validated 4096 full
Hesse only after verifying identical workspace bytes and the same optimum.
A separate full Hesse and independent directional checks are run for 8192.
Neither refinement replaces any published scan baseline.

The 8192 sample contains two negative reference weights. An explicit
`--allow-signed-quadrature` option permits the existing weighted sums to
approximate the population Poisson score integral. Negative coefficients are
retained as integration terms, not interpreted as negative Poisson means.
Replacing them with absolute or squared weights would define a different
calculation. The output records positive and negative contributions separately,
and requires finite nonnegative aggregate estimates and derivative stability.
`test_signed_score_quadrature.py` verifies exact polynomial integration with
both positive and signed rules and checks that no clipping occurs.

`report_summary_diagnostics.py` provides two pull origins and two uncertainty
views. The original-nominal view shows fitted alpha and injected alpha_A; the
Asimov-centered view shows their difference. Unconstrained qqZZ parameters
remain in their raw factor coordinates. Fifteen auxiliary sources are shown
in one shared order across both sample sizes, while the complete 130-NP
tables retain model-specific ranks, conditional impacts and correlations.
The grouped plot includes all primary regions and auxiliary groups, with
propagated and Hesse totals separately displayed. Adjacent paper-figure slots
are configured in the presentation branch. See
[STRUCTURED_SUMMARY.md](STRUCTURED_SUMMARY.md) for reproducibility and costs.
