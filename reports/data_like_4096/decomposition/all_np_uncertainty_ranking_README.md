# Complete nuisance-parameter uncertainty ranking

`all_np_uncertainty_ranking.csv` contains all 130 floating nuisance parameters for the
true and data-like 4096-event Asimov fits, side by side. It uses the saved
full observed-Hesse decomposition; no new fits or rescaling are performed.

## Ordering and columns

- `auxiliary_rank`: descending importance among the 127 constrained
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

The final 3 rows are the unconstrained qqZZ normalizations.
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
| All Gaussian auxiliaries | 0.1389450291886888 | 0.08107689093624723 |
| Observed Hesse | 1.342558063102105 | 1.065039869433231 |
| All primary + auxiliary sources | 1.376477703380673 | 1.066873717623524 |

The full source sum differs from the Hesse variance by
2.5265% (true) and
0.1722% (data-like). No source is rescaled.

Input: completed decomposition `summary.json`.
SHA256: `9d48c48622e6d119a8c5d77c378ee65b1dd0a0135993d6dded2d8c25278c7e25`.
