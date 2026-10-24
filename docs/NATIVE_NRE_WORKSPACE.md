# Native RooFit workspace for the NRE data likelihood

This description reconstructs the graph in the supplied
`HZZ_4l_offshell_v2.root`; the original production builder was not supplied.
The file was opened read-only. `inspect_native_workspace.py` exports aggregate
structure and exact coefficient formulas without exporting observed event values.
The inspection is saved with the structured-summary reports.

## Precompute neural outputs

Simulation-trained ratios estimate a process density relative to a common,
parameter-independent reference. Nominal ratios and systematic up/down
predictions are evaluated at each observed event and stored as RooDataSet
columns. The workspace fits analytic combinations of these numbers; it does
not contain a neural-network inference stage that runs during profiling.
The common reference-density factor cancels from likelihood differences.

The concrete dataset `obsData_4l_offshell` has 2074 unit-weight signal-region
rows and one weighted row for each control region: 1377, 231 and 155 observed
counts. There are 2554 stored columns: ten nominal ratios, 2540 variation
anchors (two directions times 127 NPs times ten terms), three CR observables
and the channel category. These stored columns are distinct from the thirteen
variables listed in ModelConfig's observable set. The learned shapes are
normalized using simulation; the observed events are not a normalization
reference sample.

## Build process functions and the extended signal-region sum

For a process/basis term p, the graph implements

```
V_pi(alpha) = piecewise_interp_p         # per-event anchors, interpolation code 5
R_p(alpha)  = mod_p                     # integrated rate modifier, code 4
h_pi        = r_pi * V_pi / R_p
C_p         = nominal_yield_p * kappa_p(mu, beta) * R_p
I_i         = sum_p C_p * h_pi
Lambda_SR   = sum_p C_p
```

`RooProduct` and `RooRatio` form `prod_process_p`; a `RooWrapperPdf` exposes
that function as `pdf_process_p`. `RooAddPdf("pdfSRUnc",...)` combines ten
components with extended coefficients `nfFactor_p`. A unit-valued
`dummyFactor` in the native graph maintains observable dependencies and
leaves the numerical intensity unchanged.

The wrappers return `selfNormalized() == true`. RooFit therefore trusts the
normalization supplied by the input construction rather than adding an
integration over the ratio-coordinate space. See the
[ROOT 6.34 RooWrapperPdf reference](https://root.cern.ch/doc/v634/classRooWrapperPdf.html).
In particular, an observed dataset is not used to recalibrate these shapes.

The gg basis is S, SBI and B with coefficients mu-sqrt(mu), sqrt(mu) and
1-sqrt(mu). This represents the usual signal, interference and background
combination using the stored basis predictions. EW uses EWB, EWSBI and
EWSBI10 with the formulas recorded in the inventory. The other terms are
qqZZ in three jet bins and ttV. Coefficients can be signed; the combined
intensity and expected counts must stay positive.

## Add channels and constraints

Each CR combines eight `RooUniform` one-bin components with their expected
yields in a `RooAddPdf`. A single row weighted by the observed count gives the
counting likelihood when the extended term is included. It contributes
`Lambda_c - n_c*log(Lambda_c)` up to constants.

The 127 `RooGaussian` constraints have unit width and fixed global observations
at zero in the original data file. A shared `RooProdPdf("constraintPdf",...)`
is attached to the SR and CR channel PDFs. Constraint extraction in the NLL
includes this shared factor once. `RooSimultaneous("simPdf",...,cat)` selects
SR, CR0, CR1 or CR2 for each row.

Three unconstrained qqZZ parameters are cumulative: effective jet-bin scales
are beta0, beta0*beta1 and beta0*beta1*beta2, with the same factors in the
corresponding CRs. They have no Gaussian auxiliary constraints.

The complete NLL, omitting parameter-independent constants, is

```
Lambda_SR - sum_SR log(I_i)
+ sum_CR [Lambda_c - n_c log(Lambda_c)]
+ 0.5 * sum_k (alpha_k - a_k)^2.
```

`root_workspace_validation.py` independently assembles this expression and
checks likelihood differences against RooFit.

## Serialize and fit

The `combined` RooWorkspace owns the graph, parameters and dataset.
`ModelConfig_4l_offshell` points to `simPdf` and declares observables, nuisances,
global observations and POIs. Its serialized POI set contains mu and the three
qqZZ factors; its nuisance set contains the 127 alphas. Our statistical study
scans only mu and profiles all 130 remaining parameters. The native input
bounds and fixed global observations are retained.

`root_workspace_fit.py` creates the extended NLL with the ModelConfig nuisance
and global-observation sets, model-sourced global values and the CPU evaluation
backend. The original global fit returned mu_hat=0.8731396397. This workspace
and its data were never rebuilt for the observed fit.

## Where RooSameSampleNormalizer enters

A finite-reference Asimov model needs exact integral closure on its chosen
simulation reference at every tested nuisance point. The adaptation replaces
`rhat*V/R` by `rhat*V/Z(alpha)`, with
`Z(alpha)=sum_reference u*rhat*V`, while retaining the physical coefficient
`C_p`, including its intended rate modifier R. Generate weighted SR rows, CR
counts and auxiliary observations consistently at the chosen parameters.
This changes the finite-Asimov construction, not the preserved observed-data
workspace. The data-like construction injects only observed-fit NPs; mu_A
remains one. See [ROO_SAME_SAMPLE_NORMALIZER.md](ROO_SAME_SAMPLE_NORMALIZER.md).
