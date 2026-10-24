# RooSameSampleNormalizer

`RooSameSampleNormalizer` is a scalar `RooAbsReal` that integrates one process's
interpolated event contribution over a fixed reference sample. Every current
true/data-like/pruned workspace at 4096 or 8192 uses ten instances, one per SR component, to keep each process
shape normalized as nuisance parameters change during fitting.

The implementation is in [RooSameSampleNormalizer.h](../RooSameSampleNormalizer.h)
and [RooSameSampleNormalizer.cxx](../RooSameSampleNormalizer.cxx). The builder is
[build_shifted_asimov.py](../build_shifted_asimov.py), and
[same_sample_normalizer.py](../same_sample_normalizer.py) loads its ROOT dictionary.

## What the class calculates

For process p and the selected reference events i, define

$$
u_i=\frac{w_i}{\sum_j w_j},\qquad
m_{pi}=u_i\widehat r_{pi},\qquad
V_{pi}(\boldsymbol\alpha)=\text{the event-level systematic interpolation}.
$$

The class returns

$$
Z_p^{\mathrm{raw}}(\boldsymbol\alpha)
=\sum_i m_{pi}V_{pi}(\boldsymbol\alpha).
$$

The measure m is fixed when the class is constructed. It uses the reference
weights and nominal normalized ratios, **never the generated Asimov weights**.
The up/down event anchors are also fixed. Only the alpha values change during
likelihood evaluation.

Normalizing nominal ratios and individual up/down anchors does not ensure that
the full shape remains normalized between anchors or when several nuisances
vary together. This class evaluates that missing integral at each trial point.

## Connection to the existing RooFit graph

The original process function has the form

$$
h_{pi}=\texttt{prod\_process}_p(x_i;\boldsymbol\alpha)
=\widehat r_{pi}\frac{V_{pi}(\boldsymbol\alpha)}{R_p(\boldsymbol\alpha)},
$$

where `mod_<process>` supplies the physical rate modifier R. The original
`dummyFactor` equals one and retains the observable dependencies.

The builder clones h under `raw_process_<process>` and installs:

| RooFit node | Meaning |
| --- | --- |
| `reference_Z_<process>` | This class: sum of u × rhat × V |
| `normalized_numerator_<process>` | `raw_process_<process>` × `mod_<process>` = rhat × V |
| `prod_process_<process>` | `RooRatio`: numerator / `reference_Z_<process>` |
| `nfFactor_<process>` | Existing physical yield coefficient, unchanged |

Thus the replacement shape is

$$
\overline h_{pi}=\frac{\widehat r_{pi}V_{pi}}{Z_p^{\mathrm{raw}}},
\qquad \sum_i u_i\overline h_{pi}=1.
$$

This is algebraically the proposed h / sum(u h) construction. Using the raw
interpolation integral lets the class avoid depending on the separate rate node.
The physical coefficient still contains R, so normalization uncertainties in
the supplied full-MC rates are preserved.

The ten components are S, SBI, B, EWB, EWSBI, EWSBI10, ttV, and qqZZ_0/1/2.
Each qqZZ jet component gets its own integral. The existing cumulative free
normalizations remain beta0, beta0 × beta1, and beta0 × beta1 × beta2.
Neither mu nor these three free factors is a server of the normalizer: they
enter the physical coefficients outside it.

## Interpolation convention

This implementation supports the convention in the supplied workspace:

- `PiecewiseInterpolation`, code 5 for every event nuisance;
- nominal value 1 and interpolation boundary at ±1;
- strictly positive up/down anchors;
- no positive-definite clipping of the interpolated value.

The constructor receives an ordered `RooArgList` of nuisances, an N-element
measure vector, and two flattened K × N anchor arrays. The flattened index
is `nuisance_index * N + event_index`. The builder takes all three ordered
lists directly from the existing interpolation node, checks the convention,
and checks dimensions and finite input values.

For each nuisance, code 5 uses a sixth-order polynomial between −1 and +1
and power-law extrapolation outside that interval. The full event modifier
multiplies these contributions. The optimized polynomial coefficients match
ROOT 6.34.10's
[`MathFuncs::flexibleInterpSingle`](https://github.com/root-project/root/blob/v6-34-10/roofit/roofitcore/inc/RooFit/Detail/MathFuncs.h).
The native
[`PiecewiseInterpolation::evaluate`](https://github.com/root-project/root/blob/v6-34-10/roofit/histfactory/src/PiecewiseInterpolation.cxx)
uses the same function. This is a specialized implementation, not a replacement
for arbitrary RooFit interpolation codes. Revalidate it if that convention or
the ROOT version changes.

## RooFit dependencies, derivatives, and caching

`parameters_` is a `RooListProxy` containing the nuisance parameters. RooFit
therefore sees the normalizer as dependent on alpha, and invalidates its scalar
value when those servers change. There are no event-observable servers in this
node, so evaluating successive dataset rows does not change its integral.

The normalizer is shared by every event of that process. Its fixed arrays are
stored once per class instance. ROOT may clone instances for its evaluation
backend; clones have independent transient caches.

The implementation also caches:

1. Polynomial coefficients and logarithms derived from the fixed anchors.
2. Each nuisance's last value and per-event interpolation factors.
3. The product of those factors for each reference event.

When only one nuisance changes, it updates that nuisance's factors and the
event products, then takes the weighted sum. It rebuilds the complete products
every 128 evaluations, and when an old factor is zero or nonfinite, to limit
roundoff accumulation and avoid dividing by such a factor. No anchors or
intensities are clipped. An invalid likelihood point must be handled by the
fit's existing validity checks.

There is no separate analytic derivative implementation. Minuit's numerical
derivatives reevaluate the graph at changed alpha values, including Z(alpha).
The cache is an evaluation optimization; it does **not** hold Z fixed during
minimization. Changing only mu leaves Z unchanged, as required by its definition.

`directValue()` supplies a slower validation calculation: it loops over all
reference events and nuisances using ROOT's interpolation function directly,
without the optimization caches.

## Constructing the data-like Asimov

The builder opens the preserved true-Asimov workspace read-only, takes its
existing 4096 rows, normalized reference weights, ratios and corrected anchors,
and builds a separate workspace with the normalized shapes. It then:

1. Sets mu = 1, all 127 alpha values and all three cumulative qqZZ factors to
   the observed-data global-fit nuisance values.
2. Evaluates the SR intensity I at this point and saves fixed weights a_i = u_i I_i.
3. Sets each CR count to its expectation at the same generating point.
4. Sets each Gaussian global observable `nom_<name>` to the injected
   `alpha_<name>`, and keeps these global observables constant.
5. Saves the data and `asimov_generation` snapshot. All 131 model parameters
   are free in the global fit; a profile fixes only mu.

The normalization reference remains fixed throughout the fit, while the model
normalizers change with the trial nuisances. Changing only the generating
parameters and Asimov weights can reuse this reference. Changing the selected
reference rows, reference weights, nominal ratios, or event anchors requires
rebuilding the normalizers.

For the nonnegative 4096 reference weights and positive combined intensities,
this construction makes the generating point a global likelihood minimum,
including consistently generated CR and auxiliary data. It does not guarantee
uniqueness or accurate approximation of the full-sample curvature and tails.
The nonnegative-divergence argument does not apply to the signed 8192 sample.
Negative signal/interference basis coefficients are allowed; the combined
intensity entering each logarithm must be positive.

## Serialization and loading

The current [six-workspace family](SAME_SAMPLE_WORKSPACE_FAMILY.md) uses this
same version-1 class for true, data-like and pruned data-like models at both
4096 and 8192 events. The original anchor-only true models and their fitted
results remain separately preserved. No class or serialization-version change
was needed; the builder now supports nominal generation and signed 8192 inputs.

The `ClassDefOverride(..., 1)` dictionary streams the proxy, fixed measure,
and anchor arrays into the ROOT file. Polynomial and evaluation caches are
transient and reconstructed after loading or cloning. The builder also calls
`RooWorkspace::importClassCode` so the workspace carries the custom source.
See ROOT's [workspace documentation](https://root.cern.ch/doc/v634/classRooWorkspace.html).

For the project workflow, load the dictionary explicitly before reading the
workspace. From the analysis repository, with the local ROOT environment
active and a Slurm allocation available:

```python
import ROOT
from same_sample_normalizer import load_normalizer

load_normalizer()
source = ROOT.TFile.Open("/path/to/HZZ_4l_data_like_asimov_4096.root", "READ")
w = source.Get("combined")
w.loadSnapshot("asimov_generation")
z = w.function("reference_Z_S")
print(z.getVal(), z.directValue())
# Keep source open while using the workspace.
```

The loader uses ACLiC and a cache keyed by the header/source SHA256. Set
`OFFSHELL_ROOT_BUILD` to a persistent project directory in `/work` for batch
jobs; otherwise it defaults to `/tmp/offshell_root_build`. A file lock prevents
concurrent compilation of the same library. It distinguishes a real dictionary
from an emulated class that ROOT may create while inspecting file metadata.

`root_workspace_fit.py` recognizes the metadata field
`normalization_implementation = RooSameSampleNormalizer_v1` and loads the class
before retrieving the workspace. Older workspaces continue through their
original path. ROOT can also compile the embedded class when opened without
the project loader, provided a compiler and writable build location are available.

## Validation and operational limits

[validate_shifted_asimov.py](../validate_shifted_asimov.py) checks a saved,
reopened workspace, including:

- Exact injected parameters, fixed globals, unchanged reference rows and weights.
- Regenerated SR weights and CR counts versus the saved dataset.
- Per-process integrals at the observed-nuisance injection and other generating
  points: mu = 0.5 with a luminosity shift, mu = 1.5 with an individual shape
  shift beyond an anchor, and mu = 1.5 with all observed-fit nuisances.
- Cached versus direct ROOT interpolation integrals, also after derivative updates.
- All 131 likelihood scores at the main injection, using two finite-difference
  steps and Richardson extrapolation; selected scores at the other points.
- Recovery of the injected parameters and likelihood by the global fit.

The construction checks passed with maximum absolute score 7.58 × 10⁻⁸ and
maximum process-integral error 4.67 × 10⁻¹⁵. These are construction checks;
profile convergence and agreement with HEPData are separate results.

[root_workspace_validation.py](../root_workspace_validation.py) independently
assembles the extended SR/CR/Gaussian likelihood for every saved profile point.
For the new model it also requires relative SR integral closure below 3 × 10⁻¹².
It checks positivity at the actual fitted point; a fitted nuisance vector is
not required to remain physical at every other mu.

The class has mutable caches and is used here with the single-CPU RooFit CPU
backend. The same instance has not been validated for concurrent evaluation
from multiple threads. All substantial construction, validation and fitting
runs through Slurm. The custom arrays and backend clones increase memory use;
use measured batch accounting when sizing later jobs.

### Measured evaluation overhead

An alternating true/data-like/data-like/true benchmark on one CPU of the same
Slurm allocation compared the native likelihoods at identical trial parameter
sequences. Initialization and initial cache evaluations were excluded.

| Updated parameters | True model, mean CPU ms | Data-like model, mean CPU ms | New / old |
| --- | ---: | ---: | ---: |
| One alpha | 159.548 | 161.276 | 1.011 |
| All 127 alpha | 159.577 | 196.100 | 1.229 |
| Only mu | 0.08157 | 0.08092 | 0.992 |

The sub-percent mu-only difference is timing variation, not evidence of a
speed improvement. Peak process RSS in these evaluation runs was 1.10 versus
1.63 GiB. The datasets have different fixed Asimov weights, so this compares
practical model-evaluation cost rather than isolating the C++ class alone.

The common numerical-gradient pattern changes one or a few nuisances per
evaluation. The expensive all-nuisance update occurs less often in this
minimization workflow. A complete profile also depends on how many likelihood
evaluations Minuit needs; the two models can have different curvatures and
iteration counts. Use the separately recorded single-profile and full-scan
timings to assess that cost. The scripts are
[benchmark_normalized_nll.py](../benchmark_normalized_nll.py) and
[report_normalizer_timing.py](../report_normalizer_timing.py).

A complete mu=0.5 profile on that same CPU, initialized from each model's
global-fit nuisance vector, gave:

| Model | Minimization time | Likelihood evaluations | Seconds / evaluation |
| --- | ---: | ---: | ---: |
| Preserved true Asimov | 251.99 s | 1,579 | 0.15959 |
| Data-like Asimov | 304.10 s | 1,841 | 0.16518 |

The complete profile was 20.7% slower, with 16.6% more likelihood evaluations
and 3.5% more time per evaluation. Both points converged and passed independent
likelihood validation. This is one coordinate and two different generated
models, so it is not a universal profile-time multiplier. The raw aggregate
measurements are in [timing.json](../reports/data_like_4096/timing.json).

### Minimum validation versus covariance validation

A strategy-0 fit started at an almost exact stationary point can return a
converged likelihood with an inadequate covariance estimate. This happened
for the initial mu=1 start: quality 1 and a non-positive-definite saved matrix.
The covariance validator rejected it. A displaced mu=1.4 start converged with
quality 3 and was repeated to save that result as the profile baseline.
The initial artifacts were preserved. This is separate from the normalization
or score checks, and does not justify accepting a non-positive-definite
covariance. The profile baseline was saved without an explicit full Hesse.
The subsequent requested [uncertainty decomposition](UNCERTAINTY_DECOMPOSITION.md)
computes separate full observed-Hesse matrices and validates their curvature;
those results do not replace the scan baseline.

The driver now preserves each attempt's ROOT fit result. For future multistart
runs, it prefers quality 3 when candidate NLLs differ by at most 1e-7; otherwise
the lower NLL wins. That numerical tie tolerance is much smaller than the scan
agreement tolerance. The baseline always retains the selected fit's actual
parameters and NLL.

The true-Asimov files under `results/asimov_anchor_normalized/` are preserved.
They retain the previous nominal/anchor normalization. Therefore the new
true-versus-data-like comparison includes both the nuisance shift and this
refinement of the finite-reference model.

## Small executable example

[`normalizer_usage_example.py`](../normalizer_usage_example.py) checks a
two-event constructor, both anchors, an interior nuisance point and cache
updates. Run with the project's loader, or pass `--library` with the path to
the standalone shared library. The structured presentation shows the same
example alongside compilation and read-only workspace-loading snippets.
See [STRUCTURED_SUMMARY.md](STRUCTURED_SUMMARY.md) for the full workflow.
