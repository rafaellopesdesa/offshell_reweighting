# Scientific summary and six-model resource accounting

The presentation now follows the finite-Asimov construction, its RooFit
implementation, validated likelihood results and measured computing cost.
The original research chronology is retained separately in the presentation
branch as `presentation/history-20261023.html`.

## Native workspace primer

Four slides before the custom class reconstruct the supplied observed-data
workspace: precomputed neural outputs, shape/rate factors, channel and
constraint composition, then serialization/ModelConfig. The inspection is
read-only and exports aggregate graph structure without event values. See
[NATIVE_NRE_WORKSPACE.md](NATIVE_NRE_WORKSPACE.md).

## Method and executable example

The opening reference is [Finite Asimov Sample Construction in Unbinned Neural
Simulation-Based Inference](https://arxiv.org/abs/2609.14136). The later
[hNDE study](https://arxiv.org/abs/2609.30196) supplies the shape-systematic
interpolation example. Same-sample normalization uses the fixed reference
measure for both generation and every tested parameter point. Our physical
rate modifiers remain outside the shape normalizers. See
[ROO_SAME_SAMPLE_NORMALIZER.md](ROO_SAME_SAMPLE_NORMALIZER.md).

`normalizer_usage_example.py` constructs a two-event, one-NP normalizer and
checks both anchors, an interior point and cache updates against `directValue`.
At alpha=1, Z=.4*1.2+.6*.8=.96; normalized shapes are 1.25 and 5/6. Their
reference integral is one, while an external coefficient of 100 retains a
yield of 100. With the standalone demo library built using ROOT >= 6.34.04:

```bash
python normalizer_usage_example.py --library /path/to/libRooSameSampleNormalizer.so
```

The standalone directory is on branch `demo/colab-workspace-scans`, where
`python normalizer/build.py --output build/normalizer` builds both library and
ROOT dictionary. Build and load with the ROOT runtime used to deserialize the
workspace. The source ROOT file must remain open during use.

## Likelihood and nuisance diagnostics

All six existing 61-point scans at mu=0,.05,...,3 are retained byte-for-byte,
including the continuous global-fit baseline and original HEPData values.
[The comparison report](SIX_SCAN_HEPDATA_COMPARISONS.md) gives their validation.
All current workspaces use the full RooSameSampleNormalizer. Original files,
reference selections, rate effects and pruned-NP choices remain unchanged.

The two diagnostic subjects are the **unpruned data-like Asimov models**, with
mu_A=1 and observed-data alpha/qqZZ injections, at 4096 and 8192 events.
They are distinct from the original observed dataset. For the 4096 model, the
full Hesse is reused after matching the workspace hash and optimum: the current
scan baseline differs by only 2.89e-9 in NLL and 9.16e-5 in any parameter.
The 8192 model receives a separate full Hesse. Scan baselines are not replaced.

The full-covariance propagation, signed-quadrature policy and distinction
between source responses and conditional impacts are documented in
[UNCERTAINTY_DECOMPOSITION.md](UNCERTAINTY_DECOMPOSITION.md). Results are generated
only after Hesse, parameter/likelihood identity, derivative-step stability and
independent curvature checks pass. Both source and Hesse totals are retained.
No rescaling enforces their agreement.

The pull plot can show original nominal coordinates or alpha_hat-alpha_A.
Free qqZZ factors are shown in raw coordinates because they have no Gaussian
constraint. The uncertainty plot switches between individual auxiliary
responses and all source groups. All 130 NPs are downloadable; the plot uses
15 leading Gaussian sources in a shared order across the two sizes. The
presentation reserves four adjacent slots for the author's paper figures.

## Resource accounting

`report_summary_resources.py` reads the saved per-point wall/CPU/RSS fields,
checks their agreement with the final selected profiles and original CSV
hashes, and queries Slurm accounting for all 36 actual task IDs. Slurm array
display IDs and actual job IDs are recorded separately in the local audit.
Public tables use anonymous worker numbers and omit hosts and execution
timestamps while retaining all requested resource measurements. The six plots use
common axes and contain all 61 coordinates, including mu=0 and mu=1.

| Model | Median wall min/point | Maximum wall min/point | Max process RSS GiB | Max Slurm step RSS GiB |
|---|---:|---:|---:|---:|
| True 4096 | 4.62 | 10.27 | 1.632 | 1.347 |
| Data-like 4096 | 4.55 | 10.22 | 1.633 | 1.348 |
| Data-like pruned 4096 | 2.37 | 3.72 | 1.613 | 1.538 |
| True 8192 | 9.16 | 18.71 | 2.378 | 2.550 |
| Data-like 8192 | 10.43 | 20.29 | 2.379 | 2.548 |
| Data-like pruned 8192 | 4.72 | 8.13 | 2.360 | 2.530 |

Each profile uses one CPU. Point timing measures minimization, including a
strategy retry if required, and excludes loading/initial NLL setup. Profiles
use warm starts within six chunks; a new branch can give a time spike.
Process RSS is a high-water mark through that point, not an isolated allocation
at each mu. Slurm MaxRSS has a different scope and sampling frequency and is
reported separately, without summing concurrent jobs' peaks.

There are 66 computations per model because chunk boundaries overlap. The
aggregate minimizer CPU includes all 66; whole-job CPU includes loading and
uses parent records without double-counting job steps. Scan elapsed spans the
first chunk start through the last finish, includes scheduling gaps and has
at most four concurrent chunks per model. Global fits, final validation and
separate Hesse calculations are excluded from this scan-cost comparison.

The author confirms the historical **24 hours and 500 GB were for ONE profiled
mu point**. Hardware, core count and requested-versus-measured RAM are unknown.
The plots convert 500 decimal GB to about 466 GiB and show this as a historical
context, not a controlled speedup or identical-hardware benchmark.

## Reproduce the exports

Use a suitable Slurm allocation for numerical work and plotting. From the
work-area root, after both diagnostic directories are complete:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
.venv_presentation/bin/python offshell_reweighting/report_summary_resources.py \
  --run results/asimov_scans_0_3_61 \
  --output results/summary_presentation_20261024/resources
.venv_presentation/bin/python offshell_reweighting/report_summary_diagnostics.py \
  --results results/summary_presentation_20261024/diagnostics \
  --data-fit results/data_workspace/bestfit.json \
  --output results/summary_presentation_20261024/diagnostics_report
.venv_presentation/bin/python results_slides/presentation/tools/export_summary_results.py \
  --run results/summary_presentation_20261024 --analysis offshell_reweighting
```

The presentation exporter requires complete reports and matching scan hashes,
records the generated artifact hashes, and preserves the original scan assets.
`check_summary_presentation.py` checks all slides, diagnostic switches, local
downloads/checksums, navigation, mobile deep links and the archived deck.
