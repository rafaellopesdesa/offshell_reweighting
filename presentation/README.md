# Structured scientific summary

Open `index.html` through a local HTTP server or the repository's GitHub Pages
site. The 33-slide Reveal.js deck follows four parts:

1. Finite Asimov construction, exact same-sample closure, and shape systematics.
2. The native NRE RooFit workspace: ratio columns, shapes/rates, channels and
   ModelConfig; then `RooSameSampleNormalizer`, loading and a tested constructor.
3. Six consecutive 61-point HEPData comparisons, then data-like unpruned 4096/8192
   nuisance pulls and covariance source uncertainties.
4. Wall time, CPU and memory versus μ for all six models, followed by full costs.

The earlier 35-slide research history is preserved as `history-20261023.html`,
using its original JavaScript, styles and numerical assets. The main deck uses
`summary.css` and `summary.js`. Below 700 px it becomes a continuous reading page;
both `#/slide-id` and reading-mode anchors work. Contents and keyboard navigation
are supported. All figures and results are static local files.

## Add the paper figures

Put the four images in `assets/` and edit `data/paper_figures.js`:

```javascript
window.PAPER_FIGURES = {
  pulls_4096: {src: 'assets/paper_pulls_4096.svg', alt: 'Paper nuisance pull plot', caption: 'Paper result'},
  uncertainties_4096: {src: 'assets/paper_uncertainties_4096.svg', alt: 'Paper uncertainty decomposition'},
  pulls_8192: {src: 'assets/paper_pulls_8192.svg', alt: 'Paper nuisance pull plot'},
  uncertainties_8192: {src: 'assets/paper_uncertainties_8192.svg', alt: 'Paper uncertainty decomposition'}
};
```

Unconfigured slots remain visibly reserved. The export script does not touch this
configuration. SVG, PNG or JPEG images work. No scientific numbers are inferred
from the paper images. Each figure sits beside our corresponding diagnostic.

## Refresh the summary

Run report generation in a suitable Slurm allocation. ROOT >= 6.34.04 with RooFit
is required for the covariance calculation and executable class example. The
report/export steps use Python, NumPy and Matplotlib; browser checks use Playwright.
From the work-area root, with both covariance reports independently validated:

```bash
.venv_presentation/bin/python offshell_reweighting/report_summary_resources.py --help
.venv_presentation/bin/python offshell_reweighting/report_summary_diagnostics.py \
  --results results/summary_presentation_20261024/diagnostics \
  --data-fit results/data_workspace/bestfit.json \
  --output results/summary_presentation_20261024/diagnostics_report
.venv_presentation/bin/python results_slides/presentation/tools/export_summary_results.py \
  --run results/summary_presentation_20261024 --analysis offshell_reweighting
PLAYWRIGHT_BROWSERS_PATH="$PWD/.presentation_browsers" \
  .venv_presentation/bin/python results_slides/presentation/tools/check_summary_presentation.py \
  --output results/summary_presentation_20261024/browser
```

The exporter requires both complete covariance reports and six resource reports
whose scan checksums match `family_snapshot.json`. It copies only `summary_*`
assets and records artifact/source hashes. Existing scans, baselines and HEPData
values are unchanged. The original six-model export remains
`tools/export_family_results.py`.

## Scientific and resource definitions

- Original-nominal pull view: fitted α with full Hesse error; diamonds show the
  observed-data values used to generate this Asimov. Asimov-centered view subtracts
  those centers. The three unconstrained qqZZ factors remain in raw units.
- Plot labels abbreviate electron `SIMPLIFIED_UncorrUncertainty` to `uncorr`
  and `CorrUncertainty` to `corr`; the CSV preserves every exact RooFit name.
- The CSV field `posterior_sigma` is the local post-fit Hesse error, not a
  Bayesian credible interval.
- Auxiliary source response: `V[mu,alpha]` for a unit Gaussian observation shift.
  Its square is the additive source variance. `V[mu,alpha]/sigma_alpha` is a
  different conditional diagnostic, available in the all-NP CSV.
- Primary sources: the intensity-score covariance propagated through the full V.
  The two negative 8192 quadrature weights are retained as integration coefficients,
  not interpreted as negative Poisson variances. Source/Hesse closure is measured
  and reported without rescaling. These are local errors, not profile intervals.
- Each point uses one CPU. Timing curves measure minimization and exclude loading.
  Process RSS is a high-water mark up through that point. Slurm MaxRSS is separately
  reported for whole jobs; sampling and process scope differ. Public records use
  anonymous worker numbers; detailed cluster identifiers remain in the local audit.
- 61 unique points require 66 computations per model because of overlapping chunks.
  Aggregate CPU includes all 66. Whole-job CPU includes loading; scan elapsed spans
  the first start to last finish, with scheduling gaps and at most four simultaneous
  chunks per model. Global fits, final validation and Hesse are separate work.
- The author's historical **24 h / 500 GB is per profiled μ point**, not per scan.
  Hardware, core count and RAM accounting are unspecified. It is a contextual cost
  comparison, not a controlled speedup. Decimal 500 GB is about 466 GiB.

Primary references: [finite Asimov construction](https://arxiv.org/abs/2609.14136),
[hNDE shape example](https://arxiv.org/abs/2609.30196),
[covariance propagation](https://arxiv.org/abs/2608.15312), and
[HEPData Table 2](https://doi.org/10.17182/hepdata.156805.v1/t2).
