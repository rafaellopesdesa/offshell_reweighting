# Off-shell workspace results slideshow

A static Reveal.js 6.0.2 presentation with 35 slides. After the overview,
six consecutive slides compare the current true, data-like and pruned Asimov
scans at 4096 and 8192 events with HEPData Table 2 Exp NSBI (4l, stat+syst).
Each comparison uses all 61 coordinates from 0 through 3, including mu=1,
with common axes, an signed difference panel, and original CSV/comparison
CSV/PDF downloads. All current models use RooSameSampleNormalizer.

Earlier observed, anchor-normalized, data-like, covariance and pruning studies
follow with their original figures and conventions. They remain identified as
earlier studies. The page displays its frozen snapshot time and does not query
Unity or update automatically. No event arrays or ROOT workspaces are published.

## Refresh the current six-model comparisons

After all six scan validations and the CSV bundle complete, run inside an
allocation using the analysis report generator and presentation exporter:

```bash
python /path/to/analysis/report_family_scan_comparisons.py \
  --run /path/to/results/asimov_scans_0_3_61 \
  --hepdata /path/to/results/hepdata_comparison/table_2.json \
  --output /path/to/results/asimov_scans_0_3_61/presentation_report
python presentation/tools/export_family_results.py \
  --run /path/to/results/asimov_scans_0_3_61 \
  --report /path/to/results/asimov_scans_0_3_61/presentation_report
python presentation/tools/check_presentation.py --output /path/to/browser_checks
```

The generator cross-checks the cached official HEPData JSON and CSV and uses
all 61 exact shared coordinates. It preserves the published values and the
continuous fitted likelihood baseline: no shifts, rescaling, clipping or
smoothing. Every figure uses the same x, q and residual limits. The exporter
requires all six complete validated scans and verifies the original scan CSV
hashes. `data/family_snapshot.json` records numerical results and provenance;
`data/family_scans_csv.zip` contains the complete original scan/parameter tables.

The current true files replace the old anchor-only true models in the current
comparison. Their historical results remain on the later slides. Both pruned
models retain the same 54 ranked constrained NPs and all three qqZZ factors.
The 8192 reference sample preserves its two negative quadrature weights.

## Enable GitHub Pages

In the repository's **Settings → Pages**, select **Deploy from a branch**,
choose **codex/asimov-results-slides-20261021**, and select **/ (root)**.
The root `index.html` redirects to `presentation/`; `.nojekyll` enables static
serving. No build workflow, npm install, or external runtime CDN is required.

The expected URL, once Pages is enabled, is:
`https://rafaellopesdesa.github.io/offshell_reweighting/`.
The branch does not enable Pages or change repository settings itself.

## View and navigate

Open `presentation/index.html` directly, or serve the repository root with
`python3 -m http.server 8000` and visit `/presentation/`.
Use arrow keys / Space to move, Esc for overview, F for full screen, and the
Contents button to jump to a section. Append `?print-pdf` for browser PDF export.
On small screens, the same material appears as a readable scrolling page.
The last slide links the vector figures, CSV tables, fitted parameters, and
`data/snapshot.json` with relative source paths and SHA-256 hashes.

## Refresh the numerical snapshot

With Python and matplotlib installed, run this on Unity inside an allocation:

```bash
python3 presentation/tools/export_results.py \
  --results /path/to/work-area/results \
  --output presentation
```

This reads saved JSON results only. It never opens ROOT workspaces or event
arrays, changes inputs, or launches fits. It reads the corrected experiment
under `asimov_anchor_normalized/n_SIZE/{workspace,fit}` and its completed
HEPData comparison. Both scans must contain all 31 coordinates, pass the
independent likelihood checks, and agree at overlaps before export. The
aggregate CSVs and figures reproduce the saved continuous-fit baseline.

The discussion includes the original CR and nuisance-score diagnosis, the
physical-rate-preserving anchor correction, and the before/after minima.
The old 4096 estimate is labeled as a known higher local minimum, not a valid
global baseline. Historical first-start integral discrepancies are explicitly
separated from the new selected-fit diagnostics. The new global fits saved
MIGRAD covariance estimates; their parameter errors are not full Hesse errors.

Review the regenerated tables, figures, and text before committing an update.

Reveal.js is vendored with its MIT license and source/checksum manifest under
`vendor/reveal/`. It is unmodified. See the
[official documentation](https://revealjs.com/installation/).

To validate the rendered slides after export, use Playwright with Chromium
inside an allocation:

```bash
python3 presentation/tools/check_presentation.py --output /work/path/browser_checks
```

This checks all desktop slides, local downloads, JavaScript errors, navigation,
mobile overflow, and hosting under a GitHub project subpath. It creates a
temporary localhost server and writes screenshots plus a JSON report.


## Data-like 4096 follow-up and source decomposition

The new snapshot adds the data-like injection, full nuisance-dependent
same-reference normalization, both 42-point scans through mu=2, controlled
normalizer timing, and covariance-propagation tables/ranking. The original
4096/8192 results remain as historical true-Asimov slides.

After all numerical and accounting checks pass, run:

```bash
python presentation/tools/export_shifted_results.py \
  --results /path/to/results/asimov_data_like_4096 \
  --analysis-commit ANALYSIS_COMMIT
python presentation/tools/export_decomposition.py \
  --results /path/to/results/asimov_data_like_4096/decomposition \
  --analysis-commit ANALYSIS_COMMIT
python presentation/tools/check_presentation.py --output /path/to/browser_checks
```

Both exporters require complete validated results. `decomposition_snapshot`
records observed-Hessian and propagated totals separately, with no forced
closure. Its plot ranks unit auxiliary-observation responses, not conventional
conditional nuisance impacts. Covariance CSVs and all named-nuisance diagnostics
are included in `presentation/data/`; no ROOT workspaces or event arrays are
published.


The ranking slide links [all 130 NPs in one sorted CSV](data/all_np_uncertainty_ranking.csv),
with true/data-like source responses side by side and per-model ranks.
[Column definitions](data/all_np_uncertainty_ranking_README.md) explain the
importance metric and the unconstrained qqZZ rows. The analysis exporter
`export_np_ranking.py` creates these files from the saved decomposition summary;
`export_decomposition.py` copies them into the presentation.


## Pruned data-like follow-up

The pruned model keeps the first 54 constrained NPs in the shared ranking and
all three qqZZ normalizations. The other 73 coordinates are fixed at their
injected values; the data are unchanged and all older workspaces are retained.
Four slides show the restriction, its HEPData comparison, stage accounting and
per-mu minimization time/process RSS. Each coordinate uses a fresh Slurm process.

After the complete scan and independent checks pass, export its report with:

```bash
python presentation/tools/export_pruned_results.py \
  --report /path/to/results/asimov_data_like_4096_pruned54/report \
  --analysis-commit ANALYSIS_COMMIT
```
