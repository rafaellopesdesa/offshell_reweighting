# Off-shell workspace results slideshow

A static Reveal.js 6.0.2 presentation of the observed-workspace validation and
the corrected 4096/8192 reference-sample Asimov results. The 17 slides include
the control-region checks, nuisance-score diagnosis, anchor correction,
before/after fits, completed expected-scan comparisons, and measured Slurm
execution time and memory for both sizes. The page displays
its snapshot time. It does not query
Unity or update automatically. Only aggregate fit results are included.

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
