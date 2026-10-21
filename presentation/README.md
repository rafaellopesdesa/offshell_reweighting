# Off-shell workspace results slideshow

A static Reveal.js 6.0.2 presentation of the observed-workspace validation and
the available 4096/8192 reference-sample Asimov results. The page displays its
snapshot time and labels incomplete scans as preliminary. It does not query
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
arrays, changes inputs, or launches fits. Partial profiles use markers only;
missing points are not interpolated. Profiles are rebased to the chosen saved
global minimum. Final global fits are used when available; otherwise the
two-start checkpoint is used. Both full profile coverage and successful final
validation are required before the snapshot is labeled complete.

Review the regenerated tables, figures, and text before committing an update.
The first-start normalization diagnostic is deliberately preserved and labeled
as such; it is not silently replaced with a later fit. Error columns in a
checkpoint parameter CSV are provisional until the full covariance is saved.

Reveal.js is vendored with its MIT license and source/checksum manifest under
`vendor/reveal/`. It is unmodified. See the
[official documentation](https://revealjs.com/installation/).
