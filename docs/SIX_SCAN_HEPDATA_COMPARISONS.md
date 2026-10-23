# Individual comparisons with HEPData

`report_family_scan_comparisons.py` reads the validated six-workspace scan
outputs and compares each with Table 2, **Exp NSBI (4l, stat+syst)**, from
[HEPData 156805](https://doi.org/10.17182/hepdata.156805.v1/t2).
It does not load ROOT files, change models, or perform additional fits.

The presentation orders the models as true, data-like, and pruned data-like
4096, followed by the same three 8192 models. All six current workspaces use
`RooSameSampleNormalizer`; the older anchor-only true models are retained as
historical studies on later slides. Both pruned models retain the same 54
ranked constrained nuisances and all three free qqZZ normalizations.

## Numerical conventions

Each comparison uses all 61 exact shared coordinates, `mu = i/20` for
`i = 0,...,60`, including zero and one. The generator cross-checks the cached
official HEPData JSON and CSV. The local statistic remains
`q(mu) = 2 * (profile NLL - continuous global-fit NLL)`; neither curve is
shifted to its grid minimum, rescaled, clipped, or smoothed. Lines simply
connect the supplied points. Tiny negative local q values from minimizer
precision remain in the downloads.

Every figure has common x, q, and residual limits. Its lower panel shows the
absolute difference `delta_q = q_scan - q_HEPData`. The sidebar gives the
continuous fitted mu, the RMS of delta_q across all 61 coordinates, and its
maximum absolute value and location. These are descriptive discrepancies,
not goodness-of-fit probabilities or uncertainty estimates.

Before plotting, the generator requires complete global fits, converged
merged scans, passing independent likelihood validation at all 61 points,
matching workspace hashes, status zero, finite EDM <= 1e-4, and consistency
between each CSV q and the continuous global-fit baseline. Provenance hashes
are recorded in `summary.json`. The presentation exporter additionally
requires the complete 366-point export and verifies the original CSV and
combined archive checksums.

## Reproduce the report

Run inside a Slurm allocation with Python and matplotlib:

```bash
python report_family_scan_comparisons.py \
  --run /path/to/results/asimov_scans_0_3_61 \
  --hepdata /path/to/results/hepdata_comparison/table_2.json \
  --output /path/to/results/asimov_scans_0_3_61/presentation_report
```

The matching `table_2.csv` must be next to the JSON. The report contains
unchanged original scan CSVs, per-model comparison CSVs, SVG/PDF/PNG figures,
the common HEPData CSV, a summary CSV, and a complete JSON snapshot. The
optional `--preview` flag permits private checks of already validated
subsets; incomplete reports are rejected by the presentation exporter.

On the presentation branch, run `presentation/tools/export_family_results.py`
with `--report` and `--run` pointing to the same directories. Then run
`presentation/tools/check_presentation.py --output /path/to/browser_checks`
to check desktop slides, mobile navigation, numeric labels, download links,
and that the six comparisons occur consecutively after the overview.
