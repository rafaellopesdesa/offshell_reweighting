# Off-shell HZZ reweighting

Notebooks for the off-shell HZZ reweighting project. Input data and generated inspection reports live in `My Drive/Colab Notebooks/HZZ_OffShell/`.

## Plotting defaults

Use **0.90–1.10 for every histogram ratio panel**, including future notebooks, unless a different range is explicitly requested. This preference is recorded in `AGENTS.md` and `DEFAULT_RATIO_YLIM` in `utils.py`.

## Inspect the input ntuple

[Open 00_inspect_dataset.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/00_inspect_dataset.ipynb)

1. Choose **Runtime → Run all** and authorize the Google Drive mount. A CPU runtime is sufficient.
2. The notebook reads `dataset_wDR_mean_origIndex.h5` in place, prints its schema, attributes, column labels, sample entries, and sampled diagnostics.
3. It saves `dataset_wDR_mean_origIndex.inspection.txt` and `dataset_wDR_mean_origIndex.inspection.json` beside the ntuple. Tell the assistant when these are ready, so it can read and interpret them directly from Drive.

Shapes and stored metadata are exact. Statistics use small deterministic samples, retain signed values, and do not establish full-file event counts, yields, or weight normalization. Some variable-length object payloads are deliberately left unread; the report identifies skipped samples and errors. The HDF5 input is opened read-only.

The notebook is self-contained and uses `h5py`, `numpy`, and `pandas`, which are available in standard Colab runtimes. For local use, install those packages plus `ipython`/Jupyter and change `FILE_PATH`.

## Compare the reference and target samples

[Open 01_reference_sample.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/01_reference_sample.ipynb)

Use a CPU runtime. If you ran an older notebook version, run the installation cell once, choose **Runtime → Restart session**, then **Run all**. Setup pins `matplotlib==3.10.0` and `mplhep==0.4.1`, detects stale imports after installation, and checks PNG/PDF export before loading the ntuple.

Notebook 01 uses `weight_kFact` to compare the reference (`Type` in `S`, `EWSBI10`) with the target (`Type` in `SBI`, `qqZZ`) for `NN_disc_New` and `m4l_fsr`.

Each variable produces two separate figures: **weighted yields first, then a normalized shape**. Bins are fixed to `[-1, 1]` in steps of `0.05` for `NN_disc_New`, and `[150, 2000]` GeV in steps of `50` GeV for `m4l_fsr`. Signed weights are retained, and rows with nonfinite plotted values or weights are counted and excluded. Yields outside the plotting ranges are reported separately; normalized shapes use the in-range signed weight.

Figure names end in `_yields` or `_normalized`. ATLAS styling uses mplhep; the label uses public Matplotlib text calls for compatibility with Colab's Matplotlib versions. Re-running setup downloads and reloads the updated helper.

The four-column reader and ATLAS plot formatting are in `utils.py`; sample selections and NumPy histograms remain in the notebook. PNG/PDF figures are saved to `HZZ_OffShell/plots/01_reference_sample/` in Drive. The reader is tailored to this ntuple's pandas fixed-format HDF5 layout.

## Test the stored reweighting ratios

[Open 02_reweighting_tests.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/02_reweighting_tests.ipynb)

Notebook 02 applies `r_SBI` and `r_qqZZ` to the reference events with yield scaling:

```python
w_reweighted = weight_kFact * (
    (lambda_sbi / lambda_ref) * r_SBI
    + (lambda_qqzz / lambda_ref) * r_qqZZ
)
```

Each lambda is the signed sum of original event weights for its sample after the same finite-value selection as notebook 01 and before plotting-range cuts. Stored ratios are used without renormalization. Ratios are checked only on reference events; signed values and negative event weights are retained.

The notebook checks SBI, qqZZ, and their sum against the corresponding direct samples. For both variables it produces weighted-yield and normalized-shape plots with **reweighted / direct** ratio panels, using the binning of notebook 01. It also prints weighted ratio means, full-sample yield closure, and underflow/in-range/overflow yields. PNG/PDF figures are saved to `HZZ_OffShell/plots/02_reweighting_tests/` in Drive.

Blue error bars show the reweighted-reference MC uncertainty; the grey band shows the direct-target MC uncertainty. Variances use squared event weights with fixed yield coefficients and stored ratios; the combined weights are summed before squaring. Shape uncertainties propagate the in-range normalization. Zero-direct-yield bins are omitted from ratio panels. All ratio axes use `RATIO_YLIM = (0.90, 1.10)`. Model uncertainties and correlations between reference and direct samples are not included.

Use the same setup/restart procedure as notebook 01. The selections, reweighting formula, and NumPy histograms are visible in the notebook; ratio plotting and normalized uncertainty propagation are in `utils.py`.

## Compare same-sample normalization

[Open 03_same_sample_normalization.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/03_same_sample_normalization.ipynb)

Step 1 normalizes the reference event weights and computes `Z_SBI = sum(w_ref_normalized * r_SBI)` and `Z_qqZZ = sum(w_ref_normalized * r_qqZZ)`. It divides each ratio by its own Z and prints the weighted means before and after correction. Z is computed once on the full selected reference sample, before plotting-range cuts.

Step 2 compares the original and corrected reweightings against the direct SBI, qqZZ, and combined targets. Each figure has both reweightings above and two reweighted/direct curves below, with ratio axes fixed to **0.90–1.10**. Both weighted-yield and unit-area plots use the established bins. PNG/PDF figures go to `HZZ_OffShell/plots/03_same_sample_normalization/` in Drive.

Corrected full-sample yield closure is enforced by construction. For individual components, the unit-area shapes coincide because a constant Z cancels; for the mixture, correcting the components separately can change the shape. Signed weights are retained. MC errors propagate the same-sample normalization and shared-reference component correlations, with full bin covariance used for unit-area shapes. Yield coefficients and stored ratios are held fixed, and model uncertainty and correlations with the direct sample are not included.

## Stat-only likelihood scans

[Open 04_likelihood_scans.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/04_likelihood_scans.ipynb)

Notebook 04 uses the [NSBI LHC toolkit](https://github.com/iris-hep/nsbi-lhc-toolkit/) for an unbinned extended likelihood, fitting and scanning one physical POI, mu, and drawing the likelihood curves with ATLAS styling. It uses SM Asimov data at mu = 1 and includes both shape and event yield. The scan defaults to **0–1.8 in steps of 0.025**, also the default for future notebooks.

The process scaling is `sqrt(mu)` for SBI, `mu - sqrt(mu)` for S, `1 - sqrt(mu)` for B, and 1 for qqZZ. The pinned toolkit supports multiplicative normalization factors; `utils_likelihood.py` maps the single physical mu to the three internal factors and writes the toolkit workspace. It does not implement a likelihood. The internal factors must remain tied together when reopening a workspace.

Compare the full S + EWSBI10 reference with 2,048, 16,384, 131,072, and 1,048,576 events. Subsets are nested, uniformly drawn without replacement with seed 42. Every sample, including the full reference, recalculates all four weighted Z factors using its own normalized signed reference weights. Full native process yields stay fixed, so reducing the reference count does not reduce the expected event yield. Each sample supplies its own same-sample-normalized Asimov integration weights; no MC-statistical, training, or systematic nuisance parameters are included.

A second figure compares the full-reference reweighted Asimov with a **direct SBI + qqZZ Asimov**. The direct sample uses the original `weight_kFact` weights and ratio branches evaluated on its own events. Its model uses the same full-reference Z factors; the ratios are not renormalized using the direct sample. The fitted minimum is free to move away from mu = 1, exposing finite-MC or density-ratio closure differences. Both constructions have the same physical yield.

The earlier scan failed because a full-reference event's combined intensity crossed zero at mu approximately 1.8947. The logarithm is then undefined even though the total expected yield remains positive. Before each fit, the helper now checks positivity over the entire requested interval, including any interior minimum of the quadratic in sqrt(mu). It reports the worst event and saves `domain_check.json`. The adapter keeps the fit within that interval as well. No event weights, densities, or individual process coefficients are clipped.

Use a CPU runtime and the notebook's setup/restart procedure. Setup pins matching CPU JAX/JAXlib 0.6.2 and removes installed CUDA JAX plugins to avoid Colab plugin-version conflicts. Toolkit source is pinned to commit `fc09848fc6540fd32310faebbe9db6eea7ecd17b`. The same finite-value selection as notebooks 01–03 is retained, including events outside the histogram ranges. Workspaces, corrected ratio arrays, original row positions, Z checks, domain checks, and scan tables are saved under `HZZ_OffShell/likelihoods/04_stat_only/`. The reference-size and Asimov-comparison PNG/PDF figures go to `HZZ_OffShell/plots/04_likelihood_scans/`.

## Full model with profiled qqZZ normalizations

[Open 05_full_model.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/05_full_model.ipynb)

Notebook 05 adds EWB, EWSBI, EWSBI10, and ttV to the model. It uses the supplied electroweak interpolation coefficients, keeps ttV fixed, and profiles three independent qqZZ normalizations for 0, 1, and at least 2 jets. The four physical fit parameters are `mu`, `mu_qqZZ0`, `mu_qqZZ1`, and `mu_qqZZ2`. The nuisance parameters have no Gaussian or other constraint terms; their physical domain is nonnegative, without an upper bound. All three remain free at every fixed-mu scan point. Scans cover **0–1.8 in steps of 0.025**.

The reference is S + EWSBI10. The SM Asimov is SBI + EWSBI + qqZZ + ttV. As in 04, compare the full reference with 2,048, 16,384, 131,072, and 1,048,576 nested random events, recomputing all eight same-sample Z factors for each case while holding the eight full native yields fixed. A second figure compares the full-reference reweighted Asimov with native events from the four SM processes, retaining their original weights and the full-reference model normalization.

To retain a common inclusive qqZZ yield and ratio without triple-counting the rate, compute `f_c = sum(u * rhat_qqZZ * I_c)` on each reference set. The workspace uses category yield `lambda_qqZZ * f_c` and ratio `rhat_qqZZ * I_c / f_c`. Their product is exactly `lambda_qqZZ * rhat_qqZZ * I_c`; the profiled qqZZ rate is `lambda_qqZZ * sum(mu_qqZZc * f_c)`. The fractions sum to one and are recalculated for each reference subset. The direct Asimov uses the same full-reference Z factors **and** jet fractions. Jet masks use the jets of each evaluation event.

`utils_full_likelihood.py` writes the ten-process workspace and ties its six signal/interference normfactors to mu, exposing only the four physical parameters. The toolkit performs the likelihood evaluation and profiling; the existing helper delegates curve drawing to its plotting function. A continuous nominal-intensity check precedes each fit, invalid trial intensities are rejected without clipping, and failed profile-minimization warnings stop the scan. Curves are relative to their own grid minimum. Fitted mu and all three nuisance estimates are recorded separately.

Use the same CPU setup as 04. Workspaces, arrays, Z checks, jet fractions, domain checks, fit summaries, and scan tables go to `HZZ_OffShell/likelihoods/05_full_model/`. The two PNG/PDF comparisons go to `HZZ_OffShell/plots/05_full_model/`. Native/direct fits can differ from the SM point because this is also a density-ratio closure test.

## Inspect Type and subType labels

[Open auxA.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/auxA.ipynb)

This small diagnostic lists all Type values, all subType values, and every Type/subType combination with event counts and signed finite-weight sums. It applies no event cuts, includes missing labels, counts nonfinite weights separately, and shows exact label representations so whitespace and capitalization are visible. It checks the exact `Type == "ttV"` selection used in 05 and highlights labels containing `tt` or `top` without choosing a replacement selection. It reads only the two label columns and `weight_kFact`, and saves CSV tables under `HZZ_OffShell/inspection/auxA/`.

Keep datasets, reports, and notebook outputs containing event data out of this public repository.
