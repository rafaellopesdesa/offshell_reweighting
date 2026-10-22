# Off-shell HZZ demonstration

[Open notebook 06 in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/demo/colab-workspace-scans/06_root_workspace_scans.ipynb)

Use a fresh **CPU runtime** and run the cells in order. Notebook 06 installs
ROOT 6.38.00 with RooFit (`root==0.1a10`), compiles the custom class, reads
the three preserved 4096-event workspaces, fits and scans them, and compares
with HEPData's expected four-lepton NSBI curve. No previous notebooks need
to be run first.

Upload these three files to **MyDrive/HZZ_OffShell/workspaces/06_root_scans/**:

| File | Constrained NPs | Free qqZZ factors |
| --- | ---: | ---: |
| `HZZ_4l_true_asimov_4096.root` | 127 | 3 |
| `HZZ_4l_data_like_asimov_4096.root` | 127 | 3 |
| `HZZ_4l_data_like_asimov_4096_pruned54.root` | 54 | 3 |

The notebook verifies the [workspace checksums](workspaces.json), copies the
inputs to local runtime storage and opens them read-only. No NumPy sidecars
are needed. Change `DRIVE_BASE` in the notebook to use another Drive folder.

The default continues the previous study's 42-point **[0, 2]** scan, including
0.01. For a shorter demonstration, select the pruned workspace and/or use
`GRID = [0.0, 1.0, 2.0]`. A complete run of all three can take several hours;
the notebook saves every completed point and resumes if interrupted.
It runs one workspace per subprocess to release memory between models.
Keep the same code/ROOT version when resuming, or choose a new output folder.

Results go to **MyDrive/HZZ_OffShell/likelihoods/06_root_scans/full_normalization_v1/root_6_38_00/**:

- `asimov_vs_hepdata.png` and `.pdf`: curves and residuals.
- `fit_summary.csv`: fitted mu, number of NPs, comparison and measured costs.
- Per-workspace scan and HEPData comparison CSVs.
- Per-workspace `scan.json`: resumable global/profile fits, all fitted
  parameter values, convergence and independent likelihood checks.

`scan_workspace.py` contains the RooFit likelihood and minimization. It uses
the weighted extended Asimov NLL without a SumW2 correction, preserves all
fixed/global observables, and profiles all three qqZZ factors. Each curve is
relative to its continuous global minimum. Neighboring points reuse Minuit
state within each side of the minimum. Failed fits stop the calculation;
the HEPData values are never re-zeroed or rescaled.

All three models use `RooSameSampleNormalizer` to normalize the interpolated
shape at every nuisance point, keeping physical
rate variations separate. All three generate at mu=1; the data-like inputs
shift the nuisance values, including qqZZ, CR expectations and Gaussian
auxiliaries, to the observed-data best fit. Pruning fixes 73 alphas at these
generating values. The true and data-like unpruned models now differ only in
their generating nuisance values and corresponding Asimov observations.

Use the new `HZZ_4l_true_asimov_4096.root`, not the historical
`HZZ_4l_asimov_4096.root` with anchor-only normalization. Its old scan belongs
to the old model and is preserved separately. The new output folder keeps the
updated demonstration separate from those earlier results.

The matching three **8192-event workspaces** are archived for future use;
[workspaces_8192.json](workspaces_8192.json) records their names and checksums.
They use the same class, observed nuisance injection and 54-NP selection.
The 8192 reference includes two signed negative weights, retained without
clipping. They have construction/likelihood validation, with no new fits or
scans. Notebook 06 continues to select the three 4096 files only.

See [normalizer/](normalizer/) for the standalone shared-library build and
[reference/](reference/) for the exact HEPData series and provenance.

## Earlier notebooks

Copies of the original notebooks and their helpers are included. Their
Colab helper downloads point to this branch; the calculations are unchanged.
Saved outputs and execution counters are cleared for a fresh demonstration.
They use their original input dataset in Drive, independently of notebook 06.

| Notebook | Purpose |
| --- | --- |
| [00](00_inspect_dataset.ipynb) | Inspect the HDF5 dataset |
| [01](01_reference_sample.ipynb) | Construct and inspect the reference sample |
| [02](02_reweighting_tests.ipynb) | Test reweighting against direct samples |
| [03](03_same_sample_normalization.ipynb) | Apply same-sample normalization |
| [04](04_likelihood_scans.ipynb) | Stat-only likelihood scans |
| [05](05_full_model.ipynb) | Full nominal model with free qqZZ factors |
| [06](06_root_workspace_scans.ipynb) | Native RooFit workspaces with systematics |

ROOT files, generated results and build products are excluded from this branch.
The analysis and presentation branches retain the full studies and results.

## Validation

The standalone library builds with ROOT 6.34.10 and with the notebook's pinned
ROOT 6.38.00 wheel. The original demonstration's global fits, representative
profiles, resume and exports were tested on Unity; its 126 saved profile vectors
agree across ROOT versions within 1.5e-11. Those old true-model profiles validate
the historical model, not the replacement. The replacement and the three 8192
files pass serialization, normalization, generating-score and independent NLL
checks; no replacement-model fit/scan is claimed. All six files are also checked
with the notebook's ROOT wheel. Mount and authorize Drive in your own Colab session.
