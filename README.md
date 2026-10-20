# Off-shell HZZ reweighting

Notebooks for the off-shell HZZ reweighting project. Input data and generated inspection reports live in `My Drive/Colab Notebooks/HZZ_OffShell/`.

## Inspect the input ntuple

[Open 00_inspect_dataset.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/00_inspect_dataset.ipynb)

1. Choose **Runtime → Run all** and authorize the Google Drive mount. A CPU runtime is sufficient.
2. The notebook reads `dataset_wDR_mean_origIndex.h5` in place, prints its schema, attributes, column labels, sample entries, and sampled diagnostics.
3. It saves `dataset_wDR_mean_origIndex.inspection.txt` and `dataset_wDR_mean_origIndex.inspection.json` beside the ntuple. Tell the assistant when these are ready, so it can read and interpret them directly from Drive.

Shapes and stored metadata are exact. Statistics use small deterministic samples, retain signed values, and do not establish full-file event counts, yields, or weight normalization. Some variable-length object payloads are deliberately left unread; the report identifies skipped samples and errors. The HDF5 input is opened read-only.

The notebook is self-contained and uses `h5py`, `numpy`, and `pandas`, which are available in standard Colab runtimes. For local use, install those packages plus `ipython`/Jupyter and change `FILE_PATH`.

## Compare the reference and target samples

[Open 01_reference_sample.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/01_reference_sample.ipynb)

Run all cells on a CPU runtime. Notebook 01 uses `weight_kFact` to compare the reference (`Type` in `S`, `EWSBI10`) with the target (`Type` in `SBI`, `qqZZ`) for `NN_disc_New` and `m4l_fsr`.

Each variable produces two separate figures: **weighted yields first, then a normalized shape**. Bins are fixed to `[-1, 1]` in steps of `0.05` for `NN_disc_New`, and `[150, 2000]` GeV in steps of `50` GeV for `m4l_fsr`. Signed weights are retained, and rows with nonfinite plotted values or weights are counted and excluded. Yields outside the plotting ranges are reported separately; normalized shapes use the in-range signed weight.

Figure names end in `_yields` or `_normalized`. ATLAS styling uses mplhep; the label uses public Matplotlib text calls for compatibility with Colab's Matplotlib versions. Re-running setup downloads and reloads the updated helper.

The four-column reader and ATLAS plot formatting are in `utils.py`; sample selections and NumPy histograms remain in the notebook. PNG/PDF figures are saved to `HZZ_OffShell/plots/01_reference_sample/` in Drive. The reader is tailored to this ntuple's pandas fixed-format HDF5 layout.

Keep datasets, reports, and notebook outputs containing event data out of this public repository.
