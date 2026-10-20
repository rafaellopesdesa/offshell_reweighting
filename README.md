# Off-shell HZZ reweighting

Notebooks for the off-shell HZZ reweighting project. Input data and generated inspection reports live in `My Drive/Colab Notebooks/HZZ_OffShell/`.

## Inspect the input ntuple

[Open 00_inspect_dataset.ipynb in Colab](https://colab.research.google.com/github/rafaellopesdesa/offshell_reweighting/blob/main/00_inspect_dataset.ipynb)

1. Choose **Runtime → Run all** and authorize the Google Drive mount. A CPU runtime is sufficient.
2. The notebook reads `dataset_wDR_mean_origIndex.h5` in place, prints its schema, attributes, column labels, sample entries, and sampled diagnostics.
3. It saves `dataset_wDR_mean_origIndex.inspection.txt` and `dataset_wDR_mean_origIndex.inspection.json` beside the ntuple. Tell the assistant when these are ready, so it can read and interpret them directly from Drive.

Shapes and stored metadata are exact. Statistics use small deterministic samples, retain signed values, and do not establish full-file event counts, yields, or weight normalization. Some variable-length object payloads are deliberately left unread; the report identifies skipped samples and errors. The HDF5 input is opened read-only.

The notebook is self-contained and uses `h5py`, `numpy`, and `pandas`, which are available in standard Colab runtimes. For local use, install those packages plus `ipython`/Jupyter and change `FILE_PATH`.

Keep datasets, reports, and notebook outputs containing event data out of this public repository.
