# Fit the observed ROOT workspace

`root_workspace_fit.py` reads an existing ROOT file and fits its native
extended likelihood with RooFit/Minuit2. It does not rebuild the model or
renormalize the supplied density ratios. ROOT >= 6.34.04 with RooFit,
RooStats, and PyROOT is required. The tested runtime is ROOT 6.34.10 with
Python 3.12.3 on Ubuntu 24.04 x86_64.

`verify_root_runtime.py --output-dir /work/path/to/runtime_check` provides a
small seeded Gaussian fit and workspace save/load check for a new ROOT
installation, without using the observed dataset.

Run the numerical work inside a Slurm allocation on Unity. Activate the
project ROOT environment before running the commands below. Keep the input
file read-only and all generated outputs outside this public repository.
The command-line driver requires `SLURM_JOB_ID` for fitting and scanning.

```bash
python3 root_workspace_fit.py fit \
  --input /path/to/read-only/HZZ_4l_offshell_v2.root \
  --output /work/path/to/data_fit --hesse

python3 root_workspace_fit.py scan \
  --input /path/to/read-only/HZZ_4l_offshell_v2.root \
  --output /work/path/to/data_fit --range 0.5 1.5 --points 41

python3 root_workspace_fit.py finalize \
  --input /path/to/read-only/HZZ_4l_offshell_v2.root \
  --output /work/path/to/data_fit --range 0.5 1.5

python3 root_workspace_validation.py \
  --input /path/to/read-only/HZZ_4l_offshell_v2.root \
  --bestfit /work/path/to/data_fit/bestfit.json \
  --output /work/path/to/data_fit/likelihood_validation.json
```

The default workspace, ModelConfig, and dataset are `combined`,
`ModelConfig_4l_offshell`, and `obsData_4l_offshell`; each can be overridden
on the command line. ROOT's CPU evaluation backend is used: the code-generation
backend in ROOT 6.34.10 does not support the workspace's `RooWrapperPdf` nodes.

For two-CPU allocations, the scan can start from `--bestfit
/work/path/to/data_fit/bestfit_checkpoint.json` while the final covariance
calculation finishes on the other CPU. After both jobs finish, `finalize`
checks that all requested points converged, recomputes likelihood ratios
against the final `bestfit.json`, and writes the final tables and plots.
Separate scan chunks can be combined with `finalize --scan-parts
/path/to/first/scan_points.json /path/to/second/scan_points.json`. The union
must cover the requested grid exactly. Overlapping points must agree within
1e-4 in `2*NLL`; the lower NLL is retained. Never write two running scans to
the same output directory.

## Statistical conventions

The unconditional fit floats `mu`, the 127 systematic `alpha_*` parameters,
and the three free `qqZZ` normalization parameters. It preserves all input
parameter bounds and fixed global observables. The same constraint PDF is
shared across channels; RooFit extracts it once. Control-region integer
count weights are preserved and no SumW2 correction is applied.

The three `qqZZ` parameters are cumulative in this ROOT workspace. Effective
jet-bin scales are `mu_qqZZ`, `mu_qqZZ * mu_qqZZ_1`, and
`mu_qqZZ * mu_qqZZ_1 * mu_qqZZ_2`. The best-fit JSON saves these derived scales
as well as the original parameters. This differs from the independent-factor
parameterization used in notebook 05.

Global fits start at `mu=1` and `mu=1.4`, with the other parameters reset to
their input values, and use Minuit2/Migrad strategy 1. `--hesse` additionally
calculates the full covariance matrix. The profile fixes only `mu`, leaving
130 parameters free, and walks outward from the global minimum using warm
starts. It uses strategy 0, with a strategy-1 retry on failure. All minimizations
must have status zero, finite NLL, and EDM <= 1e-4. The requested convergence
tolerance is 0.01 (Minuit's observed EDM target is 1e-5 for this NLL).

The vertical axis is `2 * (NLL_profile(mu) - NLL_global_min)`. The unconditional
minimum uses the original parameter range, not the plot interval. A restricted
scan is not sufficient to report confidence endpoints outside that interval.
The symmetric Hesse error, if quoted, is a local curvature error and is not a
substitute for an asymmetric profile-likelihood interval.

## Saved outputs

- `bestfit.json`: full precision `mu_hat`, all 127 `alpha_hat` values, all 131
  floating parameters with errors and bounds, effective `qqZZ` scales,
  convergence diagnostics, ROOT version, and input SHA-256.
- `bestfit_parameters.csv`: one row per floating parameter.
- `bestfit.root`: the `RooFitResult`, parameter snapshot, covariance and
  correlation matrices.
- `postfit_workspace.root`: a separate workspace file containing the `bestfit`
  snapshot; the original input is never opened for writing.
- `scan_points.json`: every profile point and its complete fitted parameter
  vector, written after each minimization. `--resume` continues an interrupted
  scan against the same baseline.
- `likelihood_scan.csv`, `scan_summary.json`, and PNG/PDF/ROOT plots.

The validation script assembles the extended likelihood explicitly from the
stored event-intensity functions, control-region yields, and Gaussian penalties.
It compares likelihood *differences*, allowing the parameter-independent
constants in RooFit's NLL. It also checks continuous positivity on [0.5,1.5]
using the quadratic dependence on `sqrt(mu)`, at the tested fixed nuisance
settings. It does not claim positivity everywhere in the nuisance-parameter
domain; RooFit must reject invalid minimizer trial points.

Do not commit input data, fitted workspaces, or generated result files to this
public repository.

## Compare the observed scan with HEPData

The matching published result is the `Obs NSBI` series in
[HEPData 156805 v1, Table 2](https://doi.org/10.17182/hepdata.156805.v1/t2),
corresponding to Figure 11(b): four-lepton only, including statistical and
systematic uncertainties. Table 1 / Figure 11(a) repeats this same series.
The other tables include different observables or channel combinations.
Equation (16) of [the paper](https://arxiv.org/html/2412.01548v2#S7) defines
the same statistic as the saved RooFit scan. The published confidence belts
use a Neyman construction; `q=1` is not its 68% confidence belt.

Download the public table JSON and CSV into an output directory outside the
repository, then run the comparison with the ROOT environment activated:

```bash
curl --fail --location 'https://www.hepdata.net/record/data/156805/1759753/1' \
  --output /work/path/to/comparison/table_2.json
curl --fail --location 'https://www.hepdata.net/record/data/156805/1759752/1' \
  --output /work/path/to/comparison/table_1.json

python3 compare_hepdata_scan.py \
  --scan /work/path/to/data_fit/likelihood_scan.csv \
  --bestfit /work/path/to/data_fit/bestfit.json \
  --reference /work/path/to/comparison/table_2.json \
  --crosscheck-table /work/path/to/comparison/table_1.json \
  --output /work/path/to/comparison
```

Optionally pass the official Table 2 CSV download as `--reference-csv` to
verify that its `Obs NSBI` values match the JSON. The script compares only
shared mu coordinates, with no interpolation, extrapolation, vertical shift,
or fitted rescaling. It writes a point-by-point CSV, a JSON summary with
input checksums and absolute residual metrics, and PNG/PDF/ROOT plots.
The lower panel shows the signed difference, local minus HEPData. Its gray
band is the half-step of three-decimal rounding, not a statistical error.
Relative residuals are avoided because the statistic vanishes at its minimum.

## Build reference-sample Asimov workspaces

`root_asimov_inputs.py`, `build_root_asimov.py`, and `validate_root_asimov.py`
construct and check ROOT workspaces from the supplied reference-sample NumPy
arrays. Run all three inside a Slurm allocation with the ROOT environment
activated. Staging opens inputs with `mmap_mode="r"`, extracts the selected
systematic slices once, and saves compact arrays under `/work`. It does not
read any `_full` arrays or modify the source directory.

```bash
python3 root_asimov_inputs.py \
  --source /project/pi_rclsa_umass_edu/Inference_data_RefRecreationTest \
  --output /work/path/to/asimov/staged_reference --maximum 8192 --seed 42 --workers 2

python3 build_root_asimov.py \
  --template /path/to/read-only/HZZ_4l_offshell_v2.root \
  --staged /work/path/to/asimov/staged_reference \
  --size 4096 --output /work/path/to/asimov/n_4096

python3 validate_root_asimov.py \
  --input /work/path/to/asimov/n_4096/HZZ_4l_asimov_4096.root \
  --staged /work/path/to/asimov/staged_reference \
  --output /work/path/to/asimov/n_4096/validation.json
```

Repeat the build and validation with size 8192 and a separate `n_8192` output
directory. The builder refuses to overwrite an existing final workspace.
Neither construction nor validation performs a fit or profile scan.
The staging worker count must fit within the allocation. Completed slices
are checkpointed and checksum-verified when resuming an interrupted staging
command with the same selection settings. The final staging manifest is
written only after every required input has been extracted successfully.

The samples are prefixes of one seeded random permutation, matching notebook
05. Each prefix independently computes signed reference weights
`u=w/sum(w)`, process factors `Z_p=sum(u*r_p)`, and corrected ratios
`rhat_p=r_p/Z_p`. The reference is the already-selected S+EWSBI10 sample in
the supplied arrays. No events are discarded or selected using a Type label.
Negative MC weights remain signed.

Native inclusive yields come from `yields_array_SR.npy`. The qqZZ jet-bin
fractions are computed separately on each subset, and their yields sum to
the fixed inclusive yield. At the nominal SM point the SR event weights are
`u * (lambda_SBI*rhat_SBI + lambda_EWSBI*rhat_EWSBI +
lambda_qqZZ*rhat_qqZZ + lambda_ttV*rhat_ttV)`. The 4096/8192 counts refer to
SR quadrature events, not physical expected yields. Each combined dataset
also contains three rows weighted by the CR MC expectations.

The ROOT model preserves the observed template's function/PDF/variable names,
classes, ModelConfig sets, parameter bounds, Gaussian constraints, and
interpolation codes. It retains cumulative qqZZ factors. The new files use
the NumPy inputs for nominal yields, event variations, and rate variations;
these inputs differ from the observed template in some numerical values.
The CR process yields use double-precision `hist_combined.npy`, checked
against the float32 `yields_array_CR.npy`. The generator records all 34 yield
replacements and the bijective mapping of 127 systematic names.

Event-level up/down arrays already include their rate changes. The template
has `nu=lambda*rate_modifier` and `shape=rhat*event_modifier/rate_modifier`,
so their product retains the supplied event variation. Full-MC rate factors
come from `list_weight_*` in the SR and `arr_ratio_*_CR` in the CRs. They are
not recalibrated on each subset. Nominal same-sample normalization therefore
guarantees nominal yield and free-parameter score closure; it does not force
every nuisance score to vanish on a finite subset.

Outputs for each size:

- `HZZ_4l_asimov_N.root`: workspace `combined`, PDF `simPdf`, ModelConfig
  `ModelConfig_4l_offshell`, and Asimov data under the compatible
  `obsData_4l_offshell` name. The three `data_CR*` datasets also contain Asimov
  counts. The file contains no observed data rows.
- `asimov_generation`: workspace snapshot of POIs, nuisance parameters, and
  global observables at the generation point.
- `construction.json` and ROOT `asimov_metadata`: input provenance, selection,
  normalization factors, fractions, yields, and model-structure checks.
- `event_rows.npy`, `reference_weights_normalized.npy`, `asimov_weights.npy`,
  `jet_category.npy`, and `nominal_ratios.npz`: exact quadrature construction.
- `validation.json`: serialization and column checks, independent weighted
  likelihood comparisons, continuous nominal positivity on [0.5,1.5], and
  nominal score closure for mu and the three qqZZ parameters.

The validator compares all 127 single-nuisance up/down endpoints at mu=1.3,
where every signal/interference basis coefficient is nonzero, plus nominal,
mu endpoint, and cumulative qqZZ settings. It checks rates separately and
compares native RooFit NLL differences against direct NumPy event intensities,
CR Poisson terms, and Gaussian penalties. Signed Asimov weights are retained
in the weighted log likelihood; no SumW2 correction is applied.

The current `root_workspace_fit.py` driver deliberately rejects non-Poisson
weighted datasets. Before profiling these Asimov workspaces, add an explicit
Asimov-data mode with the appropriate weighted-likelihood conventions; do not
silently remove the observed-data guard. Keep all generated arrays, ROOT
files, and numerical result manifests outside the public repository.
