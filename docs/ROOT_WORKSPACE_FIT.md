# ROOT workspace fits and Asimov validation

The corrected experiment is documented in
[Anchor-normalized results](ANCHOR_NORMALIZED_RESULTS.md), with downloadable
aggregate tables under [reports/anchor_normalized](../reports/anchor_normalized/).

After the corrected profiles, independent validation, and comparison finish,
collect measured Slurm execution time and memory, then export the aggregate
results (including parameter values, scans, validation and resource CSVs):

```bash
python3 collect_anchor_resources.py --results /work/path/results/asimov_anchor_normalized
python3 report_anchor_scans.py --results /work/path/results/asimov_anchor_normalized \
  --output reports/anchor_normalized
```

The collector requires completed jobs and final peak-memory accounting.
It saves raw Slurm accounting and its checksum. Per-phase elapsed time spans
the first start to the last finish; CPU time sums parent-job TotalCPU without
counting job steps again. Memory is the largest individual job's peak RSS,
not summed concurrent memory. Every job requests one CPU and 4 GiB, with
four profile chunks per sample running in parallel. The global stage includes
construction, both starts, and checks; finalization includes numerical scan
validation. Measurements come from different nodes and are not a controlled
hardware benchmark. See [sacct field definitions](https://slurm.schedmd.com/sacct.html).

## Anchor-normalized profile follow-up (2026-10-22)

After both new global fits finish, submit their profiles with:

```bash
python3 submit_anchor_scans.py --results /work/path/results/asimov_anchor_normalized \
  --reference /work/path/results/hepdata_comparison --setup /work/path/setup_root.sh
```

The submitter requires both original global jobs to be COMPLETED with exit 0
and independent global validation to pass. It does not interrupt them or
resume the canceled nominal-only scans. It requires committed source and
records its commit, script checksums, exact commands, and job dependencies.
Four overlapping chunks per sample cover 0.01, 0.05, ..., 1.50 with all 130
nuisance parameters profiled. Each requests one CPU and 4 GB; 4096 chunks
request two hours and 8192 chunks three hours. Dependent jobs merge all 31
coordinates, check overlaps, validate every native likelihood independently,
and compare with the official Table 2 Exp NSBI values at 30 shared coordinates.

New fit directories are `n_SIZE/fit/`; construction outputs are in
`n_SIZE/workspace/`. The comparison tool accepts `--sample-subdirectory fit`
and `--normalization-mode anchors`, and requires successful final profile
validation. Anchor correction preserves supplied physical rate factors. It
does not impose normalization throughout nonlinear nuisance interpolation.
The unconditional results contain MIGRAD covariance estimates, with no
explicit full Hesse calculation; do not describe their errors as Hesse errors.

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
file read-only and generated ROOT files, event arrays, and full job outputs
outside this public repository; publish only the reviewed aggregate reports.
The command-line driver requires `SLURM_JOB_ID` for fitting and scanning.

For the two-sample Asimov workflow, `report_asimov_fits.py --fits RESULTS
--status-doc PROJECT/docs/CURRENT_STATUS.md` records the completed numerical
comparison in the project status document. Run it after provenance succeeds;
it requires the completed run summary and unchanged-input check. It preserves
earlier entries and explicitly leaves visual plot inspection for a manual check.
The comparison also saves a separate plot for each sample so that a large
discrepancy in one curve does not obscure the other on the combined plot.

If a physical profile point improves the saved global minimum, the scan
raises an error rather than accepting a negative likelihood-ratio statistic.
Independently validate that point, then seed a new unconditional fit with
`--initial-parameters seed.json --starts MU`. The seed must contain the exact
input SHA-256 and a `parameters` object mapping all 131 names to values.
Use a separate output directory to preserve the earlier local minimum.
A successful covariance or agreement between two starts does not establish
global optimality; all profiles must still be checked against the new minimum.

`diagnose_asimov_stationarity.py` checks whether the nominal generation point
is stationary in a serialized weighted-reference Asimov workspace. It reads
both combined and separate control-region counts, compares them with nominal
model yields, checks the fixed global observables, evaluates native NLL
derivatives for all 131 parameters, and independently checks the largest
systematic derivatives against the SR integral mismatch. It also fits mu and
the three qqZZ factors from two starts while fixing all systematic parameters
to zero. Run it inside a Slurm allocation; it writes JSON only:

```bash
python3 diagnose_asimov_stationarity.py \
  --input /work/path/to/HZZ_4l_asimov_8192.root \
  --output /work/path/to/nominal_stationarity.json
```

For reference weights `u_i`, intensity `I_i(theta)`, expected SR yield
`nu(theta)`, and Asimov weights `a_i = u_i I_i(theta0)`, the SR NLL score at
generation is `d nu/d theta - sum_i u_i d I_i/d theta`. Asimov Poisson counts
and zero Gaussian pulls contribute zero score at that point. Nominal
same-sample normalization closes the mu and qqZZ directions with systematic
parameters fixed. Keeping the supplied systematic factors unnormalized on
the subset does not close their derivatives. Consequently this construction
does not enforce an exact stationary nominal point when all systematics are
profiled. The diagnostic does not renormalize these variations or modify the
model. A future change that enforces nuisance-dependent integral closure
would be a different construction prescription; stationarity alone would
still not certify a global minimum, particularly with signed weights.

For the separately authorized anchor-normalization experiment, pass
`--normalize-systematic-anchors` to `build_root_asimov.py`. For each subset,
SR process (including each qqZZ jet bin), nuisance and direction, it replaces
the supplied mixed shape/rate factor `v` by `v * R / sum(u*rhat*v)`, with `R`
from `list_weight_up/down_<process>.npy`. The correction preserves the intended
varied rate. It leaves the nominal Asimov weights, CRs, full-MC rate factors,
interpolation codes, and fit bounds unchanged. New workspaces must use a
separate output directory. `anchor_normalization.json` records all 2540
correction factors and their original integrals and target rates.

`validate_root_asimov.py` independently recomputes the corrections, checks
all serialized anchors and 258 native-versus-explicit likelihood settings,
and records remaining nominal nuisance scores. This checks anchor closure;
it does not assume closure between anchors or when several nuisances vary.

`run_anchor_asimov.slurm` builds, validates and fits one subset from mu=1 and
1.4, floating all 131 parameters. It uses Migrad strategy 0 with
`--fit-eps 0.0001`, without a full Hesse calculation or any profile scan.
`validate_anchor_global_fit.py` then checks both fitted parameter vectors
against explicit intensities, verifies saved artifacts, and compares
unchanged nominal inputs with the earlier workspace. Report the Migrad
covariance as an estimate, not a validated full-Hesse uncertainty.
After both jobs finish, `report_anchor_global_fits.py --results OUTPUT`
requires both independent validations and saves `summary.json` and
`global_fit_comparison.csv`, including mu displacement, starting-point
agreement, maximum nuisance pulls, and numerical validation precision.

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

Keep input event data, fitted workspaces, and full job outputs outside the
public repository. Reviewed aggregate figures, tables, and provenance can
be exported for the results report and slideshow.

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

## Fit the two Asimov workspaces and compare expected scans

Use `--weighted-asimov` explicitly. This mode requires the construction
metadata and generation snapshot and preserves signed quadrature weights.
Without that flag the observed-data guard still rejects non-Poisson weights.
Both fits use the same native likelihood, parameter bounds, two independent
starts, and full covariance calculation as the observed baseline.
For these larger Asimov fits use `--fit-strategy 0 --hesse`: strategy 0
avoids computing a full covariance at each preliminary start, and the final
explicit Hesse uses strategy 1. The EDM/status checks remain unchanged, and
the final covariance must have quality 3. The default global strategy remains
1 for compatibility with the observed baseline. Each completed starting-point
fit also writes its full parameter vector as a recovery checkpoint.

For the requested [0,1.5] plot, use **0.01 as the proxy for zero** and then
the HEPData coordinates 0.05, 0.10, ..., 1.50. Never evaluate a scan at zero.
An explicit `--grid` separates evaluated coordinates from the plot range:

```bash
grid=(0.01 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50
      0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
      1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50)

python3 root_workspace_fit.py fit --weighted-asimov --hesse --fit-strategy 0 \
  --input /work/path/to/asimov/n_4096/HZZ_4l_asimov_4096.root \
  --output /work/path/to/asimov_fits/n_4096

python3 root_workspace_fit.py scan --weighted-asimov \
  --input /work/path/to/asimov/n_4096/HZZ_4l_asimov_4096.root \
  --output /work/path/to/asimov_fits/n_4096 \
  --range 0 1.5 --grid "${grid[@]}"

python3 root_workspace_fit.py finalize --weighted-asimov \
  --input /work/path/to/asimov/n_4096/HZZ_4l_asimov_4096.root \
  --output /work/path/to/asimov_fits/n_4096 \
  --range 0 1.5 --grid "${grid[@]}"

python3 root_workspace_validation.py --weighted-asimov --range 0.01 1.5 \
  --input /work/path/to/asimov/n_4096/HZZ_4l_asimov_4096.root \
  --bestfit /work/path/to/asimov_fits/n_4096/bestfit.json \
  --scan-points /work/path/to/asimov_fits/n_4096/scan_points.json \
  --output /work/path/to/asimov_fits/n_4096/likelihood_validation.json
```

Repeat for 8192. For longer work, `run_root_fit.slurm`,
`run_root_scan.slurm`, and `run_root_finalize.slurm` accept the ROOT setup
script, repository directory, input workspace, output directory, and their
remaining documented arguments. Supply the current Unity account, partition,
CPU count, memory, walltime, and dependencies to `sbatch`; these are not
hard-coded in the scripts. Each process uses one CPU. Overlapping chunks
allow independent convergence checks before merging. The finalizer expects
chunk directories named `scan_*` under that sample's output directory.

`release_root_scans.py --run /work/path/to/asimov_fits` is an optional
lightweight Slurm controller for runs with saved `fit_jobs.json`,
`scan_jobs.json`, and `finalize_jobs.json` submission manifests. When both
global starts have produced `bestfit_checkpoint.json`, it replaces only held,
pending scan jobs with jobs using that checkpoint. It leaves fits and
already-active scans untouched. Finalization depends on both the completed fit and all new
scan jobs, so the final covariance and minimum still govern the result.
Use a small workflow allocation for this controller; its polling interval is
30 seconds. Submission changes are recorded in `scan_release_history.json`.

The validator independently reconstructs every saved profile-point likelihood
and verifies the saved covariance and snapshot. It records continuous
intensity minima at fixed profile nuisance settings across the requested
interval, while requiring positivity at each actual profile point. Positivity
at one profile setting does not establish positivity at every nuisance value.
For profile validation, the global minimum and the saved profile points are
the required physical configurations. A fixed fitted nuisance vector can
become unphysical at another mu, which is reported as a domain diagnostic;
the fit must reprofile the nuisances there. The six nominal/perturbed
fixed-setting checks remain available when no `--scan-points` is supplied.

If a profile warm start has a non-finite NLL after changing mu, the driver
tries nuisance vectors at 1/2, 1/4, and 1/8 of the displacement from the
nominal vector, then nominal itself. It resets the numerical steps and
minimizer covariance after such a fallback and profiles all 130 parameters
normally. The saved `warm_start_fraction` records this choice. This only
chooses a physical starting point: no density, bound, or likelihood term is
clipped or changed.

When a weighted-Asimov input has its construction-side
`reference_weights_normalized.npy` alongside it, validation also compares
the finite-reference integral `sum(u_i * intensity_i)` with the model's
full-MC expected SR yield at every checked setting. The sidecar path and
checksum are recorded. Nominal same-sample normalization closes this integral
exactly; supplied systematic variations can cause a mismatch at fitted
nuisance values. This is a diagnostic of the quadrature, not a normalization
correction or an additional likelihood term. Signed reference weights are
retained.

```bash
python3 compare_asimov_hepdata.py \
  --fits /work/path/to/asimov_fits \
  --reference /work/path/to/hepdata/table_2.json \
  --crosscheck-table /work/path/to/hepdata/table_1.json \
  --reference-csv /work/path/to/hepdata/table_2.csv \
  --output /work/path/to/asimov_fits/hepdata_comparison
```

The comparison selects **Exp NSBI**, the four-lepton stat+syst series in
Table 2 / Figure 11(b). It verifies the duplicate table and CSV, writes both
local curves and their signed residuals in PNG/PDF/ROOT, and saves per-sample
and between-sample CSVs plus a checksum-bearing JSON summary. Only the 30
shared coordinates from 0.05 to 1.50 enter reference residuals. The published
point at zero remains at zero, separate from the local 0.01 proxy. Neither
curve is shifted, interpolated, or rescaled. The statistics-only series uses
a different nuisance prescription and is not the counterpart of these fits.

Finite-sample nuisance score closure is not guaranteed by nominal-only
same-sample normalization. Preserve any fitted shifts and report them. The
provided NumPy yields/systematics also differ from the observed ROOT template;
these two nested subsets alone cannot separate input-model differences from
quadrature effects or estimate a distribution of finite-sample fluctuations.

The `root_workspace_fit.py` driver requires explicit `--weighted-asimov`
for these signed, weighted datasets; the observed-data guard remains active
without that flag. Keep event arrays, ROOT files, and full job manifests
outside the public repository. The validated aggregate reports under
`reports/anchor_normalized/` and the presentation's figures and tables are
the reviewable public results requested for this experiment.
