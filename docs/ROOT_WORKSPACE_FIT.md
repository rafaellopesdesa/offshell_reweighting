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
