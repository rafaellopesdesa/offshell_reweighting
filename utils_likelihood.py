"""Workspace formatting and the physical-mu adapter for notebook 04.

All likelihood evaluation, fitting, scanning, and curve drawing use
nsbi_common_utils. This file contains no likelihood implementation.
"""

import json
from pathlib import Path

import numpy as np


TOOLKIT_COMMIT = "fc09848fc6540fd32310faebbe9db6eea7ecd17b"
PROCESSES = ("SBI", "S", "B", "qqZZ")


def write_workspace(directory, yields, metadata):
    """Write the toolkit's JSON workspace using the prepared .npy inputs.

    The toolkit's model accepts this dictionary directly. Its YAML/ROOT builder
    is unnecessary because our yields and density ratios already come from HDF5.
    Paths are absolute so the workspace can be loaded from any working directory
    in the same Drive mount. Initial factors define mu=1 Asimov data.
    """
    directory = Path(directory).resolve()
    samples = []
    for name in PROCESSES:
        samples.append({
            "name": name,
            "data": [float(yields[name])],
            "ratios": str(directory / f"r_{name}.npy"),
            "modifiers": [] if name == "qqZZ" else [{
                "name": f"scale_{name}", "type": "normfactor", "data": None,
            }],
        })

    workspace = {
        "version": "1.0.0",
        "channels": [{
            "name": "off_shell", "type": "unbinned",
            "weights": str(directory / "asimov_weights.npy"),
            "samples": samples,
        }],
        "measurements": [{
            "name": "stat_only",
            "config": {
                "poi": "scale_SBI",
                "parameters": [
                    {"name": "scale_SBI", "inits": [1.0]},
                    {"name": "scale_S", "inits": [0.0]},
                    {"name": "scale_B", "inits": [0.0]},
                ],
            },
        }],
        "metadata": {
            **metadata,
            "toolkit_commit": TOOLKIT_COMMIT,
            "physical_poi": "mu",
            "asimov_mu": 1.0,
            "scaling": {"SBI": "sqrt(mu)", "S": "mu - sqrt(mu)",
                        "B": "1 - sqrt(mu)", "qqZZ": "1"},
            "parameter_adapter": "utils_likelihood.mu_likelihood",
        },
    }
    path = directory / "workspace.json"
    path.write_text(json.dumps(workspace, indent=2, allow_nan=False) + "\n")
    return path


def mu_likelihood(model, mu_bounds=None):
    """Tie the toolkit's three normfactors to ONE physical parameter, mu.

    The pinned toolkit only supports multiplicative normfactors. These are
    coordinates of the physical curve, not three independently fitted parameters.
    Return the toolkit's own -2 log L unchanged; inference subtracts its minimum.
    Numerical differentiation by the toolkit's fitter avoids a sqrt derivative at
    mu=0. Reject points outside the checked interval. The toolkit fitter does
    not pass workspace bounds to Minuit, so the adapter enforces this domain
    during the fit as well as during the scan.
    """
    names, _ = model.get_model_parameters()
    expected = {"scale_SBI", "scale_S", "scale_B"}
    if set(names) != expected or len(names) != 3:
        raise ValueError(f"Unexpected toolkit parameters: {names}")
    if mu_bounds is None:
        mu_bounds = model.workspace.get("metadata", {}).get("scan_range", (0.0, np.inf))
    lower, upper = mu_bounds
    if not 0 <= lower < upper:
        raise ValueError(f"Invalid mu interval: {mu_bounds}")

    def nll(parameters):
        mu = float(parameters[0])
        if not np.isfinite(mu) or not lower <= mu <= upper:
            return np.inf
        root_mu = np.sqrt(mu)
        factors = {"scale_SBI": root_mu, "scale_S": mu - root_mu,
                   "scale_B": 1.0 - root_mu}
        return float(model.model(np.array([factors[name] for name in names])))

    return nll


def _quadratic_minimum(a, b, c, lower, upper):
    """Minimum of a*t**2 + b*t + c on a closed interval, element by element."""
    a, b, c = np.broadcast_arrays(a, b, c)
    left = (a * lower + b) * lower + c
    right = (a * upper + b) * upper + c
    minimum = np.minimum(left, right)
    location = np.where(left <= right, lower, upper)
    vertex = np.divide(-b, 2 * a, out=np.zeros_like(a), where=a > 0)
    vertex_value = (a * vertex + b) * vertex + c
    use_vertex = ((a > 0) & (lower < vertex) & (vertex < upper)
                  & (vertex_value < minimum))
    return np.where(use_vertex, vertex_value, minimum), np.where(use_vertex, vertex, location)


def check_workspace_domain(workspace_path, mu_bounds):
    """Check positive total intensities over the ENTIRE requested interval.

    This is input validation, not a likelihood calculation. With t=sqrt(mu),
    the intensity is a quadratic, so endpoints plus any interior minimum
    suffice. Check in chunks to keep memory small. Negative MC weights and
    negative individual process coefficients remain allowed and unchanged.
    """
    workspace_path = Path(workspace_path)
    workspace = json.loads(workspace_path.read_text())
    samples = workspace["channels"][0]["samples"]
    yields = {sample["name"]: float(sample["data"][0]) for sample in samples}
    ratios = {sample["name"]: np.load(sample["ratios"], mmap_mode="r", allow_pickle=False)
              for sample in samples}
    n = len(ratios["SBI"])
    if n == 0 or any(r.shape != (n,) for r in ratios.values()):
        raise ValueError("Workspace ratios must be nonempty aligned 1D arrays.")
    if not 0 <= mu_bounds[0] < mu_bounds[1] or not np.isfinite(mu_bounds).all():
        raise ValueError(f"Invalid scan interval: {mu_bounds}")
    lower, upper = np.sqrt(mu_bounds)
    worst, worst_mu, worst_position, invalid = np.inf, None, None, 0
    for start in range(0, n, 262144):
        stop = min(start + 262144, n)
        terms = {name: yields[name] * ratios[name][start:stop] for name in PROCESSES}
        if any(not np.isfinite(term).all() for term in terms.values()):
            raise ValueError("Nonfinite ratio or yield in the workspace inputs.")
        a = terms["S"]
        b = terms["SBI"] - terms["S"] - terms["B"]
        c = terms["B"] + terms["qqZZ"]
        values, roots = _quadratic_minimum(a, b, c, lower, upper)
        invalid += int(np.count_nonzero(values <= 0))
        position = int(np.argmin(values))
        if values[position] < worst:
            worst = float(values[position])
            worst_mu = float(roots[position] ** 2)
            worst_position = start + position

    rate, _ = _quadratic_minimum(
        yields["S"], yields["SBI"] - yields["S"] - yields["B"],
        yields["B"] + yields["qqZZ"], lower, upper,
    )
    report = {"mu_low": float(mu_bounds[0]), "mu_high": float(mu_bounds[1]),
              "asimov_events": n, "min_intensity": worst, "mu_at_min_intensity": worst_mu,
              "worst_event_position": worst_position, "nonpositive_events": invalid,
              "min_total_yield": float(rate)}
    rows_path = workspace_path.parent / "event_rows.npy"
    if rows_path.exists():
        rows = np.load(rows_path, mmap_mode="r", allow_pickle=False)
        report["worst_input_row"] = int(rows[worst_position])
    (workspace_path.parent / "domain_check.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    if invalid or rate <= 0:
        raise ValueError(
            f"{workspace_path.parent.name}: nonpositive combined intensity on {mu_bounds}; "
            f"{invalid} affected events, minimum {worst:.6g} at mu={worst_mu:.6g}. "
            "See domain_check.json. No events or densities were clipped."
        )
    return report


def plot_likelihood_scans(scan_arrays, labels, mu_bounds):
    """Draw toolkit curves with ATLAS formatting and a useful likelihood scale."""
    import matplotlib.pyplot as plt
    from nsbi_common_utils import inference

    fig, ax = plt.subplots(figsize=(9, 6.5))
    styles = ["solid", "dotted", "dashed", "dashdot", (0, (5, 2))]
    colors = ["black", "#D55E00", "#0072B2", "#009E73", "#CC79A7"]
    inference.plot_NLL_scans(
        parameter_label=r"$\mu$", list_scan_points=[scan[0] for scan in scan_arrays],
        list_nll_values=[scan[1] for scan in scan_arrays], list_labels=labels,
        list_linestyles=styles[:len(labels)], list_colors=colors[:len(labels)], ax=ax,
    )
    # The toolkit's 2/3-sigma lines should not flatten a short, low-NLL scan.
    ax.set_ylim(0, max(1.2, 1.15 * max(float(np.max(scan[1])) for scan in scan_arrays)))
    for text in ax.texts:
        text.set_clip_on(True)
    ax.set(xlim=mu_bounds, ylabel=r"$t_\mu=-2\Delta\log L$")
    ax.text(0.02, 1.06, "ATLAS", transform=ax.transAxes,
            fontsize=17, fontweight="bold", fontstyle="italic")
    ax.text(0.19, 1.06, "Simulation Internal", transform=ax.transAxes, fontsize=15)
    ax.text(0.02, 1.01, r"Stat-only · Asimov $\mu=1$ · Same-sample normalization",
            transform=ax.transAxes, fontsize=12)
    ax.legend(loc="best", fontsize=11, frameon=False)
    fig.tight_layout()
    return fig, ax
