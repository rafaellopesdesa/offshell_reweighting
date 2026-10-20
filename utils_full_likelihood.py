"""Workspace and parameter adapters for notebook 05; no likelihood implementation.

The NSBI toolkit evaluates the extended likelihood and performs all fits/scans.
Only workspace formatting, physical coefficients, and input checks live here.
"""

import json
import warnings
from pathlib import Path

import numpy as np

from utils_likelihood import TOOLKIT_COMMIT, _quadratic_minimum


PROCESSES = ("SBI", "S", "B", "EWB", "EWSBI", "EWSBI10", "qqZZ", "ttV")
SAMPLES = tuple(p for p in PROCESSES if p != "qqZZ") + ("qqZZ0", "qqZZ1", "qqZZ2")
PARAMETERS = ("mu", "mu_qqZZ0", "mu_qqZZ1", "mu_qqZZ2")


def coefficients(mu):
    """Physical basis coefficients, including negative interference terms."""
    t = np.sqrt(mu)
    root10 = np.sqrt(10.0)
    d = -10.0 + root10
    return {
        "SBI": t, "S": mu - t, "B": 1.0 - t,
        "EWB": ((1.0 - root10) * mu + 9.0 * t - 10.0 + root10) / d,
        "EWSBI": (root10 * mu - 10.0 * t) / d,
        "EWSBI10": (-mu + t) / d,
        "ttV": 1.0,
    }


def write_workspace(directory, yields, jet_fractions, metadata):
    """Encode each masked qqZZ component with its integral and unit-mean ratio.

    data = lambda_qqZZ * f_c and ratio = (r_qqZZ/Z_qqZZ) * I_c / f_c.
    Their product is exactly lambda_qqZZ * (r_qqZZ/Z_qqZZ) * I_c.
    The three yields sum to the inclusive qqZZ yield, not three times it.
    """
    directory = Path(directory).resolve()
    fractions = np.asarray(jet_fractions, dtype=float)
    if fractions.shape != (3,) or not np.isfinite(fractions).all() or (fractions <= 0).any():
        raise ValueError("All three reference qqZZ jet fractions must be finite and positive.")
    np.testing.assert_allclose(fractions.sum(), 1.0, atol=1e-10, rtol=0)
    nominal = coefficients(1.0)
    samples, parameters = [], []
    for name in SAMPLES:
        is_qq = name.startswith("qqZZ")
        value = yields["qqZZ"] * fractions[int(name[-1])] if is_qq else yields[name]
        factor = f"mu_{name}" if is_qq else f"scale_{name}"
        modifiers = [] if name == "ttV" else [{"name": factor, "type": "normfactor", "data": None}]
        samples.append({"name": name, "data": [float(value)],
                        "ratios": str(directory / f"r_{name}.npy"), "modifiers": modifiers})
        if modifiers:
            parameters.append({"name": factor, "inits": [1.0 if is_qq else float(nominal[name])]})
    workspace = {
        "version": "1.0.0",
        "channels": [{"name": "off_shell", "type": "unbinned",
                      "weights": str(directory / "asimov_weights.npy"), "samples": samples}],
        "measurements": [{"name": "full_stat_only", "config": {
            "poi": "scale_SBI", "parameters": parameters,
        }}],
        "metadata": {
            **metadata, "toolkit_commit": TOOLKIT_COMMIT,
            "physical_parameters": list(PARAMETERS), "asimov_parameters": [1.0] * 4,
            "parameter_adapter": "utils_full_likelihood.physical_likelihood",
            "qqZZ_jet_fractions": fractions.tolist(),
            "qqZZ_inclusive_yield": float(yields["qqZZ"]),
            "qqZZ_categories": ["n_jets == 0", "n_jets == 1", "n_jets >= 2"],
            "scaling": {
                "SBI": "sqrt(mu)", "S": "mu - sqrt(mu)", "B": "1 - sqrt(mu)",
                "EWB": "((1-sqrt(10))*mu + 9*sqrt(mu) - 10 + sqrt(10))/(-10+sqrt(10))",
                "EWSBI": "(sqrt(10)*mu - 10*sqrt(mu))/(-10+sqrt(10))",
                "EWSBI10": "(-mu + sqrt(mu))/(-10+sqrt(10))", "ttV": "1",
                "qqZZ0": "mu_qqZZ0", "qqZZ1": "mu_qqZZ1", "qqZZ2": "mu_qqZZ2",
            },
            "nuisance_constraints": "none; qqZZ normalizations are nonnegative with no upper bound",
        },
    }
    path = directory / "workspace.json"
    path.write_text(json.dumps(workspace, indent=2, allow_nan=False) + "\n")
    return path


def physical_likelihood(model, mu_bounds=None):
    """Expose mu and three independent, unpenalized qqZZ normalizations.

    The six internal signal/interference factors stay tied to mu. The toolkit
    computes -2 log L, including the rate and excluding any constraint penalty.
    Nonnegative qqZZ factors are a physical domain, not a Gaussian constraint.
    Invalid intensities are rejected, never clipped. Numerical differentiation
    by the toolkit avoids introducing a singular analytic derivative at mu=0.
    """
    names, _ = model.get_model_parameters()
    expected = {f"scale_{p}" for p in PROCESSES if p not in ("qqZZ", "ttV")}
    expected.update(PARAMETERS[1:])
    if len(names) != len(expected) or set(names) != expected:
        raise ValueError(f"Unexpected toolkit parameters: {names}")
    if model.num_unconstrained_param != len(names):
        raise ValueError("Notebook 05 must have no constrained workspace parameters.")
    if mu_bounds is None:
        mu_bounds = model.workspace["metadata"]["scan_range"]
    lower, upper = mu_bounds
    if not 0 <= lower < upper or not np.isfinite(mu_bounds).all():
        raise ValueError(f"Invalid mu interval: {mu_bounds}")

    def nll(parameters):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (4,):
            raise ValueError("Expected [mu, mu_qqZZ0, mu_qqZZ1, mu_qqZZ2].")
        mu = parameters[0]
        if (not np.isfinite(parameters).all() or not lower <= mu <= upper
                or (parameters[1:] < 0).any()):
            return np.inf
        factors = {f"scale_{p}": value for p, value in coefficients(mu).items() if p != "ttV"}
        factors.update(zip(PARAMETERS[1:], parameters[1:]))
        value = float(model.model(np.array([factors[name] for name in names])))
        return value if np.isfinite(value) else np.inf

    return nll


def _intensity_polynomial(terms):
    """Quadratic in sqrt(mu), at nominal qqZZ normalizations, for validation."""
    root10 = np.sqrt(10.0)
    d = -10.0 + root10
    a = terms["S"] + ((1-root10)*terms["EWB"] + root10*terms["EWSBI"] - terms["EWSBI10"]) / d
    b = terms["SBI"] - terms["S"] - terms["B"] + (9*terms["EWB"] - 10*terms["EWSBI"] + terms["EWSBI10"]) / d
    c = terms["B"] + terms["EWB"] + terms["ttV"] + sum(terms[f"qqZZ{k}"] for k in range(3))
    return a, b, c


def check_workspace_domain(workspace_path, mu_bounds):
    """Check the complete interval at nominal qqZZ factors, including vertices.

    This establishes a valid nominal starting curve. Profiling can visit other
    qqZZ factors; physical_likelihood rejects their nonfinite likelihoods too.
    """
    path = Path(workspace_path)
    workspace = json.loads(path.read_text())
    samples = workspace["channels"][0]["samples"]
    yields = {s["name"]: float(s["data"][0]) for s in samples}
    ratios = {s["name"]: np.load(s["ratios"], mmap_mode="r", allow_pickle=False) for s in samples}
    n = len(ratios["SBI"])
    if n == 0 or any(r.shape != (n,) for r in ratios.values()):
        raise ValueError("Workspace ratios must be nonempty aligned 1D arrays.")
    if not 0 <= mu_bounds[0] < mu_bounds[1] or not np.isfinite(mu_bounds).all():
        raise ValueError(f"Invalid scan interval: {mu_bounds}")
    lower, upper = np.sqrt(mu_bounds)
    worst, worst_mu, worst_position, invalid = np.inf, None, None, 0
    for start in range(0, n, 262144):
        stop = min(start + 262144, n)
        terms = {p: yields[p] * ratios[p][start:stop] for p in SAMPLES}
        if any(not np.isfinite(t).all() for t in terms.values()):
            raise ValueError("Nonfinite ratio or yield in the workspace inputs.")
        values, roots = _quadratic_minimum(*_intensity_polynomial(terms), lower, upper)
        invalid += int(np.count_nonzero(values <= 0))
        position = int(np.argmin(values))
        if values[position] < worst:
            worst, worst_mu, worst_position = float(values[position]), float(roots[position]**2), start+position
    rate, _ = _quadratic_minimum(*_intensity_polynomial(yields), lower, upper)
    report = {"mu_low": float(mu_bounds[0]), "mu_high": float(mu_bounds[1]),
              "qqZZ_factors": [1.0]*3, "asimov_events": n,
              "min_intensity": worst, "mu_at_min_intensity": worst_mu,
              "worst_event_position": worst_position, "nonpositive_events": invalid,
              "min_total_yield": float(rate)}
    rows = np.load(path.parent / "event_rows.npy", mmap_mode="r", allow_pickle=False)
    report["worst_input_row"] = int(rows[worst_position])
    (path.parent / "domain_check.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if invalid or rate <= 0:
        raise ValueError(
            f"{path.parent.name}: nominal full-model intensity is nonpositive on {mu_bounds}; "
            f"{invalid} events, minimum {worst:.6g} at mu={worst_mu:.6g}. "
            "Inspect domain_check.json before fitting. No densities were clipped."
        )
    return report


def profile_scan(fitter, mu_bounds, size):
    """Delegate profiling to the toolkit, treating failed minimizations as errors.

    The pinned toolkit discards iminuit's per-point success flags, but iminuit
    also warns when a profile minimization fails. Do not silently plot those.
    """
    from iminuit.util import IMinuitWarning

    with warnings.catch_warnings():
        warnings.simplefilter("error", IMinuitWarning)
        points, values = fitter.perform_profile_scan(
            parameter_name="mu", bound_range=mu_bounds, size=size,
        )
    points, values = np.asarray(points), np.asarray(values)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite profiled likelihood; inspect the model domain and fit output.")
    return points, values
