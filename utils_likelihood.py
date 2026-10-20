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


def mu_likelihood(model):
    """Tie the toolkit's three normfactors to ONE physical parameter, mu.

    The pinned toolkit only supports multiplicative normfactors. These are
    coordinates of the physical curve, not three independently fitted parameters.
    Return the toolkit's own -2 log L unchanged; inference subtracts its minimum.
    Numerical differentiation by the toolkit's fitter avoids a sqrt derivative at
    mu=0. Reject mu<0 rather than clipping or taking its absolute value.
    """
    names, _ = model.get_model_parameters()
    expected = {"scale_SBI", "scale_S", "scale_B"}
    if set(names) != expected or len(names) != 3:
        raise ValueError(f"Unexpected toolkit parameters: {names}")

    def nll(parameters):
        mu = float(parameters[0])
        if not np.isfinite(mu) or mu < 0:
            return np.inf
        root_mu = np.sqrt(mu)
        factors = {"scale_SBI": root_mu, "scale_S": mu - root_mu,
                   "scale_B": 1.0 - root_mu}
        return float(model.model(np.array([factors[name] for name in names])))

    return nll
