"""Small I/O and plotting helpers for the off-shell HZZ notebooks."""

from io import BytesIO

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import tables


def load_columns(path, columns, key="dataset"):
    """Read selected columns from this ntuple's pandas fixed-format HDF5 group.

    Numeric columns are sliced directly, avoiding a full 66-column DataFrame.
    PyTables must deserialize the complete object block containing Type/subType;
    only the requested columns are retained. Use this with the trusted ntuple.
    """
    result = {}
    with tables.open_file(path, mode="r") as handle:
        group = handle.get_node("/" + key.strip("/"))
        for block_id in range(group._v_attrs.nblocks):
            names = handle.get_node(group, f"block{block_id}_items").read().astype(str).tolist()
            wanted = [name for name in columns if name in names]
            if not wanted:
                continue
            node = handle.get_node(group, f"block{block_id}_values")
            if isinstance(node, tables.VLArray):
                block = node[0]
                for name in wanted:
                    result[name] = block[:, names.index(name)].copy()
                del block
            else:
                for name in wanted:
                    result[name] = node[:, names.index(name)]
    missing = set(columns) - result.keys()
    if missing:
        raise KeyError(f"Missing columns: {sorted(missing)}")
    return {name: result[name] for name in columns}


def plot_comparison(reference_hist, target_hist, bins, *, xlabel, ylabel):
    """Draw two already-binned histograms; selection and weighting stay in the notebook."""
    fig, ax = plt.subplots(figsize=(8, 6))
    hep.histplot(reference_hist, bins=bins, ax=ax, histtype="step", yerr=False,
                 color="#0072B2", linewidth=2, label="Reference: S + EWSBI10")
    hep.histplot(target_hist, bins=bins, ax=ax, histtype="step", yerr=False,
                 color="#D55E00", linewidth=2, linestyle="--", label="Target: SBI + qqZZ")
    ax.set(xlabel=xlabel, ylabel=ylabel, xlim=(bins[0], bins[-1]))
    low = min(0.0, float(np.min(reference_hist)), float(np.min(target_hist)))
    high = max(0.0, float(np.max(reference_hist)), float(np.max(target_hist)))
    span = high - low if high > low else 1.0
    ax.set_ylim(low - 0.06 * span if low < 0 else 0, high + 0.38 * span)
    if low < 0:
        ax.axhline(0, color="0.5", linewidth=0.8)
    # mplhep's label helper uses a private text-layout API that differs across
    # Colab's Matplotlib versions. Use public Matplotlib text for this label.
    ax.text(0.04, 0.95, "ATLAS", transform=ax.transAxes, va="top",
            fontsize=17, fontweight="bold", fontstyle="italic")
    ax.text(0.21, 0.95, "Simulation Internal", transform=ax.transAxes,
            va="top", fontsize=15)
    ax.legend(loc="upper right", fontsize=13, frameon=False)
    fig.tight_layout()
    return fig, ax


def check_plot_exports():
    """Check the plotting stack and both export backends before reading the ntuple."""
    fig, _ = plot_comparison([0.8, 1.3], [1.7, 0.4], [0.0, 1.0, 2.0],
                             xlabel="Export check", ylabel="Export check")
    try:
        for extension in ("png", "pdf"):
            with BytesIO() as buffer:
                fig.savefig(buffer, format=extension, bbox_inches="tight")
    finally:
        plt.close(fig)


def normalize_histogram(hist, variance, bins):
    """Return a unit-area density and its diagonal delta-method variance.

    The normalization is the signed in-range sum. The input variance is sum(w^2)
    in each bin; its covariance with the normalization is included.
    """
    hist = np.asarray(hist, dtype=float)
    variance = np.asarray(variance, dtype=float)
    total = hist.sum()
    if total == 0:
        raise ValueError("Cannot normalize a histogram with zero signed yield.")
    fraction = hist / total
    other_variance = np.maximum(variance.sum() - variance, 0.0)
    density_variance = ((1.0 - fraction)**2 * variance
                        + fraction**2 * other_variance) / (total * np.diff(bins))**2
    return fraction / np.diff(bins), density_variance


def plot_reweighting_comparison(reweighted_hist, target_hist, reweighted_variance,
                               target_variance, bins, *, xlabel, ylabel,
                               title, ratio_ylim=None):
    """Overlay histograms and show reweighted/direct with separate MC errors.

    Points show reweighted-reference uncertainty with the denominator fixed.
    The grey band shows the direct-target uncertainty around unity. Bins with
    zero direct yield have an undefined ratio and are omitted, never set to 1.
    """
    reweighted_hist = np.asarray(reweighted_hist, dtype=float)
    target_hist = np.asarray(target_hist, dtype=float)
    reweighted_error = np.sqrt(np.asarray(reweighted_variance, dtype=float))
    target_error = np.sqrt(np.asarray(target_variance, dtype=float))
    bins = np.asarray(bins, dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])
    fig, (ax, ratio_ax) = plt.subplots(
        2, 1, figsize=(8, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06},
    )

    # With bare arrays, yerr=None draws steps without inferring Poisson errors.
    # The signed-weight MC errors are drawn explicitly below.
    hep.histplot(target_hist, bins=bins, ax=ax, histtype="step", yerr=None,
                 color="black", linewidth=1.7, label="Direct target")
    hep.histplot(reweighted_hist, bins=bins, ax=ax, histtype="step", yerr=None,
                 color="#0072B2", linewidth=2, linestyle="--",
                 label="Reweighted reference")
    lower, upper = target_hist - target_error, target_hist + target_error
    ax.fill_between(bins, np.r_[lower, lower[-1]], np.r_[upper, upper[-1]],
                    step="post", color="0.5", alpha=0.25, linewidth=0,
                    label="Target MC stat.")
    ax.errorbar(centers, reweighted_hist, yerr=reweighted_error, fmt="none",
                color="#0072B2", elinewidth=1, alpha=0.8)
    ax.set(ylabel=ylabel, xlim=(bins[0], bins[-1]))
    ax.set_title(title, fontsize=16, pad=12)
    low = min(0.0, float(lower.min()), float((reweighted_hist - reweighted_error).min()))
    high = max(0.0, float(upper.max()), float((reweighted_hist + reweighted_error).max()))
    span = high - low if high > low else 1.0
    ax.set_ylim(low - 0.06 * span if low < 0 else 0, high + 0.40 * span)
    if low < 0:
        ax.axhline(0, color="0.5", linewidth=0.8)
    ax.text(0.04, 0.96, "ATLAS", transform=ax.transAxes, va="top",
            fontsize=17, fontweight="bold", fontstyle="italic")
    ax.text(0.23, 0.96, "Simulation Internal", transform=ax.transAxes,
            va="top", fontsize=13)
    ax.legend(loc="upper right", fontsize=11, frameon=False)

    valid = target_hist != 0
    ratio = np.full_like(target_hist, np.nan)
    ratio_error = np.full_like(target_hist, np.nan)
    target_relative_error = np.full_like(target_hist, np.nan)
    ratio[valid] = reweighted_hist[valid] / target_hist[valid]
    ratio_error[valid] = reweighted_error[valid] / np.abs(target_hist[valid])
    target_relative_error[valid] = target_error[valid] / np.abs(target_hist[valid])
    lower, upper = 1.0 - target_relative_error, 1.0 + target_relative_error
    ratio_ax.fill_between(bins, np.r_[lower, lower[-1]], np.r_[upper, upper[-1]],
                          step="post", color="0.5", alpha=0.25, linewidth=0)
    ratio_ax.axhline(1, color="black", linestyle="--", linewidth=1)
    ratio_ax.errorbar(centers[valid], ratio[valid], yerr=ratio_error[valid],
                      fmt="o", markersize=3, color="#0072B2", elinewidth=1)
    ratio_ax.set(xlabel=xlabel, ylabel="Reweighted\n/ direct")
    ratio_ax.yaxis.label.set_size(15)
    ratio_ax.tick_params(axis="both", labelsize=13)
    ratio_ax.grid(axis="y", alpha=0.2)
    if ratio_ylim is not None:
        ratio_ax.set_ylim(*ratio_ylim)
    if not np.all(valid):
        ratio_ax.text(0.02, 0.94, f"Zero-target bins omitted: {np.count_nonzero(~valid)}",
                       transform=ratio_ax.transAxes, va="top", fontsize=9)
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.12, top=0.93)
    return fig, (ax, ratio_ax)
