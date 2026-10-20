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
