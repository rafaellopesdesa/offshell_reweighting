# Plotting preference

Use a ratio-axis range of **0.90 to 1.10** for all histogram ratio panels, including future notebooks, unless Rafael explicitly requests a different range. The shared default is `DEFAULT_RATIO_YLIM` in `utils.py`.

Keep sample selections, weights, and NumPy histogram calculations readable in the notebooks. Put lengthy plot formatting in `utils.py`.

# Likelihood normalization

Use same-sample normalization for every process in likelihood studies. Compute
Z with the normalized signed weights of the exact reference event set used in
that workspace, and recompute Z for each random subset. Keep full native process
yields fixed when reducing the reference event count, unless Rafael explicitly
asks to change the physical yield.
