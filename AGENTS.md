# Plotting preference

Use a ratio-axis range of **0.90 to 1.10** for all histogram ratio panels, including future notebooks, unless Rafael explicitly requests a different range. The shared default is `DEFAULT_RATIO_YLIM` in `utils.py`.

Keep sample selections, weights, and NumPy histogram calculations readable in the notebooks. Put lengthy plot formatting in `utils.py`.

# Likelihood normalization

Use **mu in [0, 1.8]** as the default likelihood scan interval, including future
notebooks, unless Rafael explicitly requests another range. Keep the continuous
intensity-domain checks; the range alone does not guarantee positivity.

Use same-sample normalization for every process in likelihood studies. Compute
Z with the normalized signed weights of the exact reference event set used in
that workspace, and recompute Z for each random subset. Keep full native process
yields fixed when reducing the reference event count, unless Rafael explicitly
asks to change the physical yield.

Distinguish the Asimov events from the model's normalization reference. When
comparing a direct native-process Asimov with the full-reference reweighted Asimov,
keep the full-reference Z factors for both models and retain the original
direct event weights. Do not recalibrate Z on the direct Asimov sample.

For the full model, also retain the full-reference qqZZ jet fractions in the
direct comparison. Split qqZZ using the shared inclusive yield and ratio:
f_c = sum_ref(u * rhat_qqZZ * I_c), yield_c = lambda_qqZZ * f_c,
ratio_c = rhat_qqZZ * I_c / f_c. Recompute f_c for every reference subset;
the inclusive native qqZZ yield stays fixed. Profile all three free qqZZ
normalizations, with no nuisance constraint terms, even in stat-only scans.

# Ntuple process labels

Use **subType only** to identify physical processes in
`dataset_wDR_mean_origIndex.h5`: `subType == process_name` for every process.
Do not select or filter processes using Type, and do not request Type in analysis
notebooks. In particular, qqZZ is `subType == "qqZZ"` and ttV is
`subType == "ttV"`; they are disjoint without additional conditions.

Apply these exact masks consistently to the reference, native yields, direct
targets/Asimov events, and comparisons. The reference is subType in S and
EWSBI10. Notebooks 01–04 exclude ttV; 05 includes it separately. Keep the stored
branch names `r_qqZZ` and `r_ttV` and retain the direct-Asimov yield consistency
check. Saved workspace metadata must identify subType as the process label.
The auxiliary label census may still inspect both Type and subType.
