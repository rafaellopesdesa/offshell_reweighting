# Consistent 4096 and 8192 workspace family

All six current workspaces use the unchanged version-1
`RooSameSampleNormalizer`. Each has ten instances, one for each SR process
component, including the three conditional qqZZ jet components. The class
normalizes the fully interpolated shape on the fixed reference sample at
every nuisance point. The full-MC physical rate variations remain separate.

The original true 4096 and 8192 models normalized only nominal/up/down anchors.
They and their fitted results remain in `results/asimov_anchor_normalized/`.
The replacement true files have `true_asimov` in their filenames and are
separate artifacts. Their old profile curves must not be assigned to them.
The existing unpruned and pruned data-like 4096 files already used the class;
their bytes are unchanged.

## Generation and pruning

Every model generates at mu=1. True Asimov generation has all alpha=0 and all
three qqZZ factors=1. Data-like generation injects the observed-data global
fit's alpha and cumulative qqZZ values, with CR counts and Gaussian auxiliary
observations generated consistently at that point. The model's fixed
reference measure is distinct from the generated Asimov observation weights.

The 8192 reference selection is the original nested extension of the 4096
selection. Normalization and qqZZ fractions are specific to each size. No
events are redrawn and no physical yields are rescaled with sample size.
Its two negative reference weights, and resulting negative observation
weights, are preserved. Generating-point stationarity can be checked for
this signed quadrature; it does not inherit the nonnegative-weight guarantee
of a global likelihood minimum.

Both pruned versions retain the same first 54 constrained NPs from
`reports/data_like_4096/decomposition/all_np_uncertainty_ranking.csv`, together
with all three free qqZZ factors. The remaining 73 alphas are fixed at the
data-like generating values. Their fixed shape/rate factors remain in the
model. The 8192 selection uses this existing ranking; no new uncertainty
ranking has been estimated for 8192.

| Variant, at either size | Constrained NPs | Free qqZZ factors | Global floating parameters |
| --- | ---: | ---: | ---: |
| True | 127 | 3 | 131 |
| Data-like | 127 | 3 | 131 |
| Data-like, pruned | 54 | 3 | 58 |

## Files and use

The work-area archive is `results/asimov_same_sample_normalized/`, relative
to the Codex work area. `workspaces/` contains the six ROOT files ready to
copy to Drive; `workspace_family.json` and `SHA256SUMS` identify every file.
The published catalog and validation summary are under
`reports/same_sample_workspace_family/` in this analysis branch.

Load the class dictionary before retrieving the ROOT workspace, as described
in [the normalizer documentation](ROO_SAME_SAMPLE_NORMALIZER.md). The reference
measure and interpolation anchors are embedded; no NumPy sidecars are needed
for fitting. The analysis construction/validation scripts use sidecars for
additional checks and provenance.

The demo branch selects the updated three 4096 files in `workspaces.json`.
It records the 8192 archive in a separate `workspaces_8192.json`; those files
are not included in notebook 06's default scan. Its results directory has a
new `full_normalization_v1` component to preserve old checkpoints/results.

## Reproduction and validation

Run inside a Slurm allocation with ROOT >=6.34.04. For each size, use the
corresponding preserved anchor workspace as `ANCHOR` and a new output directory:

```bash
python build_shifted_asimov.py --source "$ANCHOR" --generation nominal --output "$TRUE_OUTPUT"
python build_shifted_asimov.py --source "$ANCHOR" --generation data-like \
  --observed-bestfit "$OBSERVED_BESTFIT" --output "$DATA_LIKE_OUTPUT"
python prune_data_like_asimov.py --input "$DATA_LIKE_WORKSPACE" --keep 54 \
  --ranking reports/data_like_4096/decomposition/all_np_uncertainty_ranking.csv \
  --output "$PRUNED_OUTPUT"
```

The builders refuse to overwrite existing ROOT files. Construction validation
checks serialization, exact reference preservation, consistent CR/auxiliary
observations, closure between and beyond individual anchors, and numerical
derivatives of all floating parameters at generation. The normalizer's cached
value is checked against its independent direct interpolation sum.
Independent intensity calculations test RooFit likelihood differences and
positivity over [0,2] at selected fixed nuisance settings. Pruning checks exact
agreement with the corresponding restriction of the parent model.

`audit_asimov_normalizers.py` checks actual class instances and graph connections.
`validate_asimov_model_pair.py` verifies that each true/data-like pair now has
identical shapes, physical rates and CR expectations at common parameter
points, and that the corrected true model preserves central observations.
No new minimizations, Hesse calculations, or likelihood scans are part of
this archival update.
