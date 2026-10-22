# RooSameSampleNormalizer

This is the unchanged version-1 class used to serialize all three current
4096 workspaces and their archived 8192 counterparts. At every nuisance point it integrates the interpolated process
shape over the fixed, same-sample reference measure. Physical rate modifiers
remain separate; the Asimov observation weights do not define this measure.
The historical anchor-only true-Asimov file is preserved outside this manifest.

With ROOT >= 6.34.04 (including RooFit and PyROOT) and a C++ compiler:

```bash
python normalizer/build.py --output build/normalizer
```

The build script uses ROOT's compiler flags and rootcling to build
`build/normalizer/libRooSameSampleNormalizer.so` and its I/O dictionary.
Keep the generated dictionary files alongside the library. Rebuild with
the ROOT installation used for the fit; do not transfer the Unity binary
to Colab. The notebook builds locally in `/content`, not Google Drive.

Load the library **before reading the workspace object**:

```python
import ROOT
assert ROOT.gSystem.Load("build/normalizer/libRooSameSampleNormalizer.so") >= 0
source = ROOT.TFile.Open("HZZ_4l_data_like_asimov_4096.root", "READ")
workspace = source.Get("combined")
```

The reference measure and interpolation anchors are embedded in the ROOT
file, so no `.npy` files are needed. ROOT dictionary generation is described
in the [ROOT manual](https://root.cern/manual/io_custom_classes/).
