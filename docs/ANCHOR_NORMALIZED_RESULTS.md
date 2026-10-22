# Asimov workspaces with normalized systematic anchors

The corrected 4096 and 8192 workspaces return unconditional signal strengths
**1.001444013** and **1.001201170**, respectively. All 127 systematic nuisance
parameters and all three cumulative qqZZ factors float. Both original global
jobs finished uninterrupted before the new profile jobs were submitted.

## Diagnosis and correction

The three control regions were already Asimov: their counts are
1244.3595984565, 250.2201341248, and 178.5237118586, agreeing with their
nominal model expectations within 2.3e-13. Fixing all systematic parameters
to zero returns mu approximately 0.99996 and 0.99998. The displacement in
the original fully profiled fits is associated with a nonzero systematic
score at the generating point, where

`d NLL_A / d theta = d nu / d theta - sum_i u_i d I_i / d theta`.

For example, the old 8192 HOQCD_1Jet score is 1.43052737. The independent
expected-yield derivative minus integrated-intensity derivative gives
4.62077486 - 3.19024750 = 1.43052736. The Asimov control regions and zero
Gaussian pulls contribute zero score at generation.

Each supplied event systematic variation contains both shape and rate.
For each subset, SR component, nuisance and up/down direction, the correction is

`C_pk = sum_i u_i rhat_pi v_pki`,
`v_corrected_pki = v_pki R_pk / C_pk`.

The target `R_pk` is the physical full-sample rate factor from
`list_weight_up/down_<process>.npy`; it is not forced to one. The three
conditional qqZZ jet components are treated independently. Each workspace
passes 2540 anchor-integral checks and 258 native-versus-explicit likelihood
settings. The event selections, nominal ratios, signed weights, yields,
Asimov control regions, constraints, interpolation and bounds are unchanged.

## Completed unconditional fits

| Reference events | Selected mu hat | Percent above one | Maximum absolute alpha | EDM | NLL spread between starts |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4096 | 1.001444013 | 0.144401% | 0.0497735 | 2.54e-8 | 5.48e-8 |
| 8192 | 1.001201170 | 0.120117% | 0.0486170 | 2.46e-8 | 3.16e-8 |

Starts at mu=1 and 1.4 all converge with status 0. The saved MIGRAD covariance
estimates have quality 3 and positive eigenvalues; no explicit full Hesse
calculation was performed. Parameter CSV errors are MIGRAD estimates.
The selected fitted-point SR integral mismatch is 3.20e-6 and 2.93e-6
relative, versus +18.60% and +0.92% at the original first converged starts.
The old 4096 mu=0.12618 solution is a known higher local minimum: a later
profile lowered NLL by 86.162, so that old baseline and curve are invalid.
The old 8192 unconditional estimate was 1.104358886. All original outputs
remain preserved separately.

## Remaining normalization limitation

The largest absolute systematic score at generation decreases from 1.87048
to 0.17141 for 4096, and from 1.43053 to 0.15515 for 8192. Correcting anchors
does not normalize nonlinear interpolation or simultaneous nuisance
variations at every point. Native PiecewiseInterpolation code 5 and rate
FlexibleInterpVar code 4 are retained. This explains why exact stationarity
is not guaranteed by this experiment.

The distinction between anchor closure and complete nuisance-dependent
normalization is discussed in [arXiv:2609.30196, Sections 4.2.1 and 5.5](https://arxiv.org/html/2609.30196#S5.SS5).
The paper's shape-only example motivates the diagnosis; here physical rate
variations must be preserved. The 8192 sample retains two negative event
weights, so positive-weight global-maximum proofs require separate care.

## Profile convention

Evaluate mu=0.01, then 0.05, 0.10, ..., 1.50. Profile all 130 nuisance
parameters and normalize q to the saved continuous unconditional minimum.
Never evaluate the likelihood at exactly zero. Compare with
[HEPData 156805 v1, Table 2, Exp NSBI](https://doi.org/10.17182/hepdata.156805.v1/t2)
at the 30 shared coordinates, without shifting, rescaling or interpolating
the reference. The published zero point is displayed at zero, separately
from our 0.01 point; the published sampled minimum is left unchanged.

## Completed profiles and HEPData residuals

Both samples converged at all 31 profile coordinates. No profile point
improves the saved unconditional minimum. At the 30 coordinates shared
with HEPData, the unshifted differences are:

| Reference events | Maximum absolute Delta q | RMS Delta q | Mu at largest residual | Local q(0.01) |
| --- | ---: | ---: | ---: | ---: |
| 4096 | 0.09051853 | 0.04375547 | 0.20 | 1.99505423 |
| 8192 | 0.16462302 | 0.05993946 | 0.05 | 1.89321439 |

For these particular nested selections, the corrected 4096 curve is closer
to HEPData. This does not establish that smaller reference samples are
generally more accurate. Between the two local curves, the largest absolute
difference is 0.10183984 at mu=0.01, and the RMS difference is 0.02758312 over
all 31 coordinates. The local endpoint has no exact published counterpart
and is excluded from the HEPData residual statistics.

Finite sampling, incomplete nuisance-dependent normalization, and differences
between the supplied NumPy physics inputs and the published model can all
affect residuals. These two nested samples do not isolate those effects.
The scans are likelihood-ratio curves, not the publication's Neyman belts.

## Independent profile validation

Both merged profiles pass all 31 independent native-versus-explicit
likelihood evaluations and fitted-point positivity checks. The original
continuous global-fit baseline is retained in both cases.

| Reference events | Maximum profile EDM | Maximum overlap Delta q | Maximum Delta-NLL disagreement | Minimum fitted intensity | Maximum relative SR integral mismatch |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4096 | 6.84e-6 | 2.42e-5 | 1.53e-10 | 0.39204 | 4.55e-6 |
| 8192 | 7.14e-6 | 2.23e-5 | 2.62e-10 | 0.22817 | 4.32e-6 |

The three overlapping chunk coordinates are mu=0.4, 0.8 and 1.2; every
overlap agrees within the required 1e-4 in q. Workspace and numerical-source
checksums match the launch manifests, and both old workspaces remain unchanged.

The global workflow used source commit
[`07edd50`](https://github.com/rafaellopesdesa/offshell_reweighting/commit/07edd50330768454416f8fc42cc2815040c23f8e),
and the profile/validation/comparison jobs used
[`e2a3ed3`](https://github.com/rafaellopesdesa/offshell_reweighting/commit/e2a3ed3de56c54faddb95d480eba3998f828d6fd).

## Measured execution time and memory

All 13 jobs completed successfully. Slurm accounting gives:

| Reference events | Stage | Elapsed | Total CPU minutes | Largest individual job peak RSS |
| --- | --- | ---: | ---: | ---: |
| 4096 | Construction + two global fits + checks | 24m11s | 23.95 | 1.13 GiB |
| 8192 | Construction + two global fits + checks | 47m49s | 47.28 | 1.98 GiB |
| 4096 | Four parallel profile chunks | 49m03s | 145.90 | 1.19 GiB |
| 8192 | Four parallel profile chunks | 96m44s | 287.15 | 1.05 GiB |
| 4096 | Merge, figures + independent profile validation | 4m26s | 4.37 | 0.93 GiB |
| 8192 | Merge, figures + independent profile validation | 7m37s | 7.53 | 1.35 GiB |

Each job requested one CPU and 4 GiB. Elapsed time spans the first start to
the last finish in a stage; TotalCPU sums the jobs without double-counting
their steps. Peak RSS is the largest measured individual job peak, not the
simultaneous aggregate memory. Different CPU nodes were used, so this is
not a controlled hardware benchmark.

The global stage includes construction from already staged inputs, both
unconditional starts, and their checks, without a full Hesse. ROOT
installation, earlier input staging, and waiting before each stage's first
start are excluded. The shared final HEPData comparison took 7 seconds
elapsed and is excluded from the per-sample rows. Field definitions follow
[Slurm sacct](https://slurm.schedmd.com/sacct.html).

## Downloadable results

The [aggregate report directory](../reports/anchor_normalized/) contains
the full-precision summary, both scan tables, all selected global-fit
parameters, independent per-point validation, and measured resource tables.
Its source hashes identify the saved inputs; raw events, ROOT workspaces,
and full Slurm logs remain in the work area. Parameter errors are MIGRAD
estimates. The [presentation branch](https://github.com/rafaellopesdesa/offshell_reweighting/tree/codex/asimov-results-slides-20261021)
provides the comparison figures and discussion as a 17-slide Reveal.js deck.

## Separate exact-zero 4096 test

At the user's subsequent request, the corrected 4096 workspace was profiled
at **exactly mu=0**, with all 130 nuisances free and the original global-fit
baseline unchanged. The usual fitting procedure converged:

- **q(0) = 2.8262024415**, from profile NLL -19979.032433097178 and
  global NLL -19980.44553431794.
- MIGRAD and fit status 0, EDM 3.48e-6, covariance quality 3, and zero
  invalid NLL evaluations. No explicit full Hesse was performed.
- Independent native/explicit Delta-NLL agreement 1.02e-10; all fitted
  intensities positive, with minimum 19.16396. Relative SR integral
  mismatch is 4.51e-6.
- HEPData's expected q(0) is 2.954, giving local minus published
  **-0.1277975585**. The earlier local q(0.01)=1.9950542283 is a different
  coordinate; the rise to exact zero is 0.8311482132.

The profile started from the saved global nuisance vector. The largest
conditional absolute alpha is 0.17298. Numerical source
[`177916f`](https://github.com/rafaellopesdesa/offshell_reweighting/commit/177916fbbbf58fb680c95687f3b129179b836284)
adds only an explicit zero opt-in and a separate batch workflow; it changes
neither the native likelihood nor the fit bounds. The original workspace,
global artifacts and completed 31-point scan pass their unchanged-file checks.

The Slurm job took **8m20s elapsed**, **481.48 CPU seconds**, and **1.35 GiB
peak RSS**, including loading and independent validation. It requested one
CPU and 2 GiB. The minimization itself took 438.19 seconds.

Download the [zero-point summary](../reports/anchor_normalized/zero_4096.json)
and [conditional nuisance values](../reports/anchor_normalized/zero_4096_parameters.csv).
The point remains a separate exploratory result; the 31-point tables,
30-coordinate comparison statistics and slideshow snapshot above are unchanged.
