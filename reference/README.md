# Published expected scan

`hepdata_expected_4l.csv` is the **Exp NSBI** (statistical + systematic)
series from [HEPData 156805 v1, Table 2](https://doi.org/10.17182/hepdata.156805.v1/t2),
four-lepton channel, corresponding to Figure 11(b). It is not the observed
series, the stat-only series, or the combination with another channel.

The 41 points in the demonstration interval, mu = 0, 0.05, ..., 2, were extracted from the cached official
JSON and checked exactly against the official CSV download. The q values
are copied without a baseline shift or rescaling (q(0)=2.954, q(2)=1.316).
Only identical coordinates enter residuals; there is no HEPData value at 0.01.

The frozen file keeps the demonstration independent of HEPData availability.
