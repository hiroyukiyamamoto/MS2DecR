# diaPASEF common-A deconvolution note

This note records the small Fig.6 test implemented in
`dev/test_diapasef_deconv.R`.

## Model

For the extracted diaPASEF demo matrices, the fragment columns are shared:

- `X`: RT x fragment chromatogram matrix
- `Y`: IM x fragment mobilogram matrix
- `C`: model chromatograms
- `M`: model mobilograms
- `A`: shared component spectra

With fixed `C` and `M`, the shared spectrum estimate is the least-squares
solution:

```text
argmin_A ||X - C A||_F^2 + ||Y - M A||_F^2

A = solve(C'C + M'M + lambda I) (C'X + M'Y)
```

This is the direct extension of the MS-DIAL-like model chromatogram fitting
idea to chromatogram and mobilogram matrices that share the same fragment
spectrum `A`.

## Demo Data And Expected Labels

The local demo data are the diaTracer Fig.6 phosphopeptide example:

- `S[167]PSPPDGSPAATPEIR`
- `SPSPPDGS[167]PAATPEIR`

The extracted fragment matrix columns are labelled as `S1_phospho` and
`S8_phospho`, and `data/fig6_opentimsr_minimal/.../fig6_site_localizing_fragments.tsv`
contains the site-localizing fragment m/z list used as the practical answer key
for this test.

## Results From Current Script

Command:

```powershell
& 'C:\Program Files\R\R-4.4.0\bin\Rscript.exe' MS2DecR-main\dev\test_diapasef_deconv.R
```

Outputs are written to:

```text
MS2DecR-main/dev/diapasef_deconv_test/
```

Current summary:

- MS-DIAL-like fixed `C/M` common-A fitting assigns the two components to the
  expected `S1_phospho` and `S8_phospho` groups.
- The site purity is about `0.797` for S1 and `0.810` for S8.
- Joint non-negative MCR-ALS lowers the residual from `58063.17` to `33874.18`,
  but the spectra are more mixed: about `0.745` purity for S1 and `0.633` for
  S8.
- Starting MCR-ALS from the model chromatograms/mobilograms reaches the same
  solution, so this is not only an initialization issue in this small test.

## Interpretation

The common-A least-squares equation is valid and works well when informative
model chromatograms and mobilograms are fixed. In this Fig.6 matrix, however,
unconstrained non-negative MCR-ALS optimizes reconstruction error more strongly
than site-localizing spectral purity. Several fragments around m/z 780-782 are
nearly shared between the two phospho-site labels, so a low-residual ALS
solution can mix those signals.

For this demo, the better practical result is therefore the MS-DIAL-like route:
extract or select reliable model chromatograms and mobilograms first, then solve
the shared-spectrum least-squares problem.

Potential next improvements for MCR-ALS:

- add retention-time and ion-mobility peak-shape constraints;
- anchor selected site-localizing fragments to components;
- add sparsity or selectivity penalties on `A`;
- use explicit component apex ordering in RT and IM.
