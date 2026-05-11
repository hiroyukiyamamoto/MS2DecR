# MS2DecR

`MS2DecR` is an R package for processing and deconvolving MS/MS spectra from
DIA/SWATH-style mass spectrometry data. The main workflow is ICA/ALS-based
deconvolution for resolving overlapping MS/MS components. Another important
feature is a simplified MS-DIAL-like deconvolution workflow, which reconstructs
component spectra from overlapping fragment chromatograms using model
chromatograms and least-squares fitting. The package also provides MS2 matrix
generation and spectral matching utilities.

## Features

- **ICA/ALS-based deconvolution** for resolving overlapping MS/MS components
  with an independent component analysis and ALS workflow.
- **Simplified MS-DIAL-like deconvolution** using model chromatograms and
  least-squares fragment fitting for DIA/SWATH MS2 deconvolution.
- **MS2 matrix generation** from DIA/SWATH data by precursor m/z, retention
  time range, and m/z bin size.
- **Spectral matching** against reference spectra with configurable mass
  tolerance and similarity functions.
- **Example data** for testing both the ICA/ALS workflow and the MS-DIAL-like
  workflow without loading a raw `mzML` file.

## Installation

Install the package from a local source directory:

```r
install.packages("devtools")
devtools::install_local("path/to/MS2DecR")
```

`MS2DecR` depends on Bioconductor packages such as `MSnbase`, `xcms`, and
`MsBackendMsp`. If they are not already installed, install them with:

```r
install.packages("BiocManager")
BiocManager::install(c("MSnbase", "xcms", "MsBackendMsp"))
```

## Quick Start

The package includes an `msdial` example dataset that can be used to try the
main ICA/ALS-based deconvolution workflow immediately.

```r
library(MS2DecR)

data(msdial)

res <- deconvICA(
  Y0 = msdial$Y,
  com = 5,
  lambda = 0.1,
  maxiter = 50,
  wid = 5,
  gap = 1
)

component_chromatograms <- res$ALS$C
component_spectra <- res$ALS$A

matplot(
  msdial$rt,
  component_chromatograms,
  type = "l",
  lty = 1,
  xlab = "Retention time",
  ylab = "Intensity"
)
```

## MS-DIAL-like Deconvolution

Alongside the ICA/ALS-based main workflow, `MS2DecR` also provides a simplified
MS-DIAL-like deconvolution workflow. This is another important feature of the
package. It is intended for DIA/SWATH MS2 data where fragment chromatograms
from multiple compounds overlap within the same precursor isolation window.

The MS-DIAL-like workflow follows the same broad idea as MS-DIAL-style
deconvolution:

1. Detect candidate fragment chromatogram peaks.
2. Group fragment peaks by nearby apex scans.
3. Select representative model chromatograms around the focused peak.
4. Fit each fragment chromatogram as a combination of the selected models.
5. Reconstruct deconvoluted component spectra from the fitted coefficients.

The implementation is intentionally lightweight and transparent. It does not
aim to reproduce the full MS-DIAL algorithm exactly. Instead, it makes the
model-chromatogram-based idea available inside an R package so that the model
selection, fitting, and reconstructed spectra can be inspected directly.

The two main functions are:

- `msdial_model()`: detects fragment peak candidates, groups them, and builds
  model chromatograms named `M1`, `M2`, and `M3`.
- `msdial_deconv()`: fits fragment chromatograms using those model
  chromatograms and returns deconvoluted spectra.

The included `msdial` dataset is a precomputed example object prepared for this
workflow. It contains an MS2 intensity matrix, m/z values, retention times, and
a focused peak definition, so the workflow can be tested without loading a raw
`mzML` file.

```r
library(MS2DecR)

data(msdial)

# Step 1: build model chromatograms
msdial <- msdial_model(
  msdial,
  int_min = 300,
  detect_fwhm = 5,
  ideal_min = 0.95,
  sharp_min = 0.02,
  scan_gap = 1L
)

msdial$result$model$triplet_summary

# Step 2: deconvolve fragment spectra using the selected models
msdial <- msdial_deconv(msdial)

msdial$result$deconv$spectra
```

## ICA/ALS Deconvolution

The main deconvolution workflow is implemented by `deconvICA()`. Use this
function when you already have an MS2 intensity matrix with retention-time
scans in rows and m/z bins in columns.

```r
library(MS2DecR)

data(msdial)

res <- deconvICA(
  Y0 = msdial$Y,
  com = 5,
  lambda = 0.1,
  maxiter = 50,
  wid = 5,
  gap = 1
)

dim(res$ALS$C)  # component chromatograms
dim(res$ALS$A)  # component spectra
```

## Processing DIA/SWATH Data

For raw DIA/SWATH data, first read the data with `MSnbase`, then generate an
MS2 matrix for a precursor m/z and retention-time window. The resulting matrix
can be passed to the ICA/ALS workflow.

```r
library(MS2DecR)
library(MSnbase)

dia_data <- readMSData("example.mzML", mode = "onDisk")

ms2_matrix <- generate_ms2_matrix(
  dia_data = dia_data,
  premz = 290.1387,
  rt_range = c(170, 190),
  bin_size = 0.01
)

res <- deconvICA(ms2_matrix$Y, com = 5, lambda = 0.1, maxiter = 50)
```

For a higher-level workflow, use `process_ms2_data()` on a single target peak or
`process_ms_data()` on a peak list. These functions use the ICA/ALS-based
deconvolution workflow internally.

```r
peak <- list(mz = 290.1387, rtmin = 170, rtmax = 190)

processed <- process_ms2_data(
  swath_data = dia_data,
  peak = peak,
  bin_size = 0.01,
  com = 5,
  lambda = 0.1,
  maxiter = 50,
  wid = 5,
  gap = 1
)
```

## Main Functions

- `generate_ms2_matrix()`: filter, bin, and format DIA/SWATH MS2 data.
- `deconvICA()`: perform ICA initialization and ALS-based deconvolution.
- `performICA()`: run ICA for component initialization.
- `optimizeALS()`: refine component profiles and spectra with ALS.
- `detectPeaks()`: detect chromatographic peaks in component profiles.
- `msdial_model()`: build simplified MS-DIAL-like model chromatograms.
- `msdial_deconv()`: deconvolve spectra with the selected model chromatograms.
- `match_spectra()`: calculate similarity to reference library spectra.
- `process_ms2_data()`: process and deconvolve one target MS2 region.
- `process_ms_data()`: run a full workflow over a peak list.

## License

This package is licensed under CC BY-NC-ND 4.0.
