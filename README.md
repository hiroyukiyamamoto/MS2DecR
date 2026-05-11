# MS2DecR

`MS2DecR` is an R package for processing and deconvolving MS/MS spectra from
DIA/SWATH-style mass spectrometry data. It provides utilities for building MS2
intensity matrices, resolving overlapping fragment chromatograms with
ICA/ALS-based deconvolution, and comparing deconvoluted spectra with reference
spectral libraries.

## Features

- **MS2 matrix generation** from DIA/SWATH data by precursor m/z, retention
  time range, and m/z bin size.
- **ICA/ALS deconvolution** for resolving overlapping MS/MS components.
- **Simplified MS-DIAL-like workflow** using model chromatograms and
  least-squares fragment fitting.
- **Spectral matching** against reference spectra with configurable mass
  tolerance and similarity functions.
- **Example data** for testing the MS-DIAL-like workflow without loading a raw
  mzML file.

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

The package includes an `msdial` example dataset that can be used to test the
simplified MS-DIAL-like workflow.

```r
library(MS2DecR)

data(msdial)

msdial <- msdial_model(msdial)
msdial$result$model$triplet_summary

msdial <- msdial_deconv(msdial)
deconvoluted_spectra <- msdial$result$deconv$spectra

matplot(
  t(deconvoluted_spectra),
  type = "h",
  lty = 1,
  xlab = "m/z bin",
  ylab = "Intensity"
)
```

## ICA/ALS Deconvolution

Use `deconvICA()` directly when you already have an MS2 intensity matrix with
retention-time scans in rows and m/z bins in columns.

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
MS2 matrix for a precursor m/z and retention-time window.

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
`process_ms_data()` on a peak list.

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
