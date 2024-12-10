# MS2DecR

`MS2DecR` is an R package designed for advanced processing and deconvolution of MS/MS spectral data. By integrating Independent Component Analysis (ICA) and Alternating Least Squares (ALS), it provides tools for resolving complex MS/MS spectra and identifying compounds through spectral matching.

## Features

- **MS2 Data Preprocessing**: Filtering of MS2 data by isolation window, retention time, and collision energy.
- **Spectral Deconvolution**: Resolves overlapping spectra into component signals using ICA and ALS.
- **Spectral Matching**: Matches deconvoluted spectra against spectral libraries for compound identification.
- **Customizable Parameters**: Offers flexibility for fine-tuning processing and deconvolution.

## Installation

You can install `MS2DecR` from source using the `devtools` package:

```r
# Install from source
devtools::install_local("path_to_package_directory")
