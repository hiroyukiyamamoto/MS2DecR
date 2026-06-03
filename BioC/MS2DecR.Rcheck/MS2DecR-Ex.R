pkgname <- "MS2DecR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('MS2DecR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("deconvICA")
### * deconvICA

flush(stderr()); flush(stdout())

### Name: deconvICA
### Title: ICA-Based Deconvolution of MS2 Data
### Aliases: deconvICA
### Keywords: ICA deconvolution

### ** Examples

library(MS2DecR)

data(msdial)

# The msdial dataset is a precomputed MS-DIAL-like example derived from
# the HILIC positive SWATH example corresponding to the original MS-DIAL paper.
res <- deconvICA(
  Y0 = msdial$Y,
  com = 5,
  lambda = 0.1,
  maxiter = 50,
  wid = 5,
  gap = 1
)

res$com
dim(res$ALS$C)
dim(res$ALS$A)



cleanEx()
nameEx("detectPeaks")
### * detectPeaks

flush(stderr()); flush(stdout())

### Name: detectPeaks
### Title: Peak Detection in MS Spectral Data
### Aliases: detectPeaks
### Keywords: peak detection MS SWATH

### ** Examples

C <- matrix(runif(100), nrow = 20, ncol = 5)
A <- matrix(runif(50), nrow = 5, ncol = 10)
result <- detectPeaks(C, A, wid = 3, gap = 5)



cleanEx()
nameEx("generate_ms2_matrix")
### * generate_ms2_matrix

flush(stderr()); flush(stdout())

### Name: generate_ms2_matrix
### Title: Generate MS2 Data Matrix
### Aliases: generate_ms2_matrix
### Keywords: MS2 DIA binning spectral

### ** Examples

## Not run: 
##D # Load DIA data
##D dia_data <- MSnbase::readMSData("example.mzML", mode = "onDisk")
##D 
##D # Define parameters
##D premz <- 290.1387
##D rt_range <- c(170, 190)
##D bin_size <- 0.01
##D 
##D # Generate MS2 data matrix
##D ms2_data <- generate_ms2_matrix(dia_data, premz, rt_range, bin_size)
##D 
##D # Access the results
##D mz <- ms2_data$mz
##D Y <- ms2_data$Y
##D rt <- ms2_data$rt
##D 
##D # Inspect the results
##D print(dim(Y))  # Dimensions of the intensity matrix
##D print(mz[1:10])  # First 10 m/z values
##D print(rt[1:10])  # First 10 retention times
## End(Not run)



cleanEx()
nameEx("match_spectra")
### * match_spectra

flush(stderr()); flush(stdout())

### Name: match_spectra
### Title: Match Experimental Spectra Against a Reference Library
### Aliases: match_spectra
### Keywords: spectra matching MS2

### ** Examples

if (requireNamespace("Spectra", quietly = TRUE) &&
    requireNamespace("S4Vectors", quietly = TRUE)) {
sp <- Spectra::Spectra(S4Vectors::DataFrame(
  msLevel = 2L,
  precursorMz = 500.2,
  mz = I(list(c(100.1, 101.1, 102.1))),
  intensity = I(list(c(10, 20, 30)))
))
mz <- c(100.1, 101.1, 102.1)
intensity <- c(10, 20, 30)
target_mz <- 500.2
score <- match_spectra(mz, intensity, sp, target_mz)
score
}



cleanEx()
nameEx("matchedfilter")
### * matchedfilter

flush(stderr()); flush(stdout())

### Name: matchedfilter
### Title: Apply a Gaussian Matched Filter
### Aliases: matchedfilter

### ** Examples

x <- dnorm(seq(-3, 3, length.out = 50))
res <- matchedfilter(x, fwhm = 3)
plot(res$filtered_signal, type = "l")



cleanEx()
nameEx("msdial")
### * msdial

flush(stderr()); flush(stdout())

### Name: msdial
### Title: MS-DIAL-like example data for simplified MS2 deconvolution
### Aliases: msdial
### Keywords: datasets

### ** Examples

data(msdial)
str(msdial, max.level = 1)



cleanEx()
nameEx("msdial_deconv")
### * msdial_deconv

flush(stderr()); flush(stdout())

### Name: msdial_deconv
### Title: Perform simplified MS-DIAL-like MS2 deconvolution
### Aliases: msdial_deconv
### Keywords: deconvolution mass spectrometry

### ** Examples

## Not run: 
##D data(msdial)
##D msdial <- msdial_model(msdial)
##D msdial <- msdial_deconv(msdial)
##D msdial$result$deconv$spectra
## End(Not run)



cleanEx()
nameEx("msdial_model")
### * msdial_model

flush(stderr()); flush(stdout())

### Name: msdial_model
### Title: Build simplified MS-DIAL-like model chromatograms
### Aliases: msdial_model
### Keywords: models mass spectrometry

### ** Examples

## Not run: 
##D data(msdial)
##D msdial <- msdial_model(msdial)
##D msdial$result$model$triplet_summary
## End(Not run)



cleanEx()
nameEx("optimizeALS")
### * optimizeALS

flush(stderr()); flush(stdout())

### Name: optimizeALS
### Title: Alternating Least Squares (ALS) Optimization for Deconvolution
### Aliases: optimizeALS
### Keywords: ALS optimization

### ** Examples

# Example usage:
X <- matrix(runif(5000, 0, 1), nrow = 100, ncol = 50)
C <- matrix(runif(500, 0, 1), nrow = 100, ncol = 5)
result <- optimizeALS(X, C, com = 5, lambda = 0.1, maxiter = 50)



cleanEx()
nameEx("performICA")
### * performICA

flush(stderr()); flush(stdout())

### Name: performICA
### Title: Independent Component Analysis (ICA) for Spectral Data
### Aliases: performICA
### Keywords: ICA dimensionality reduction spectral

### ** Examples

# Example usage:
Y0 <- matrix(runif(5000, 0, 1), nrow = 100, ncol = 50)
result <- performICA(Y0, com = 5)



cleanEx()
nameEx("process_ms2_data")
### * process_ms2_data

flush(stderr()); flush(stdout())

### Name: process_ms2_data
### Title: Process MS2 Data for Deconvolution
### Aliases: process_ms2_data
### Keywords: deconvolution MS2 SWATH

### ** Examples

## Not run: 
##D swath_data <- MSnbase::readMSData("example.mzML", mode = "onDisk")
##D peak <- list(mz = 500.2, rtmin = 300, rtmax = 400)
##D result <- process_ms2_data(swath_data, peak)
## End(Not run)



cleanEx()
nameEx("process_ms_data")
### * process_ms_data

flush(stderr()); flush(stdout())

### Name: process_ms_data
### Title: Full Workflow for Processing MS Data
### Aliases: process_ms_data
### Keywords: MS SWATH processing

### ** Examples

## Not run: 
##D result <- process_ms_data(swath_data, sp)
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
