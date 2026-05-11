process_ms_data <- function(swath_data, sp, bin_size = 0.01, com = 5, lambda = 0.1,
                            maxiter = 50, wid = 5, gap = 5, snthresh = 5,
                            noise = 100, ppm = 10, peakwidth = c(3, 30)) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required for process_ms_data(). Please install xcms.")
  }

  # MS1 peak picking
  cwp <- xcms::CentWaveParam(
    snthresh = snthresh,
    noise = noise,
    ppm = ppm,
    peakwidth = peakwidth
  )
  swath_data <- xcms::findChromPeaks(swath_data, param = cwp)
  peaklist <- xcms::chromPeaks(swath_data)

  # Process each MS1 peak and evaluate raw / deconvoluted spectra.
  results <- lapply(seq_len(nrow(peaklist)), function(i) {
    print(paste("Processing peak:", i))
    peak <- peaklist[i, ]

    ms2_result <- process_ms2_data(
      swath_data,
      peak,
      bin_size = bin_size,
      com = com,
      lambda = lambda,
      maxiter = maxiter,
      wid = wid,
      gap = gap
    )
    if (is.null(ms2_result)) {
      return(list(R = NA, Q = NA))
    }

    deconv <- ms2_result$deconv
    mz <- ms2_result$mz
    Q <- sapply(seq_len(deconv$com), function(com_idx) {
      spec <- deconv$ALS$A[com_idx, ]
      match_spectra(mz[spec != 0], spec[spec != 0], sp, peak["mz"])
    })

    index_min <- which.min(abs(swath_data@featureData@data$retentionTime - peak["rt"]))
    MS2_sample <- swath_data[index_min][[1]]
    R <- match_spectra(MS2_sample@mz, MS2_sample@intensity, sp, peak["mz"])

    list(R = max(R, na.rm = TRUE), Q = max(Q, na.rm = TRUE))
  })

  R_values <- sapply(results, function(res) res$R)
  Q_values <- sapply(results, function(res) res$Q)

  list(R = R_values, Q = Q_values)
}
