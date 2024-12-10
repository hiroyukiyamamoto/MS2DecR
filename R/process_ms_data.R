process_ms_data <- function(swath_data, sp, bin_size = 0.01, com = 5, lambda = 0.1, maxiter = 50, wid = 5, gap = 5, snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3, 30)) {
  
  # MS1ピッキング
  cwp <- CentWaveParam(snthresh = snthresh, noise = noise, ppm = ppm, peakwidth = peakwidth)
  swath_data <- findChromPeaks(swath_data, param = cwp)
  peaklist <- chromPeaks(swath_data)
  
  # MS2処理とスペクトルマッチング
  results <- lapply(1:nrow(peaklist), function(i) {
    print(paste("Processing peak:", i))
    
    # ピーク情報
    peak <- peaklist[i, ]
    
    # MS2データ処理とデコンボリューション
    ms2_result <- process_ms2_data(swath_data, peak, bin_size = bin_size, com = com, lambda = lambda, maxiter = maxiter, wid = wid, gap = gap)
    if (is.null(ms2_result)) {
      return(list(R = NA, Q = NA))
    }
    
    # デコンボリューションあり
    deconv <- ms2_result$deconv
    mz <- ms2_result$mz
    Q <- sapply(1:deconv$com, function(com) {
      spec <- deconv$ALS$A[com, ]
      match_spectra(mz[spec != 0], spec[spec != 0], sp, peak["mz"])
    })
    
    # デコンボリューションなし
    index_min <- which.min(abs(swath_data@featureData@data$retentionTime - peak["rt"]))
    MS2_sample <- swath_data[index_min][[1]]
    R <- match_spectra(MS2_sample@mz, MS2_sample@intensity, sp, peak["mz"])
    
    return(list(R = max(R, na.rm = TRUE), Q = max(Q, na.rm = TRUE)))
  })
  
  # 類似度分布を返す
  R_values <- sapply(results, function(res) res$R)
  Q_values <- sapply(results, function(res) res$Q)
  
  return(list(R = R_values, Q = Q_values))
}
