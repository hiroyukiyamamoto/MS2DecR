process_ms2_data <- function(swath_data, peak, bin_size = 0.01, com = 5, lambda = 0.1, maxiter = 50, wid = 5, gap = 5) {
  # Isolation windowでMS2データをフィルタリング
  swath_data_iwindow <- filterIsolationWindow(swath_data, mz = peak["mz"])
  swath_data_MS2_RT <- MSnbase::filterRt(swath_data_iwindow, rt = c(peak["rtmin"], peak["rtmax"]))
  
  # MS/MSデータが十分にない場合はNULLを返す
  if (length(swath_data_MS2_RT) < com) return(NULL)
  
  # MS2データをビン化
  y <- MSnbase::bin(swath_data_MS2_RT, binSize = bin_size)
  mz <- y[[1]]@mz
  Y <- do.call(rbind, lapply(y, function(ys) ys@intensity))
  
  # デコンボリューションを実行
  deconv <- deconvICA(Y, com = com, lambda = lambda, maxiter = maxiter, wid = wid, gap = gap)
  
  # デコンボリューション結果とm/zを返す
  return(list(deconv = deconv, mz = mz))
}
