generate_ms2_matrix <- function(dia_data, premz, rt_range, bin_size = 0.01) {
  # MS2データのフィルタリング
  swath_data <- MSnbase::filterMsLevel(dia_data, msLevel = 2L)
  swath_data <- MSnbase::filterIsolationWindow(swath_data, mz = premz)
  swath_data_MS2_RT <- MSnbase::filterRt(swath_data, rt = rt_range)

  # MS2データのビニング
  binned_data <- MSnbase::bin(swath_data_MS2_RT, binSize = bin_size)

  # m/z 軸の取得
  mz <- binned_data[[1]]@mz

  # MS2データを行列形式に整形
  Y <- NULL
  for (i in 1:length(swath_data_MS2_RT)) {
    Y <- rbind(Y, binned_data[[i]]@intensity)
  }

  # Retention Time (RT) 軸の取得
  mt <- swath_data_MS2_RT@featureData@data$retentionTime

  # 結果をリストとして返す
  return(list(mz = mz, Y = Y, rt = mt))
}
