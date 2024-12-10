# 必要なライブラリと関数をロード
library(MSnbase)
library(loadings)
library(PSMatch)

source("C:/Users/hyama/Documents/msinfo/deconvICA.R")
source("C:/Users/hyama/Documents/msinfo/performICA.R")
source("C:/Users/hyama/Documents/msinfo/process_ms2_data.R")
source("C:/Users/hyama/Documents/msinfo/match_spectra.R")
source("C:/Users/hyama/Documents/msinfo/detectPeaks.R")
source("C:/Users/hyama/Documents/msinfo/optimizeALS.R")

# DIAファイルの読み込み
dia_file <- "C:/Users/hyama/Documents/msinfo/CN20161108_SAM_SPECTER_NB2p81_11.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

# MS2データを特定のm/zとRT範囲でフィルタリングして処理
premz1 <- 391.6932
rt_range <- c(0, 500)

# MS2データ処理
processed_data <- process_ms2_data(
  swath_data = dia_data,
  peak = list(mz = premz1, rtmin = rt_range[1], rtmax = rt_range[2]),
  bin_size = 0.01,
  com = 5,
  lambda = 0.01,
  maxiter = 50,
  wid = 5,
  gap = 3
)

# デコンボリューション結果の取得
deconv <- processed_data$deconv
mz0 <- processed_data$mz

# クロマトグラムのプロット
par(mfrow = c(deconv$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv$com) {
  plot(processed_data$deconv$ALS$C[, i], type = "l", main = paste("Component", i))
}

# スペクトルプロット
par(mfrow = c(deconv$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv$com) {
  plot(mz0, deconv$ALS$A[i, ], xlim = c(100, 1000), type = "h", main = paste("Component", i))
}

# スペクトルマッチングの計算
# GFSASSAR
frag1 <- calculateFragments("GFSASSAR", type = c("b", "y"))
r1 <- sapply(1:deconv$com, function(i) {
  match_spectra(
    mz = mz0,
    intensity = deconv$ALS$A[i, ],
    sp = frag1$mz,
    target_mz = frag1$mz,
    ppm = 10,
    fun = "dotproduct"
  )
})

# GFSANSAR
frag2 <- calculateFragments("GFSANSAR", type = c("b", "y"))
r2 <- sapply(1:deconv$com, function(i) {
  match_spectra(
    mz = mz0,
    intensity = deconv$ALS$A[i, ],
    sp = frag2$mz,
    target_mz = frag2$mz,
    ppm = 10,
    fun = "dotproduct"
  )
})

# 結果の表示
R <- cbind(r1, r2)
print(R)
