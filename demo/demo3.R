# 必要なライブラリと関数をロード
library(MSnbase)
library(loadings)

source("C:/Users/hyama/Documents/msinfo/deconvICA.R")
source("C:/Users/hyama/Documents/msinfo/performICA.R")
source("C:/Users/hyama/Documents/msinfo/process_ms2_data.R")
source("C:/Users/hyama/Documents/msinfo/match_spectra.R")
source("C:/Users/hyama/Documents/msinfo/detectPeaks.R")
source("C:/Users/hyama/Documents/msinfo/optimizeALS.R")

# DIAファイルの読み込み
dia_file <- "C:/Users/hyama/Documents/msinfo/corrdec/QC1.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

# MS2データをRT範囲でフィルタリング
rt_range <- c(8.8 * 60, 9.1 * 60)
swath_data <- MSnbase::filterMsLevel(dia_data, msLevel = 2L)
swath_data_MS2_RT <- MSnbase::filterRt(swath_data, rt = rt_range)

# MS2データを処理
processed_data <- process_ms2_data(
  swath_data = swath_data,
  peak = list(mz = NULL, rtmin = rt_range[1], rtmax = rt_range[2]),
  bin_size = 0.01,
  com = 5,
  lambda = 0.1,
  maxiter = 50,
  wid = 5,
  gap = 3
)

# デコンボリューション結果の取得
deconv <- processed_data$deconv
mz <- processed_data$mz
Y0 <- deconv$ALS$C

# クロマトグラムのプロット
plot_chromatograms(Y = Y0, mt = swath_data_MS2_RT@featureData@data$retentionTime, ylim = c(0, 1300000))

# デコンボリューションの結果のプロット
par(mfrow = c(deconv$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv$com) {
  plot(mz, deconv$ALS$A[i, ], xlim = c(50, 350), type = "h", main = paste("Component", i))
}

# スペクトルマッチングの計算
load("C:/Users/hyama/Documents/msinfo/corrdec/NAC.RData")  # NACデータ
load("C:/Users/hyama/Documents/msinfo/corrdec/GLN.RData")  # GLNデータ

r1 <- sapply(1:deconv$com, function(i) {
  match_spectra(
    mz = mz,
    intensity = deconv$ALS$A[i, ],
    sp = GLN,
    target_mz = GLN[, 1],
    ppm = 10,
    fun = "dotproduct"
  )
})

r2 <- sapply(1:deconv$com, function(i) {
  match_spectra(
    mz = mz,
    intensity = deconv$ALS$A[i, ],
    sp = NAC,
    target_mz = NAC[, 1],
    ppm = 10,
    fun = "dotproduct"
  )
})

R <- cbind(r1, r2)
print(R)

# 高い相関を持つデータを選択して再デコンボリューション
R_correlation <- apply(Y0, 2, function(y) cor(deconv$ALS$C[, 1], y))
index <- which(R_correlation > 0.9)
Y_selected <- Y0[, index]

deconv2 <- deconvICA(Y_selected, lambda = 0.1, maxiter = 50, wid = 5, gap = 5)

# 再デコンボリューションの結果をプロット
par(mfrow = c(deconv2$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv2$com) {
  plot(mz[index], deconv2$ALS$A[i, ], xlim = c(50, 350), type = "h", main = paste("Component", i))
}

# 再スペクトルマッチングの計算
r1_re <- sapply(1:deconv2$com, function(i) {
  match_spectra(
    mz = mz[index],
    intensity = deconv2$ALS$A[i, ],
    sp = GLN,
    target_mz = GLN[, 1],
    ppm = 10,
    fun = "dotproduct"
  )
})

r2_re <- sapply(1:deconv2$com, function(i) {
  match_spectra(
    mz = mz[index],
    intensity = deconv2$ALS$A[i, ],
    sp = NAC,
    target_mz = NAC[, 1],
    ppm = 10,
    fun = "dotproduct"
  )
})

R_re <- cbind(r1_re, r2_re)
print(R_re)
