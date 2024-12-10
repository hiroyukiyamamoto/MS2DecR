# 必要なライブラリと関数をロード
library(MSnbase)
library(loadings)

source("C:/Users/hyama/Documents/msinfo/deconvICA.R")
source("C:/Users/hyama/Documents/msinfo/process_ms2_data.R")
source("C:/Users/hyama/Documents/msinfo/match_spectra.R")

# DIAファイルの読み込み
dia_file <- "C:/Users/hyama/Documents/msinfo/phosphoproteome/CN20160706_P100_Plate34_PC3_T3_P-0034_H01_acq_01.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

# MS2データをRT範囲でフィルタリング
premz <- 939.68
rt_range <- c(1872, 1914)
swath_data <- MSnbase::filterMsLevel(dia_data, msLevel = 2L)
swath_data <- MSnbase::filterIsolationWindow(swath_data, mz = premz)
swath_data_MS2_RT <- MSnbase::filterRt(swath_data, rt = rt_range)

# MS2データの処理
processed_data <- process_ms2_data(
  swath_data = swath_data,
  peak = list(mz = premz, rtmin = rt_range[1], rtmax = rt_range[2]),
  bin_size = 0.01,
  com = 5,
  lambda = 0.1,
  maxiter = 50,
  wid = 5,
  gap = 5
)

# デコンボリューション結果の取得
deconv <- processed_data$deconv
mz <- processed_data$mz

# クロマトグラムのプロット
par(mfrow = c(deconv$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv$com) {
  plot(
    x = swath_data_MS2_RT@featureData@data$retentionTime,
    y = deconv$ALS$C[, i],
    type = "l",
    main = paste("Component", i),
    xlab = "Retention Time",
    ylab = "Intensity"
  )
}

# スペクトルのプロット
par(mfrow = c(deconv$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv$com) {
  plot(
    x = mz,
    y = deconv$ALS$A[i, ],
    type = "h",
    xlim = c(300, 1500),
    main = paste("Spectrum", i),
    xlab = "m/z",
    ylab = "Intensity"
  )
}

# スペクトルマッチングの実行
file <- "C:/Users/hyama/Documents/msinfo/phosphoproteome/library.csv"
L <- read.csv(file)

match_and_print <- function(deconv, mz, L, isomer, fun = "common") {
  indexL <- which(L$isomer == isomer & L$annot != "None")
  scores <- sapply(1:deconv$com, function(i) {
    match_spectra(
      mz = mz,
      intensity = deconv$ALS$A[i, ],
      sp = L[indexL, ],
      target_mz = L[indexL, 1],
      ppm = 10,
      fun = fun
    )
  })
  print(scores)
}

# S4613マッチング
match_and_print(deconv, mz, L, isomer = "S4613")

# S4618マッチング
match_and_print(deconv, mz, L, isomer = "S4618")

# 高相関データの選択と再デコンボリューション
cor_threshold <- 0.9
R <- apply(deconv$ALS$C, 2, function(y) cor(y, deconv$ALS$C[, 1]))
index <- which(R > cor_threshold)
Y_selected <- deconv$ALS$C[, index]

deconv2 <- deconvICA(Y_selected, lambda = 0.1, maxiter = 50, wid = 5, gap = 5)

# 再デコンボリューションの結果プロット
par(mfrow = c(deconv2$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv2$com) {
  plot(
    x = swath_data_MS2_RT@featureData@data$retentionTime,
    y = deconv2$ALS$C[, i],
    type = "l",
    main = paste("Component", i, "after re-deconvolution"),
    xlab = "Retention Time",
    ylab = "Intensity"
  )
}

par(mfrow = c(deconv2$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv2$com) {
  plot(
    x = mz[index],
    y = deconv2$ALS$A[i, ],
    type = "h",
    xlim = c(300, 1500),
    main = paste("Spectrum", i, "after re-deconvolution"),
    xlab = "m/z",
    ylab = "Intensity"
  )
}

# 再スペクトルマッチング
match_and_print(deconv2, mz[index], L, isomer = "S4613")
match_and_print(deconv2, mz[index], L, isomer = "S4618")
