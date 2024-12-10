# 必要なライブラリの読み込み
library(MS2DecR)

# パラメータ設定
premz <- 290.1387           # ターゲットの m/z 値
rt_range <- c(170, 190)     # RT の範囲
bin_size <- 0.01            # ビンサイズ
lambda <- 0.1               # ALS の正則化パラメータ
maxiter <- 50               # ALS の最大繰り返し回数
wid <- 5                    # ALS のウィンドウサイズ
com <- 5                    # ICA の成分数

# ファイルパスを取得
file_path <- system.file("extdata", "HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.rds", package = "MS2DecR")

# データを読み込む
load(file_path)

# MS2データから行列形式を生成
ms2_data <- generate_ms2_matrix(
  dia_data = dia_data,
  premz = premz,
  rt_range = rt_range,
  bin_size = bin_size
)

mz <- ms2_data$mz
Y <- ms2_data$Y
rt <- ms2_data$rt

# デコンボリューションの実行
deconv <- deconvICA(Y, lambda = lambda, maxiter = maxiter, wid = wid, com = com)

# デコンボリューション結果の可視化
# クロマトグラム (C 行列: 成分の時間プロファイル)
par(mfrow = c(deconv$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv$com) {
  plot(mt, deconv$ALS$C[, i], type = "l", main = paste("Component", i))
}

# スペクトル (A 行列: 成分の質量スペクトル)
par(mfrow = c(deconv$com, 1), mar = c(2, 2, 2, 2))
for (i in 1:deconv$com) {
  plot(mz, deconv$ALS$A[i, ], type = "h", xlim = c(50, 300), main = paste("Component", i))
}
