rm(list=ls(all=TRUE))

library(MS2DecR)
library(MASS)   # ginv 用

# パラメータ設定
premz   <- 290.1387           # ターゲットの m/z 値
rt_range <- c(170, 190)       # RT の範囲（秒）
bin_size <- 0.01              # ビンサイズ
lambda   <- 0.1               # ALS の正則化パラメータ（ICA側で使用）
maxiter  <- 50                # ALS の最大繰り返し回数（ICA側）
wid      <- 5                 # ALS のウィンドウサイズ（ICA側）
com      <- 5                 # ICA の成分数

# ファイルパスを取得（MS2DecR に付属のデモデータ）
file_path <- system.file(
  "extdata",
  "HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.rds",
  package = "MS2DecR"
)

# データを読み込む（dia_data が入っている想定）
load(file_path)

# MS2データから行列形式を生成
ms2_data <- generate_ms2_matrix(
  dia_data = dia_data,
  premz    = premz,
  rt_range = rt_range,
  bin_size = bin_size
)

mz <- ms2_data$mz
Y  <- ms2_data$Y   # 行=時間, 列=m/z
rt <- ms2_data$rt

detectPeaks2 <- function(Y, rt, int_min = 1000, wid = 5, gap = 3) {
  
  # 変動があり、かつ最大強度が int_min 以上の列だけ対象
  index <- which(apply(Y, 2, sd, na.rm = TRUE) != 0 &
                   apply(Y, 2, max, na.rm = TRUE) >= int_min)
  
  peaks <- list()
  k <- 1L
  
  for (i in index) {
    x <- Y[, i]
    
    # matched filter
    filtered_signal <- matchedfilter(x, wid)[[1]]
    
    # 平滑化前後でピークトップ位置があまりズレていない列だけ採用
    if (abs(which.max(filtered_signal) - which.max(x)) <= gap) {
      
      apex_idx  <- which.max(x)
      apex_int  <- x[apex_idx]
      
      peaks[[k]] <- data.frame(
        mz_id    = i,
        apex_idx = apex_idx,
        apex_rt  = rt[apex_idx],
        apex_int = apex_int
      )
      k <- k + 1L
    }
  }
  
  if (!length(peaks)) {
    return(data.frame())
  }
  
  out <- do.call(rbind, peaks)
  rownames(out) <- NULL
  out
}

# ② ピックされたピークを RT でグルーピングして
#    代表 XIC を model クロマトにする関数
group_ms2_peaks_mf <- function(peaks, Y, rt,
                               rt_tol   = 2,  # [sec] apex RT の差でクラスタリング
                               max_comp = 5) {
  if (is.null(peaks) || nrow(peaks) == 0)
    stop("No peaks to group. 'peaks' is empty.")
  
  # RT で並べ替え
  o     <- order(peaks$apex_rt)
  peaks <- peaks[o, ]
  
  # 隣の apex RT 差が rt_tol を超えるところでグループを切る
  gaps   <- c(Inf, diff(peaks$apex_rt))
  grp_id <- cumsum(gaps > rt_tol)
  peaks$group <- grp_id
  
  groups <- split(seq_len(nrow(peaks)), peaks$group)
  
  # グループの中心 RT で並べ替え
  grp_center <- tapply(peaks$apex_rt, peaks$group, mean)
  ord_grp    <- order(grp_center)
  groups     <- groups[ord_grp]
  
  n_comp  <- min(length(groups), max_comp)
  n_time  <- nrow(Y)
  C       <- matrix(0, nrow = n_time, ncol = n_comp)
  rt_center <- numeric(n_comp)
  
  for (k in seq_len(n_comp)) {
    idxs <- groups[[k]]
    
    # グループ内で apex_int が最大のピークを代表にする
    best_row <- idxs[which.max(peaks$apex_int[idxs])]
    
    j_best   <- peaks$mz_id[best_row]
    i_best   <- peaks$apex_idx[best_row]
    
    xic_best <- Y[, j_best]
    xic_best <- xic_best / max(xic_best, na.rm = TRUE)
    
    C[, k]       <- xic_best
    rt_center[k] <- rt[i_best]
  }
  
  list(
    C         = C,          # time × n_comp のモデルクロマト
    rt_center = rt_center,  # 各コンポーネントの中心 RT
    peaks     = peaks       # group 情報付きのピーク一覧（確認用）
  )
}

peaks <- detectPeaks2(Y,rt)

# 次に、RT でグルーピングしてモデルクロマトを作る
model <- group_ms2_peaks_mf(
  peaks, Y, rt,
  rt_tol   = 0.5,   # MT の近さ（秒）
  max_comp = 2
)

plot(model$C[,1], type="b")
plot(model$C[,2], type="b")