library(MSnbase)
library(loadings)
library(xcms)
library(MsBackendMsp)

# --- データの準備 ---
# DIAファイルの読み込み
dia_file <- "C:/Users/hyama/Documents/msinfo/HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.mzML"
swath_data <- readMSData(dia_file, mode = "onDisk")

# MS/MSスペクトルライブラリ読み込み
file_msp <- "C:/Users/hyama/Documents/msinfo/MSMS_Public_EXP_Pos_VS17.msp"
sp <- Spectra(file_msp, source = MsBackendMsp())

# --- 関数の呼び出し ---
results <- process_ms_data(swath_data, sp)

# --- 結果の可視化 ---
hist(results$R[!is.na(results$R)], breaks = 50, main = "R Histogram", xlab = "Similarity (R)")
hist(results$Q[!is.na(results$Q)], breaks = 50, main = "Q Histogram", xlab = "Similarity (Q)")
