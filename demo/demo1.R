# ---------
#   ICA
# ---------
library(MSnbase)
library(loadings)

source("C:/Users/hyama/Documents/msinfo/dev/deconvICA.R")
source("C:/Users/hyama/Documents/msinfo/dev/detectPeaks.R")
source("C:/Users/hyama/Documents/msinfo/dev/matchedfilter.R")
source("C:/Users/hyama/Documents/msinfo/dev/optimizeALS.R")
source("C:/Users/hyama/Documents/msinfo/dev/performICA.R")

dia_file <- "C:/Users/hyama/Documents/msinfo/HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

premz <- 290.1387
swath_data <- MSnbase::filterMsLevel(dia_data, msLevel=2L)
swath_data <- MSnbase::filterIsolationWindow(swath_data, mz=premz)
swath_data_MS2_RT <- MSnbase::filterRt(swath_data, rt=c(170,190))
# MS1で時間幅を決める

y <- MSnbase::bin(swath_data_MS2_RT, binSize=0.01)

mz <- y[[1]]@mz

Y <- NULL;mt <- NULL
for(i in 1:length(swath_data_MS2_RT)){
  Y <- rbind(Y,y[[i]]@intensity)
  mt <- c(mt,y[[i]]@rt)
}
Y0 <- Y

### デコンボリューション
deconv <- deconvICA(Y0, wid=3)

### 結果の可視化
# データのクロマトグラムを重ね合わせ
plot(mt,Y[,1],type="l",xlim=c(170,190), ylim=c(0,6650))
for(i in 1:ncol(Y)){
 if(max(Y[,i])>100){
   par(new=T)
   plot(mt,Y[,i],type="l",xlim=c(170,190), ylim=c(0,6650), xaxt="n", yaxt="n")
 }
}
# デコンボリューションの結果
# ICAその1, ICAその2, ALS, RALSの4つのパターン

# ICAその1
com <- 5
icacom1 <- ica(t(Y0),com)

par(mfrow=c(2,3), mar=c(2,2,2,2))
for(i in 1:com){
  plot(mt,icacom1$M[,i], type="l")
}

# ICAその2
icacom2 <- ica(Y0,com)

par(mfrow=c(2,3), mar=c(2,2,2,2))
for(i in 1:com){
  plot(mt,icacom2$S[,i], type="l")
}

# ALS
deconv <- deconvICA(Y0, lambda=0, maxiter=50, wid=0.5)

par(mfrow=c(2,3), mar=c(2,2,2,2))
for(i in 1:com){
  plot(mt,deconv$ALS$C[,i], type="l")
}

# RALS
deconv1 <- deconvICA(Y0, lambda=0.1, maxiter=50, wid=0.5)

par(mfrow=c(2,3), mar=c(2,2,2,2))
for(i in 1:deconv$com){
  plot(mt,deconv$ALS$C[,i], type="l")
}

# 最終結果
deconv2 <- deconvICA(Y0, lambda=0.1, maxiter=50, wid=5, gap=1)

par(mfrow=c(2,3), mar=c(2,2,2,2))
for(i in 1:deconv2$com){
  plot(mt,deconv2$ALS$C[,i], type="l")
}

# スペクトル(ALS)
deconv3 <- deconvICA(Y0, lambda=0, maxiter=50, wid=5)

par(mfrow=c(2,1), mar=c(2,2,2,2))
for(i in 1:deconv3$com){
  plot(mt,deconv3$ALS$C[,i], type="l")
}

par(mfrow=c(2,1), mar=c(2,2,2,2))
for(i in 1:deconv3$com){
  plot(mz,deconv3$ALS$A[i,], xlim=c(50,300),type="h")
}

# スペクトル (RALS)
deconv4 <- deconvICA(Y0, lambda=0.1, maxiter=50, wid=5)

par(mfrow=c(2,1), mar=c(2,2,2,2))
for(i in 1:deconv4$com){
  plot(mt,deconv4$ALS$C[,i], type="l")
}

par(mfrow=c(2,1), mar=c(2,2,2,2))
for(i in 1:deconv4$com){
  plot(mz,deconv4$ALS$A[i,], xlim=c(50,300),type="h")
}

# # --------------------------------------------------------
#
# # -------------------------
# #
# # library(MSnbase)
# # library(MsBackendMsp)
# #
# # file_msp <- "C:/Users/hyama/Documents/msinfo/MSMS_Public_EXP_Pos_VS17.msp"
# # sp <- Spectra(file_msp, source = MsBackendMsp())
# #
# # index1 <- which(sp@backend@spectraData$name=="Metoclopramide")
# # index2 <- which(sp@backend@spectraData$name=="Norcocaine")
# #
# # # # -------------------------------------
# # ### Metoclopramide
# #
# R <- NULL
#
# mzmax <- NULL
# for(k in 1:length(index1)){
#   mz_msp_j <- sp@backend@spectraData$mz[[index1[k]]]
#   int_msp_j <- sp@backend@spectraData$intensity[[index1[k]]]
#   mzmax[k] <- mz_msp_j[which.max(int_msp_j)]
# }
#
#
# #sp
#
# r <- NULL
# for(i in 1:deconv$com){
#
#   mz_msp_i <- mz
#   S <- A[i,]
#
#   int_msp_i <- A[i,]
#   int_msp_i <- S
#
#   mz_msp_j <- sp@backend@spectraData$mz[[k]]
#   int_msp_j <- sp@backend@spectraData$intensity[[k]]
#
#   s1 <- new("Spectrum2", mz=mz_msp_i, intensity=int_msp_i)
#   s2 <- new("Spectrum2", mz=mz_msp_j, intensity=int_msp_j)
#
#   r <- compareSpectra(s1, s2, fun="dotproduct") # 類似度尺度の検討
#   print(r)
# }
#
# plot(mz_msp_j,int_msp_j, type="h")
#
#
# #
# # # k=3,4
# #
# # # Norcocaine
# k <- 2
# for(i in 1:deconv$com){
#
#    mz_msp_i <- mz
#    S <- A[i,]
#
#    int_msp_i <- A[i,]
#    int_msp_i <- S
#
#    mz_msp_j <- sp@backend@spectraData$mz[[index2[k]]]
#    int_msp_j <- sp@backend@spectraData$intensity[[index2[k]]]
#
#    s1 <- new("Spectrum2", mz=mz_msp_i, intensity=int_msp_i)
#    s2 <- new("Spectrum2", mz=mz_msp_j, intensity=int_msp_j)
#
#    r <- compareSpectra(s1, s2, fun="dotproduct") # 類似度尺度の検討
#    print(r)
#
#  }
#
# # i=3の時に、0.4285269
#
# ### データのクロマトグラムを重ね合わせ
# # plot(mt/60,Y[,1],type="l",xlim=c(2.83,3.16), ylim=c(0,6650))
# # for(i in 1:ncol(Y)){
# #   if(max(Y[,i])>100){
# #     par(new=T)
# #     plot(mt/60,Y[,i],type="l",xlim=c(2.83,3.16), ylim=c(0,6650), xaxt="n", yaxt="n")
# #   }
# # }
#
# ### デコンボリューション結果のクロマトグラムを重ね合わせ
# for(i in 1:ncol(deconv$ALS$C)){
#   plot(mt/60,deconv$ALS$C[,i], type="l", xlim=c(2.83,3.16), ylim=c(0,9300))
#   par(new=T)
# }
#
# #
# # plot(mt/60,deconv$ALS$C[,1], type="l", lty=2, xlim=c(2.83,3.16), ylim=c(0,9300))
# # par(new=T)
# # plot(mt/60,deconv$ALS$C[,2], type="l", lwd=2, xlim=c(2.83,3.16), ylim=c(0,9300), xaxt="n", yaxt="n")
# # par(new=T)
# # plot(mt/60,deconv$ALS$C[,3], type="l", lwd=2, col="red", xlim=c(2.83,3.16), ylim=c(0,9300), xaxt="n", yaxt="n")
# # par(new=T)
# # plot(mt/60,deconv$ALS$C[,4], type="l", lty=2, xlim=c(2.83,3.16), ylim=c(0,9300), xaxt="n", yaxt="n")
# # par(new=T)
# # plot(mt/60,deconv$ALS$C[,5], type="l", lty=2, xlim=c(2.83,3.16), ylim=c(0,9300), xaxt="n", yaxt="n")
#
#
