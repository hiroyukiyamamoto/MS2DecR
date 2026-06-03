library(MS2DecR)

data(msdial)

res <- deconvICA(
  Y0 = msdial$Y,
  com = 3,
  lambda = 0.1,
  maxiter = 20,
  wid = 5,
  gap = 1
)

matplot(msdial$rt, res$ALS$C, type = "l", lty = 1,
        xlab = "Retention time", ylab = "Intensity")
