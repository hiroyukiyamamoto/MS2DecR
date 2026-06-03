library(MS2DecR)

data(msdial)

res <- deconvICA(msdial$Y, com = 3, lambda = 0.1, maxiter = 20, wid = 5, gap = 1)

image(t(res$ALS$A), xlab = "Component", ylab = "m/z bin")
