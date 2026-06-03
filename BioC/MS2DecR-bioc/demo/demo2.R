library(MS2DecR)

data(msdial)

modelled <- msdial_model(
  msdial,
  int_min = 300,
  detect_fwhm = 5,
  ideal_min = 0.95,
  sharp_min = 0.02,
  scan_gap = 1L
)

plot(modelled$result$model$model_triplet$M2$vector, type = "l",
     xlab = "Scan", ylab = "Intensity")
