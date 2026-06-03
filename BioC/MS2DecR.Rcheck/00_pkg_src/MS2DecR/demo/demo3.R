library(MS2DecR)

data(msdial)

msdial <- msdial_model(msdial)
msdial <- msdial_deconv(msdial)

spectra <- msdial$result$deconv$spectra
plot(as.numeric(colnames(spectra)), spectra["M2", ], type = "h",
     xlab = "m/z", ylab = "Intensity")
