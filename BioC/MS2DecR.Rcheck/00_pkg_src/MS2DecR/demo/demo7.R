library(MS2DecR)

if (requireNamespace("Spectra", quietly = TRUE) &&
    requireNamespace("S4Vectors", quietly = TRUE)) {
  ref <- S4Vectors::DataFrame(
    msLevel = 2L,
    precursorMz = 500.2,
    mz = I(list(c(100.1, 101.1, 102.1))),
    intensity = I(list(c(10, 20, 30)))
  )
  sp <- Spectra::Spectra(ref)

  score <- match_spectra(
    mz = c(100.1, 101.1, 102.1),
    intensity = c(12, 18, 33),
    sp = sp,
    target_mz = 500.2
  )
  print(score)
}
