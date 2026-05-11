.deconvolve_ms2dec <- function(Y, model_triplet, focused_peak, mz) {
  models <- list(model_triplet$M1, model_triplet$M2, model_triplet$M3)
  names(models) <- c("M1", "M2", "M3")
  models <- models[!vapply(models, is.null, logical(1))]
  if (!"M2" %in% names(models)) stop("Focused model M2 was not assigned.")

  model_matrix <- do.call(cbind, lapply(models, `[[`, "vector"))
  colnames(model_matrix) <- names(models)
  window <- focused_peak$left_idx:focused_peak$right_idx
  design <- cbind(model_matrix[window, , drop = FALSE], seq_along(window), 1)

  spectra <- matrix(
    0,
    nrow = ncol(model_matrix),
    ncol = ncol(Y),
    dimnames = list(colnames(model_matrix), formatC(mz, format = "f", digits = 3))
  )

  for (j in seq_len(ncol(Y))) {
    beta <- MASS::ginv(crossprod(design)) %*% crossprod(design, Y[window, j])
    for (nm in colnames(model_matrix)) {
      coef <- beta[match(nm, colnames(design))]
      if (is.finite(coef) && coef > 0) {
        spectra[nm, j] <- coef * max(model_matrix[window, nm], na.rm = TRUE)
      }
    }
  }

  list(spectra = spectra, model_matrix = model_matrix, window = window)
}

#' Perform simplified MS-DIAL-like MS2 deconvolution
#'
#' This function performs a simplified MS-DIAL-like deconvolution step on an
#' \code{msdial} object that already contains model chromatograms in
#' \code{msdial$result$model}. It fits fragment chromatograms by least squares
#' using model chromatograms and simple baseline terms. This is intentionally
#' MS-DIAL-like rather than an exact reimplementation of MS-DIAL.
#'
#' The results are stored in \code{msdial$result$deconv}.
#'
#' @param msdial An object containing \code{Y}, \code{mz}, \code{focused_peak},
#'   and \code{msdial$result$model}.
#'
#' @return The input \code{msdial} object with \code{msdial$result$deconv}
#'   added.
#' @examples
#' \dontrun{
#' data(msdial)
#' msdial <- msdial_model(msdial)
#' msdial <- msdial_deconv(msdial)
#' msdial$result$deconv$spectra
#' }
msdial_deconv <- function(msdial) {
  if (!is.list(msdial) || is.null(msdial$Y) || is.null(msdial$mz) ||
      is.null(msdial$focused_peak) || is.null(msdial$result$model)) {
    stop("msdial must contain Y, mz, focused_peak, and result$model.")
  }

  deconv <- .deconvolve_ms2dec(
    Y = msdial$Y,
    model_triplet = msdial$result$model$model_triplet,
    focused_peak = msdial$focused_peak,
    mz = msdial$mz
  )

  msdial$result$deconv <- list(
    spectra = deconv$spectra,
    model_matrix = deconv$model_matrix,
    window = deconv$window
  )

  msdial
}
