# Match an experimental spectrum against reference spectra near a precursor m/z.
match_spectra <- function(mz, intensity, sp, target_mz, ppm = 10, fun = "dotproduct") {
  if (!fun %in% c("dotproduct", "cosine")) {
    stop("Only 'dotproduct' and 'cosine' are currently supported.")
  }
  if (length(mz) != length(intensity)) {
    stop("mz and intensity must have the same length.")
  }
  if (!is.numeric(mz) || !is.numeric(intensity)) {
    stop("mz and intensity must be numeric vectors.")
  }

  ref_precursor <- Spectra::precursorMz(sp)
  keep <- is.finite(ref_precursor) &
    abs(ref_precursor - target_mz) <= target_mz * ppm * 1e-6
  if (!any(keep)) {
    return(NA_real_)
  }

  sp <- sp[keep]
  ref_mz <- Spectra::mz(sp)
  ref_intensity <- Spectra::intensity(sp)
  scores <- vapply(
    seq_along(ref_mz),
    function(i) .ms2decr_dotproduct(mz, intensity, ref_mz[[i]], ref_intensity[[i]], ppm),
    numeric(1)
  )
  max(scores, na.rm = TRUE)
}

.ms2decr_dotproduct <- function(query_mz, query_intensity, ref_mz, ref_intensity, ppm) {
  query_keep <- is.finite(query_mz) & is.finite(query_intensity) & query_intensity > 0
  ref_keep <- is.finite(ref_mz) & is.finite(ref_intensity) & ref_intensity > 0
  query_mz <- query_mz[query_keep]
  query_intensity <- query_intensity[query_keep]
  ref_mz <- ref_mz[ref_keep]
  ref_intensity <- ref_intensity[ref_keep]
  if (!length(query_mz) || !length(ref_mz)) {
    return(NA_real_)
  }

  aligned_query <- numeric(length(ref_mz))
  for (i in seq_along(ref_mz)) {
    tol <- ref_mz[i] * ppm * 1e-6
    idx <- which(abs(query_mz - ref_mz[i]) <= tol)
    if (length(idx)) {
      aligned_query[i] <- max(query_intensity[idx], na.rm = TRUE)
    }
  }

  denom <- sqrt(sum(aligned_query^2)) * sqrt(sum(ref_intensity^2))
  if (!is.finite(denom) || denom == 0) {
    return(NA_real_)
  }
  sum(aligned_query * ref_intensity) / denom
}
