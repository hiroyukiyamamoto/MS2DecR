# Internal helpers for simplified MS-DIAL-like model selection
.calc_ideal_slope <- function(x, left_idx, apex_idx, right_idx) {
  if (apex_idx <= left_idx || apex_idx >= right_idx) return(0)
  score <- mean(c(
    mean(diff(x[left_idx:apex_idx]) >= 0),
    mean(diff(x[apex_idx:right_idx]) <= 0)
  ))
  if (!is.finite(score)) 0 else score
}

.calc_sharpness <- function(x, left_idx, apex_idx, right_idx) {
  peak_height <- x[apex_idx]
  width <- right_idx - left_idx + 1L
  if (!is.finite(peak_height) || peak_height <= 0 || width <= 0) return(0)
  peak_height / width / max(peak_height, 1)
}

.extract_model_peak_candidates <- function(Y, mz, focused_peak, int_min, detect_fwhm,
                                           ideal_min = 0.95, sharp_min = 0.02) {
  candidates <- list()

  for (j in seq_len(ncol(Y))) {
    chrom <- Y[, j]
    if (length(chrom) < 3L || max(chrom, na.rm = TRUE) < int_min) next

    filtered <- matchedfilter(chrom, detect_fwhm)$filtered_signal
    filtered[!is.finite(filtered)] <- 0
    positive <- filtered > 0
    if (!any(positive)) next

    starts <- which(positive & !c(FALSE, positive[-length(positive)]))
    ends <- which(positive & !c(positive[-1L], FALSE))
    if (!length(starts) || !length(ends)) next
    best_idx <- which.max(vapply(seq_along(starts), function(i) {
      max(chrom[starts[i]:ends[i]], na.rm = TRUE)
    }, numeric(1)))

    left_idx <- starts[best_idx]
    right_idx <- ends[best_idx]
    if (right_idx - left_idx + 1L < 5L) next
    if (right_idx < focused_peak$left_idx || left_idx > focused_peak$right_idx) next

    apex_idx <- left_idx + which.max(chrom[left_idx:right_idx]) - 1L
    if (chrom[apex_idx] < int_min) next

    ideal_slope <- .calc_ideal_slope(chrom, left_idx, apex_idx, right_idx)
    sharpness <- .calc_sharpness(chrom, left_idx, apex_idx, right_idx)
    if (ideal_slope < ideal_min || sharpness < sharp_min) next

    candidates[[length(candidates) + 1L]] <- data.frame(
      mz_id = j,
      mz = mz[j],
      left_idx = left_idx,
      apex_idx = apex_idx,
      right_idx = right_idx,
      apex_int = chrom[apex_idx],
      ideal_slope = ideal_slope,
      sharpness = sharpness
    )
  }

  if (!length(candidates)) return(data.frame())
  do.call(rbind, candidates)
}

.build_model_groups <- function(candidates, scan_gap = 1L) {
  if (!nrow(candidates)) return(data.frame())

  ordered <- candidates[order(candidates$apex_idx), , drop = FALSE]
  group_id <- integer(nrow(ordered))
  group_id[1] <- 1L

  if (nrow(ordered) > 1L) {
    for (i in 2:nrow(ordered)) {
      same_group <- (ordered$apex_idx[i] - ordered$apex_idx[i - 1L]) <= scan_gap
      group_id[i] <- if (same_group) group_id[i - 1L] else group_id[i - 1L] + 1L
    }
  }

  ordered$group_id <- group_id
  ordered
}

.build_model_triplet <- function(Y, grouped_candidates, focused_peak, n_time) {
  if (!nrow(grouped_candidates)) stop("No candidate peaks available.")

  split_groups <- split(grouped_candidates, grouped_candidates$group_id)
  models <- lapply(split_groups, function(group_df) {
    anchor <- group_df[which.max(group_df$apex_int), , drop = FALSE]
    idx <- anchor$left_idx[1]:anchor$right_idx[1]
    vec <- numeric(n_time)
    chrom <- Y[idx, anchor$mz_id[1]]
    if (max(chrom, na.rm = TRUE) > 0) vec[idx] <- chrom / max(chrom, na.rm = TRUE)

    list(
      group_id = anchor$group_id[1],
      left_idx = anchor$left_idx[1],
      apex_idx = anchor$apex_idx[1],
      right_idx = anchor$right_idx[1],
      representative_mz = anchor$mz[1],
      n_peaks = nrow(group_df),
      vector = vec
    )
  })

  region_apex <- vapply(models, function(x) x$apex_idx, integer(1))
  center_idx <- which.min(abs(region_apex - focused_peak$apex_idx))

  list(
    M1 = if (center_idx > 1L) models[[center_idx - 1L]] else NULL,
    M2 = models[[center_idx]],
    M3 = if (center_idx < length(models)) models[[center_idx + 1L]] else NULL
  )
}

#' Build simplified MS-DIAL-like model chromatograms
#'
#' This function performs a simplified MS-DIAL-like model selection step on an
#' \code{msdial} object. It detects one matched-filter-based peak candidate from
#' each fragment chromatogram, groups nearby apex scans, and assigns left,
#' center, and right model chromatograms as \code{M1}, \code{M2}, and
#' \code{M3}. This is intentionally MS-DIAL-like rather than an exact
#' reimplementation of MS-DIAL.
#'
#' The results are stored in \code{msdial$result$model}.
#'
#' @param msdial A list-like object such as \code{data(msdial)} containing
#'   \code{Y}, \code{mz}, \code{rt}, and \code{focused_peak}.
#' @param int_min Minimum apex intensity required for a fragment peak candidate.
#' @param detect_fwhm Matched-filter width used for peak detection within each
#'   fragment chromatogram.
#' @param ideal_min Minimum ideal slope value required for a candidate peak.
#' @param sharp_min Minimum sharpness value required for a candidate peak.
#' @param scan_gap Maximum apex scan difference allowed when merging candidates
#'   into the same group.
#'
#' @return The input \code{msdial} object with \code{msdial$result$model}
#'   added.
#' @examples
#' \dontrun{
#' data(msdial)
#' msdial <- msdial_model(msdial)
#' msdial$result$model$triplet_summary
#' }
msdial_model <- function(msdial,
                         int_min = 300,
                         detect_fwhm = 5,
                         ideal_min = 0.95,
                         sharp_min = 0.02,
                         scan_gap = 1L) {
  if (!is.list(msdial) || is.null(msdial$Y) || is.null(msdial$mz) ||
      is.null(msdial$rt) || is.null(msdial$focused_peak)) {
    stop("msdial must contain Y, mz, rt, and focused_peak.")
  }
  if (is.null(msdial$result)) msdial$result <- list()

  candidates <- .extract_model_peak_candidates(
    Y = msdial$Y,
    mz = msdial$mz,
    focused_peak = msdial$focused_peak,
    int_min = int_min,
    detect_fwhm = detect_fwhm,
    ideal_min = ideal_min,
    sharp_min = sharp_min
  )
  grouped_candidates <- .build_model_groups(candidates, scan_gap = scan_gap)
  model_triplet <- .build_model_triplet(msdial$Y, grouped_candidates, msdial$focused_peak, nrow(msdial$Y))

  triplet_summary <- do.call(rbind, lapply(c("M1", "M2", "M3"), function(name) {
    peak <- model_triplet[[name]]
    if (is.null(peak)) {
      return(data.frame(model = name, mz = NA_real_, apex_rt = NA_real_, n_peaks = NA_integer_))
    }
    data.frame(
      model = name,
      mz = peak$representative_mz,
      apex_rt = msdial$rt[peak$apex_idx],
      n_peaks = peak$n_peaks
    )
  }))

  msdial$result$model <- list(
    params = list(
      int_min = int_min,
      detect_fwhm = detect_fwhm,
      ideal_min = ideal_min,
      sharp_min = sharp_min,
      scan_gap = scan_gap
    ),
    candidates = grouped_candidates,
    triplet_summary = triplet_summary,
    model_triplet = model_triplet
  )

  msdial
}
