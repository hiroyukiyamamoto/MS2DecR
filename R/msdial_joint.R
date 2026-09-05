.msdial_joint_peak_window <- function(trace, apex_idx, frac = 0.2,
                                      min_half_width = 2L,
                                      left_limit = 1L,
                                      right_limit = length(trace)) {
  cutoff <- trace[apex_idx] * frac
  left <- apex_idx
  while (left > left_limit && trace[left] >= cutoff) left <- left - 1L
  right <- apex_idx
  while (right < right_limit && trace[right] >= cutoff) right <- right + 1L

  list(
    left_idx = max(left_limit, min(left + 1L, apex_idx - min_half_width)),
    apex_idx = apex_idx,
    right_idx = min(right_limit, max(right - 1L, apex_idx + min_half_width))
  )
}

.msdial_joint_regions_overlap <- function(candidate, selected) {
  if (!length(selected)) return(FALSE)
  any(vapply(selected, function(region) {
    !(candidate$right_idx < region$left_idx || candidate$left_idx > region$right_idx)
  }, logical(1)))
}

.msdial_joint_detect_axis_peak_regions <- function(axis, trace,
                                                   axis_range = range(axis),
                                                   detect_fwhm = 5L,
                                                   min_rel_height = 0.05,
                                                   window_frac = 0.2,
                                                   min_half_width = 2L,
                                                   max_peaks = 5L) {
  in_range <- which(axis >= min(axis_range) & axis <= max(axis_range))
  if (length(in_range) < 3L) stop("axis_range must include at least 3 points.")

  x <- numeric(length(trace))
  x[in_range] <- trace[in_range]
  filtered <- matchedfilter(x, detect_fwhm)$filtered_signal
  filtered[!is.finite(filtered)] <- 0

  positive <- filtered > 0
  starts <- which(positive & !c(FALSE, positive[-length(positive)]))
  ends <- which(positive & !c(positive[-1L], FALSE))

  raw_regions <- list()
  if (length(starts) && length(ends)) {
    for (i in seq_along(starts)) {
      left_idx <- max(starts[i], min(in_range))
      right_idx <- min(ends[i], max(in_range))
      if (right_idx - left_idx + 1L < 3L) next

      apex_idx <- left_idx + which.max(x[left_idx:right_idx]) - 1L
      raw_regions[[length(raw_regions) + 1L]] <- list(
        left_idx = left_idx,
        apex_idx = apex_idx,
        right_idx = right_idx
      )
    }
  }

  max_height <- max(x[in_range], na.rm = TRUE)
  raw_regions <- raw_regions[vapply(raw_regions, function(region) {
    x[region$apex_idx] >= max_height * min_rel_height
  }, logical(1))]
  raw_regions <- raw_regions[order(vapply(raw_regions, function(region) {
    x[region$apex_idx]
  }, numeric(1)), decreasing = TRUE)]

  selected <- list()
  for (region in raw_regions) {
    rank <- length(selected) + 1L
    frac <- window_frac[min(rank, length(window_frac))]
    narrowed <- .msdial_joint_peak_window(
      x,
      region$apex_idx,
      frac = frac,
      min_half_width = min_half_width,
      left_limit = min(in_range),
      right_limit = max(in_range)
    )
    if (.msdial_joint_regions_overlap(narrowed, selected)) next
    selected[[length(selected) + 1L]] <- narrowed
    if (length(selected) >= max_peaks) break
  }

  if (!length(selected)) {
    return(data.frame(
      rank = integer(),
      left_idx = integer(),
      apex_idx = integer(),
      right_idx = integer(),
      left_value = numeric(),
      apex_value = numeric(),
      right_value = numeric(),
      apex_intensity = numeric(),
      relative_height = numeric(),
      width = numeric()
    ))
  }

  do.call(rbind, lapply(seq_along(selected), function(i) {
    region <- selected[[i]]
    data.frame(
      rank = i,
      left_idx = region$left_idx,
      apex_idx = region$apex_idx,
      right_idx = region$right_idx,
      left_value = axis[region$left_idx],
      apex_value = axis[region$apex_idx],
      right_value = axis[region$right_idx],
      apex_intensity = x[region$apex_idx],
      relative_height = x[region$apex_idx] / max_height,
      width = axis[region$right_idx] - axis[region$left_idx]
    )
  }))
}

.msdial_joint_classify_peaks <- function(rt_peaks, im_peaks) {
  rt_n <- nrow(rt_peaks)
  im_n <- nrow(im_peaks)
  if (rt_n <= 1L && im_n <= 1L) return("single_rt_single_im")
  if (rt_n > 1L && im_n <= 1L) return("multiple_rt_single_im")
  if (rt_n <= 1L && im_n > 1L) return("single_rt_multiple_im")
  "multiple_rt_multiple_im"
}

.msdial_joint_run_axis_models <- function(mat, mz, axis, peak_regions,
                                          ideal_ladder, detect_fwhm,
                                          int_min, sharp_min) {
  lapply(seq_len(nrow(peak_regions)), function(i) {
    focused_peak <- as.list(peak_regions[i, c("left_idx", "apex_idx", "right_idx")])
    focused_peak <- lapply(focused_peak, as.integer)
    errors <- character()

    for (ideal_min in ideal_ladder) {
      obj <- list(Y = mat, mz = mz, rt = axis, focused_peak = focused_peak, result = list())
      obj <- tryCatch(
        msdial_model(
          obj,
          int_min = int_min,
          detect_fwhm = detect_fwhm,
          ideal_min = ideal_min,
          sharp_min = sharp_min
        ),
        error = function(e) {
          errors <<- c(errors, sprintf("ideal_min %.2f: %s", ideal_min, conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(obj)) {
        obj <- msdial_deconv(obj)
        return(list(
          obj = obj,
          axis_apex = axis[focused_peak$apex_idx],
          peak_rank = peak_regions$rank[i],
          ideal_min = ideal_min,
          errors = errors
        ))
      }
    }

    stop(paste(errors, collapse = "\n"))
  })
}

.msdial_joint_center_matrix <- function(fits, axis_prefix) {
  mat <- do.call(cbind, lapply(fits, function(fit) {
    fit$obj$result$deconv$model_matrix[, "M2"]
  }))
  colnames(mat) <- sprintf(
    "%s peak %d %.4f",
    axis_prefix,
    vapply(fits, `[[`, numeric(1), "peak_rank"),
    vapply(fits, `[[`, numeric(1), "axis_apex")
  )
  mat
}

.msdial_joint_candidate_matrix <- function(fits, axis, axis_prefix,
                                           model_names = c("M1", "M2", "M3")) {
  pieces <- list()
  for (fit in fits) {
    available <- intersect(model_names, colnames(fit$obj$result$deconv$model_matrix))
    for (model_name in available) {
      model <- fit$obj$result$model$model_triplet[[model_name]]
      col_name <- sprintf(
        "%s peak %d %s %.4f",
        axis_prefix,
        fit$peak_rank,
        model_name,
        axis[model$apex_idx]
      )
      pieces[[col_name]] <- fit$obj$result$deconv$model_matrix[, model_name]
    }
  }
  mat <- do.call(cbind, pieces)
  colnames(mat) <- names(pieces)
  mat
}

.msdial_joint_center_spectra <- function(fits, axis_prefix) {
  spectra <- do.call(rbind, lapply(fits, function(fit) {
    fit$obj$result$deconv$spectra["M2", ]
  }))
  rownames(spectra) <- sprintf(
    "%s peak %d %.4f",
    axis_prefix,
    vapply(fits, `[[`, numeric(1), "peak_rank"),
    vapply(fits, `[[`, numeric(1), "axis_apex")
  )
  spectra
}

.msdial_joint_candidate_spectra <- function(fits, axis, axis_prefix,
                                            model_names = c("M1", "M2", "M3")) {
  pieces <- list()
  for (fit in fits) {
    available <- intersect(model_names, rownames(fit$obj$result$deconv$spectra))
    for (model_name in available) {
      model <- fit$obj$result$model$model_triplet[[model_name]]
      row_name <- sprintf(
        "%s peak %d %s %.4f",
        axis_prefix,
        fit$peak_rank,
        model_name,
        axis[model$apex_idx]
      )
      pieces[[row_name]] <- fit$obj$result$deconv$spectra[model_name, ]
    }
  }
  spectra <- do.call(rbind, pieces)
  rownames(spectra) <- names(pieces)
  spectra
}

.msdial_joint_axis_summary <- function(fits, axis_prefix, model_names = "M2") {
  rows <- list()
  for (fit in fits) {
    available <- intersect(model_names, colnames(fit$obj$result$deconv$model_matrix))
    for (model_name in available) {
      model <- fit$obj$result$model$model_triplet[[model_name]]
      rows[[length(rows) + 1L]] <- data.frame(
        axis = axis_prefix,
        peak_rank = fit$peak_rank,
        focused_apex = fit$axis_apex,
        model = model_name,
        model_apex_idx = model$apex_idx,
        model_apex = fit$obj$rt[model$apex_idx],
        representative_mz = model$representative_mz,
        n_fragment_peaks = model$n_peaks,
        ideal_min = fit$ideal_min
      )
    }
  }
  do.call(rbind, rows)
}

.msdial_joint_cosine <- function(a, b) {
  a[!is.finite(a)] <- 0
  b[!is.finite(b)] <- 0
  denom <- sqrt(sum(a * a)) * sqrt(sum(b * b))
  if (denom <= 0) return(NA_real_)
  sum(a * b) / denom
}

#' Build MS-DIAL-like RT and ion mobility model profiles
#'
#' \code{msdial_joint_model()} applies the existing \code{\link{msdial_model}}
#' and \code{\link{msdial_deconv}} workflow independently to prepared
#' retention-time and ion-mobility matrices. It returns model chromatograms,
#' model mobilograms, and neighboring \code{M1}/\code{M2}/\code{M3} candidates
#' for later matching.
#'
#' @param X RT x fragment m/z intensity matrix.
#' @param Y Ion mobility x fragment m/z intensity matrix.
#' @param mz Numeric fragment m/z values matching columns of \code{X} and
#'   \code{Y}.
#' @param rt Numeric retention-time axis matching rows of \code{X}.
#' @param im Numeric ion-mobility axis matching rows of \code{Y}.
#' @param rt_range,im_range Axis ranges used for peak-region detection.
#' @param rt_max_peaks,im_max_peaks Maximum number of focused peak regions to
#'   keep on each axis.
#' @param rt_min_rel_height,im_min_rel_height Minimum relative peak height for
#'   focused peak-region detection.
#' @param rt_detect_fwhm,im_detect_fwhm Matched-filter widths for RT and ion
#'   mobility axes.
#' @param rt_ideal_ladder,im_ideal_ladder Candidate \code{ideal_min} values
#'   tried when calling \code{\link{msdial_model}}.
#' @param int_min,sharp_min Passed to \code{\link{msdial_model}}.
#'
#' @return A list containing center model matrices \code{C} and \code{M},
#'   candidate model matrices \code{C_candidates} and \code{M_candidates},
#'   preliminary spectra, diagnostics, and summaries.
#' @examples
#' \dontrun{
#' data(diapasef_phospho_fig6)
#' model <- msdial_joint_model(
#'   diapasef_phospho_fig6$X,
#'   diapasef_phospho_fig6$Y,
#'   diapasef_phospho_fig6$mz,
#'   diapasef_phospho_fig6$rt,
#'   diapasef_phospho_fig6$im,
#'   rt_range = c(17.6, 18.6)
#' )
#' }
msdial_joint_model <- function(X, Y, mz, rt, im,
                               rt_range = range(rt),
                               im_range = range(im),
                               rt_max_peaks = 5L,
                               im_max_peaks = 5L,
                               rt_min_rel_height = 0.02,
                               im_min_rel_height = 0.2,
                               rt_detect_fwhm = 5L,
                               im_detect_fwhm = 3L,
                               rt_ideal_ladder = c(0.95, 0.9, 0.8, 0.7, 0.6, 0.5),
                               im_ideal_ladder = c(0.6, 0.5),
                               int_min = 300,
                               sharp_min = 0.02) {
  stopifnot(nrow(X) == length(rt), nrow(Y) == length(im), ncol(X) == ncol(Y))

  rt_trace <- rowSums(X, na.rm = TRUE)
  im_trace <- rowSums(Y, na.rm = TRUE)
  rt_peaks <- .msdial_joint_detect_axis_peak_regions(
    rt, rt_trace,
    axis_range = rt_range,
    detect_fwhm = rt_detect_fwhm,
    min_rel_height = rt_min_rel_height,
    window_frac = c(0.2, 0.5),
    min_half_width = 2L,
    max_peaks = rt_max_peaks
  )
  im_peaks <- .msdial_joint_detect_axis_peak_regions(
    im, im_trace,
    axis_range = im_range,
    detect_fwhm = im_detect_fwhm,
    min_rel_height = im_min_rel_height,
    window_frac = 0.2,
    min_half_width = 4L,
    max_peaks = im_max_peaks
  )

  if (!nrow(rt_peaks)) stop("No RT peak region was detected.")
  if (!nrow(im_peaks)) stop("No IM peak region was detected.")

  rt_fits <- .msdial_joint_run_axis_models(
    X, mz, rt, rt_peaks,
    ideal_ladder = rt_ideal_ladder,
    detect_fwhm = rt_detect_fwhm,
    int_min = int_min,
    sharp_min = sharp_min
  )
  im_fits <- .msdial_joint_run_axis_models(
    Y, mz, im, im_peaks,
    ideal_ladder = im_ideal_ladder,
    detect_fwhm = im_detect_fwhm,
    int_min = int_min,
    sharp_min = sharp_min
  )

  list(
    C = .msdial_joint_center_matrix(rt_fits, "RT"),
    M = .msdial_joint_center_matrix(im_fits, "IM"),
    C_candidates = .msdial_joint_candidate_matrix(rt_fits, rt, "RT"),
    M_candidates = .msdial_joint_candidate_matrix(im_fits, im, "IM"),
    A_rt = .msdial_joint_center_spectra(rt_fits, "RT"),
    A_im = .msdial_joint_center_spectra(im_fits, "IM"),
    A_rt_candidates = .msdial_joint_candidate_spectra(rt_fits, rt, "RT"),
    A_im_candidates = .msdial_joint_candidate_spectra(im_fits, im, "IM"),
    rt_fits = rt_fits,
    im_fits = im_fits,
    diagnostics = list(
      rt_peaks = rt_peaks,
      im_peaks = im_peaks,
      classification = .msdial_joint_classify_peaks(rt_peaks, im_peaks),
      should_keep_nonmax_peaks = nrow(rt_peaks) > 1L || nrow(im_peaks) > 1L,
      traces = list(rt = rt_trace, im = im_trace)
    ),
    summary = rbind(
      .msdial_joint_axis_summary(rt_fits, "RT", "M2"),
      .msdial_joint_axis_summary(im_fits, "IM", "M2")
    ),
    candidate_summary = rbind(
      .msdial_joint_axis_summary(rt_fits, "RT", c("M1", "M2", "M3")),
      .msdial_joint_axis_summary(im_fits, "IM", c("M1", "M2", "M3"))
    )
  )
}

#' Match RT and ion mobility models by preliminary MS/MS similarity
#'
#' @param model Result from \code{\link{msdial_joint_model}}.
#' @param rt_models Use RT center models or all RT candidates.
#' @param im_models Use ion-mobility center models or all ion-mobility
#'   candidates.
#'
#' @return A list containing aligned \code{C} and \code{M} matrices, matched
#'   preliminary spectra, a similarity matrix, and an alignment table.
#' @examples
#' \dontrun{
#' joint <- msdial_joint(model, rt_models = "center", im_models = "candidates")
#' }
msdial_joint <- function(model,
                         rt_models = c("center", "candidates"),
                         im_models = c("candidates", "center")) {
  rt_models <- match.arg(rt_models)
  im_models <- match.arg(im_models)

  C <- if (rt_models == "center") model$C else model$C_candidates
  M <- if (im_models == "center") model$M else model$M_candidates
  A_rt <- if (rt_models == "center") model$A_rt else model$A_rt_candidates
  A_im <- if (im_models == "center") model$A_im else model$A_im_candidates

  if (ncol(C) != ncol(M)) {
    stop(sprintf(
      "The number of RT model chromatograms (%d) and IM model mobilograms (%d) must match for msdial_joint().",
      ncol(C),
      ncol(M)
    ))
  }

  similarity <- outer(
    seq_len(nrow(A_rt)),
    seq_len(nrow(A_im)),
    Vectorize(function(i, j) .msdial_joint_cosine(A_rt[i, ], A_im[j, ]))
  )
  dimnames(similarity) <- list(rownames(A_rt), rownames(A_im))

  available <- seq_len(ncol(similarity))
  assignment <- integer(nrow(similarity))
  for (i in seq_len(nrow(similarity))) {
    best <- available[which.max(similarity[i, available])]
    assignment[i] <- best
    available <- setdiff(available, best)
  }

  M_aligned <- M[, assignment, drop = FALSE]
  colnames(M_aligned) <- colnames(C)

  list(
    C = C,
    M = M_aligned,
    A_rt = A_rt,
    A_im = A_im[assignment, , drop = FALSE],
    similarity = similarity,
    assignment = assignment,
    alignment = data.frame(
      rt_model = rownames(similarity),
      im_model = colnames(similarity)[assignment],
      cosine_similarity = similarity[cbind(seq_len(nrow(similarity)), assignment)]
    ),
    source_model = model
  )
}

#' Estimate joint MS/MS spectra from matched RT and ion mobility models
#'
#' @param joint Result from \code{\link{msdial_joint}}.
#' @param X RT x fragment m/z intensity matrix.
#' @param Y Ion mobility x fragment m/z intensity matrix.
#' @param lambda Ridge regularization strength.
#'
#' @return A list containing aligned \code{C}, \code{M}, estimated spectra
#'   \code{A}, and fitted \code{X}/\code{Y} matrices.
#' @examples
#' \dontrun{
#' fit <- msdial_joint_deconv(joint, diapasef_phospho_fig6$X, diapasef_phospho_fig6$Y)
#' }
msdial_joint_deconv <- function(joint, X, Y, lambda = 1e-8) {
  if (is.null(joint$C) || is.null(joint$M)) {
    stop("joint must contain C and M.")
  }
  if (ncol(joint$C) != ncol(joint$M)) {
    stop("C and M must have the same number of components for joint deconvolution.")
  }
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same fragment columns.")
  }

  C <- joint$C
  M <- joint$M
  lhs <- crossprod(C) + crossprod(M) + diag(lambda, ncol(C))
  rhs <- crossprod(C, X) + crossprod(M, Y)
  A <- solve(lhs, rhs)
  A[!is.finite(A)] <- 0
  A[A < 0] <- 0
  rownames(A) <- colnames(C)
  colnames(A) <- colnames(X)

  list(C = C, M = M, A = A, fitted_X = C %*% A, fitted_Y = M %*% A)
}
