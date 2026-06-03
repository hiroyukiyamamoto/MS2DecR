rm(list = ls(all = TRUE))

library(MASS)
source("matchedfilter.R")
source("generate_ms2_matrix.R")

# Parameters ---------------------------------------------------------------
# Known precursor assumption:
# We skip MS1 peak spotting and start from a user-specified precursor m/z and RT window.
premz <- 290.1387
rt_range <- c(170, 190)
bin_size <- 0.01

# Parameters for the MS2Dec-like downstream procedure
smooth_window <- 5
segment_size <- 8
band_width <- 6
ideal_min <- 0.95
sharp_min <- 0.02
int_min <- 300
filter_fwhm <- 5

# Load demo data -----------------------------------------------------------
load("inst/extdata/HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.rds")

ms2_data <- generate_ms2_matrix(
  dia_data = dia_data,
  premz = premz,
  rt_range = rt_range,
  bin_size = bin_size
)

mz <- ms2_data$mz
Y <- ms2_data$Y
rt <- ms2_data$rt

# Utility functions --------------------------------------------------------
linear_weighted_smooth <- function(x, window = 5) {
  if (window <= 1) {
    return(x)
  }

  half <- floor(window / 2)
  weights <- c(seq_len(half + 1), seq(half, 1))
  weights <- weights / sum(weights)

  out <- numeric(length(x))
  for (i in seq_along(x)) {
    idx <- seq.int(max(1, i - half), min(length(x), i + half))
    w <- weights[(idx - i) + half + 1L]
    w <- w / sum(w)
    out[i] <- sum(x[idx] * w)
  }
  out
}

first_derivative5 <- function(x) {
  n <- length(x)
  out <- numeric(n)
  if (n < 5) {
    return(out)
  }
  for (i in 3:(n - 2)) {
    out[i] <- (x[i - 2] - 8 * x[i - 1] + 8 * x[i + 1] - x[i + 2]) / 12
  }
  out
}

second_derivative5 <- function(x) {
  n <- length(x)
  out <- numeric(n)
  if (n < 5) {
    return(out)
  }
  for (i in 3:(n - 2)) {
    out[i] <- (-x[i - 2] + 16 * x[i - 1] - 30 * x[i] + 16 * x[i + 1] - x[i + 2]) / 12
  }
  out
}

compute_peak_thresholds <- function(x) {
  # Following the paper's peak detection concept, estimate noise thresholds
  # from low-amplitude regions of amplitude change, 1st derivative, and 2nd derivative.
  amp_diff <- abs(diff(x))
  d1 <- abs(first_derivative5(x))
  d2 <- abs(second_derivative5(x))

  get_med <- function(v) {
    vmax <- max(v, na.rm = TRUE)
    if (!is.finite(vmax) || vmax <= 0) {
      return(1e-4)
    }
    med <- median(v[v <= 0.05 * vmax], na.rm = TRUE)
    if (!is.finite(med) || med <= 0) 1e-4 else med
  }

  list(
    AF = get_med(amp_diff),
    FF = get_med(d1),
    SF = get_med(d2)
  )
}

find_peak_bounds <- function(x, apex_idx, threshold) {
  left <- apex_idx
  while (left > 1) {
    if ((x[left] - x[left - 1]) < threshold) {
      break
    }
    left <- left - 1L
  }

  right <- apex_idx
  while (right < length(x)) {
    if ((x[right] - x[right + 1]) < threshold) {
      break
    }
    right <- right + 1L
  }

  c(left = left, right = right)
}

detect_chrom_peaks <- function(x) {
  # Detect chromatographic peaks on a baseline-corrected trace by using
  # derivative sign changes and noise-derived thresholds.
  if (length(x) < 7) {
    return(data.frame())
  }

  thresholds <- compute_peak_thresholds(x)
  d1 <- first_derivative5(x)
  d2 <- second_derivative5(x)

  peaks <- list()
  k <- 1L
  i <- 3L

  while (i <= length(x) - 2L) {
    cond_start <- x[i] > thresholds$AF &&
      d1[i] > thresholds$FF &&
      d1[i + 1L] > thresholds$FF

    if (!cond_start) {
      i <- i + 1L
      next
    }

    apex_candidates <- which(d1[(i + 1L):(length(x) - 2L)] <= 0)
    if (!length(apex_candidates)) {
      break
    }
    apex_idx <- i + apex_candidates[1]
    if (d2[apex_idx] >= -thresholds$SF) {
      i <- apex_idx
      next
    }

    bounds <- find_peak_bounds(x, apex_idx, thresholds$AF)
    left <- bounds[["left"]]
    right <- bounds[["right"]]

    if ((right - left + 1L) >= 5 && x[apex_idx] > 0) {
      peaks[[k]] <- data.frame(
        left_idx = left,
        apex_idx = apex_idx,
        right_idx = right,
        apex_int = x[apex_idx]
      )
      k <- k + 1L
    }
    i <- max(apex_idx + 1L, right)
  }

  if (k == 1L) {
    return(data.frame())
  }

  out <- do.call(rbind, peaks)
  rownames(out) <- NULL
  out
}

baseline_correct <- function(x, segment_size = 8, band_width = 6) {
  # Approximate the paper's baseline correction:
  # split into segments, collect local minima, keep low minima, then connect them.
  n <- length(x)
  idx <- seq_len(n)

  minima_idx <- integer(0)
  for (start in seq(1, n, by = band_width)) {
    stop_idx <- min(n, start + band_width - 1L)
    local_idx <- start:stop_idx
    minima_idx <- c(minima_idx, local_idx[which.min(x[local_idx])])
  }
  minima_idx <- sort(unique(minima_idx))

  keep_idx <- integer(0)
  for (start in seq(1, n, by = segment_size)) {
    stop_idx <- min(n, start + segment_size - 1L)
    seg_idx <- minima_idx[minima_idx >= start & minima_idx <= stop_idx]
    if (!length(seg_idx)) {
      next
    }
    seg_med <- median(x[seg_idx], na.rm = TRUE)
    keep_idx <- c(keep_idx, seg_idx[x[seg_idx] <= seg_med])
  }
  keep_idx <- sort(unique(c(1L, keep_idx, n)))

  baseline <- approx(
    x = keep_idx,
    y = x[keep_idx],
    xout = idx,
    method = "linear",
    ties = "ordered"
  )$y

  corrected <- x - baseline
  corrected[corrected < 0] <- 0

  list(corrected = corrected, baseline = baseline)
}

calc_ideal_slope <- function(x, left_idx, apex_idx, right_idx) {
  # Approximate "ideal slope" by checking monotonic increase to the apex
  # and monotonic decrease after the apex.
  if (apex_idx <= left_idx || apex_idx >= right_idx) {
    return(0)
  }

  left_diff <- diff(x[left_idx:apex_idx])
  right_diff <- diff(x[apex_idx:right_idx])

  left_score <- mean(left_diff >= 0)
  right_score <- mean(right_diff <= 0)

  score <- mean(c(left_score, right_score))
  if (!is.finite(score)) 0 else score
}

calc_sharpness <- function(x, left_idx, apex_idx, right_idx) {
  # Approximate peak sharpness from peak width on the corrected chromatogram.
  peak_height <- x[apex_idx]
  width <- right_idx - left_idx + 1L
  if (!is.finite(peak_height) || peak_height <= 0 || width <= 0) {
    return(0)
  }
  peak_height / width / max(peak_height, 1)
}

select_focused_peak <- function(total_chrom, rt) {
  # The focused precursor peak is determined from the summed MS2 chromatogram,
  # because precursor information is assumed to be known in advance.
  smoothed <- linear_weighted_smooth(total_chrom, smooth_window)
  base <- baseline_correct(smoothed, segment_size = segment_size, band_width = band_width)
  peaks <- detect_chrom_peaks(base$corrected)

  if (!nrow(peaks)) {
    apex_idx <- which.max(base$corrected)
    return(list(
      left_idx = max(1, apex_idx - 3L),
      apex_idx = apex_idx,
      right_idx = min(length(total_chrom), apex_idx + 3L),
      corrected = base$corrected,
      baseline = base$baseline
    ))
  }

  best <- which.max(peaks$apex_int)
  list(
    left_idx = peaks$left_idx[best],
    apex_idx = peaks$apex_idx[best],
    right_idx = peaks$right_idx[best],
    corrected = base$corrected,
    baseline = base$baseline
  )
}

extract_model_peak_candidates <- function(Y, mz, focused_peak,
                                          ideal_min = 0.95,
                                          sharp_min = 0.02,
                                          int_min = 300) {
  # For each fragment ion chromatogram:
  # 1. smooth
  # 2. baseline-correct
  # 3. detect chromatographic peaks
  # 4. keep peaks overlapping the focused precursor region
  # 5. retain only high-quality model peak candidates
  candidates <- list()
  corrected_list <- vector("list", ncol(Y))
  baseline_list <- vector("list", ncol(Y))
  peak_table <- list()
  p <- 1L

  for (j in seq_len(ncol(Y))) {
    raw_x <- Y[, j]
    smoothed <- linear_weighted_smooth(raw_x, smooth_window)
    bc <- baseline_correct(smoothed, segment_size = segment_size, band_width = band_width)
    corrected_list[[j]] <- bc$corrected
    baseline_list[[j]] <- bc$baseline

    peaks <- detect_chrom_peaks(bc$corrected)
    if (!nrow(peaks)) {
      next
    }

    peaks$mz_id <- j
    peaks$mz <- mz[j]
    peaks$ideal_slope <- 0
    peaks$sharpness <- 0
    peaks$overlaps_focus <- FALSE

    for (i in seq_len(nrow(peaks))) {
      peaks$ideal_slope[i] <- calc_ideal_slope(
        bc$corrected,
        peaks$left_idx[i],
        peaks$apex_idx[i],
        peaks$right_idx[i]
      )
      peaks$sharpness[i] <- calc_sharpness(
        bc$corrected,
        peaks$left_idx[i],
        peaks$apex_idx[i],
        peaks$right_idx[i]
      )
      peaks$overlaps_focus[i] <-
        peaks$right_idx[i] >= focused_peak$left_idx &&
        peaks$left_idx[i] <= focused_peak$right_idx
    }

    peak_table[[p]] <- peaks
    p <- p + 1L

    keep <- which(
      peaks$overlaps_focus &
        peaks$apex_int >= int_min &
        peaks$ideal_slope >= ideal_min &
        peaks$sharpness >= sharp_min
    )
    if (!length(keep)) {
      next
    }

    best <- keep[which.max(peaks$apex_int[keep])]
    candidates[[length(candidates) + 1L]] <- peaks[best, ]
  }

  all_peaks <- if (length(peak_table)) do.call(rbind, peak_table) else data.frame()
  candidate_df <- if (length(candidates)) do.call(rbind, candidates) else data.frame()

  list(
    candidates = candidate_df,
    all_peaks = all_peaks,
    corrected = corrected_list,
    baseline = baseline_list
  )
}

build_candidate_sharpness_trace <- function(candidates, n_time) {
  # Store sharpness at each candidate apex scan, then use this trace
  # for the second Gaussian derivative filter.
  out <- numeric(n_time)
  if (!nrow(candidates)) {
    return(out)
  }
  for (i in seq_len(nrow(candidates))) {
    out[candidates$apex_idx[i]] <- out[candidates$apex_idx[i]] + candidates$sharpness[i]
  }
  out
}

find_local_maxima <- function(x) {
  if (length(x) < 3) {
    return(integer(0))
  }
  which(x[2:(length(x) - 1L)] > x[1:(length(x) - 2L)] &
          x[2:(length(x) - 1L)] >= x[3:length(x)]) + 1L
}

choose_model_triplet <- function(candidates, matched_wave, focused_peak, all_peaks) {
  # Pick the focused model peak (M2) and its nearest left/right neighbors (M1/M3).
  # If M2 is missing, fall back to an ad hoc peak with ideal slope 1 and median sharpness.
  if (!nrow(candidates)) {
    stop("No model peak candidates passed the ideal slope filter.")
  }

  candidate_scan <- candidates$apex_idx
  maxima <- find_local_maxima(matched_wave)
  maxima <- maxima[matched_wave[maxima] > 0]

  if (!length(maxima)) {
    maxima <- candidate_scan
  }

  maxima_to_candidate <- function(idx) {
    candidate_scan[which.min(abs(candidate_scan - idx))]
  }
  maxima_scan <- unique(vapply(maxima, maxima_to_candidate, integer(1)))

  focus_center <- focused_peak$apex_idx
  center_scan <- maxima_scan[which.min(abs(maxima_scan - focus_center))]

  candidate_idx <- function(scan_idx) which.min(abs(candidate_scan - scan_idx))
  center_pos <- candidate_idx(center_scan)

  left_pos <- which(candidate_scan < candidate_scan[center_pos])
  right_pos <- which(candidate_scan > candidate_scan[center_pos])

  M1 <- if (length(left_pos)) candidates[left_pos[which.max(candidate_scan[left_pos])], ] else NULL
  M2 <- candidates[center_pos, ]
  M3 <- if (length(right_pos)) candidates[right_pos[which.min(candidate_scan[right_pos])], ] else NULL

  if (!M2$overlaps_focus) {
    focus_candidates <- all_peaks[
      all_peaks$overlaps_focus & abs(all_peaks$ideal_slope - 1) < 1e-8,
      ,
      drop = FALSE
    ]
    if (nrow(focus_candidates)) {
      med_sharp <- median(focus_candidates$sharpness, na.rm = TRUE)
      M2 <- focus_candidates[which.min(abs(focus_candidates$sharpness - med_sharp)), ]
    }
  }

  list(M1 = M1, M2 = M2, M3 = M3)
}

build_model_vector <- function(corrected_list, peak_row, n_time) {
  # Construct a model chromatogram M(n) using only the peak region
  # and normalize it to focus on shape during regression.
  vec <- numeric(n_time)
  if (is.null(peak_row) || !nrow(peak_row)) {
    return(vec)
  }
  idx <- peak_row$left_idx:peak_row$right_idx
  y <- corrected_list[[peak_row$mz_id]][idx]
  if (max(y, na.rm = TRUE) > 0) {
    vec[idx] <- y / max(y, na.rm = TRUE)
  }
  vec
}

deconvolve_ms2dec <- function(Y, model_triplet, corrected_list, focused_peak) {
  # Fit each MS/MS chromatogram by least squares using model peaks and
  # linear baseline terms, then reconstruct component-wise spectra.
  models <- list(model_triplet$M1, model_triplet$M2, model_triplet$M3)
  model_names <- c("M1", "M2", "M3")
  keep <- vapply(models, function(x) !is.null(x) && nrow(x) > 0, logical(1))

  if (!keep[2]) {
    stop("Focused model peak M2 was not assigned.")
  }

  models <- models[keep]
  model_names <- model_names[keep]

  n_time <- nrow(Y)
  model_matrix <- do.call(
    cbind,
    lapply(models, build_model_vector, corrected_list = corrected_list, n_time = n_time)
  )
  colnames(model_matrix) <- model_names

  window <- focused_peak$left_idx:focused_peak$right_idx
  design <- cbind(model_matrix[window, , drop = FALSE], seq_along(window), 1)

  target_names <- colnames(model_matrix)
  spectra_mat <- matrix(0, nrow = length(target_names), ncol = ncol(Y))
  rownames(spectra_mat) <- target_names
  colnames(spectra_mat) <- signif(mz, 5)
  coeff_matrix <- matrix(0, nrow = ncol(design), ncol = ncol(Y))

  for (j in seq_len(ncol(Y))) {
    y <- corrected_list[[j]][window]
    beta <- MASS::ginv(crossprod(design)) %*% crossprod(design, y)
    coeff_matrix[, j] <- beta
    for (nm in target_names) {
      target_beta <- beta[match(nm, colnames(model_matrix))]
      if (is.finite(target_beta) && target_beta > 0) {
        spectra_mat[nm, j] <- target_beta * max(model_matrix[window, nm], na.rm = TRUE)
      }
    }
  }

  rownames(coeff_matrix) <- colnames(design)
  colnames(coeff_matrix) <- signif(mz, 5)

  reconstructed_chrom <- lapply(target_names, function(nm) {
    model_matrix[, nm] %o% pmax(coeff_matrix[nm, ], 0)
  })
  names(reconstructed_chrom) <- target_names

  list(
    spectra = spectra_mat,
    coeff_matrix = coeff_matrix,
    model_matrix = model_matrix,
    window = window,
    reconstructed_chrom = reconstructed_chrom
  )
}

# Run ---------------------------------------------------------------------
# Step 1. Define the focused precursor region from the extracted DIA window.
total_chrom <- rowSums(Y, na.rm = TRUE)
focused_peak <- select_focused_peak(total_chrom, rt)

# Step 2. Extract model peak candidates from individual MS/MS chromatograms.
model_data <- extract_model_peak_candidates(
  Y = Y,
  mz = mz,
  focused_peak = focused_peak,
  ideal_min = ideal_min,
  sharp_min = sharp_min,
  int_min = int_min
)

sharpness_trace <- build_candidate_sharpness_trace(model_data$candidates, nrow(Y))
matched_wave <- matchedfilter(sharpness_trace, filter_fwhm)$filtered_signal

# Step 3. Choose the model peaks used for deconvolution.
model_triplet <- choose_model_triplet(
  candidates = model_data$candidates,
  matched_wave = matched_wave,
  focused_peak = focused_peak,
  all_peaks = model_data$all_peaks
)

# Step 4. Perform least-squares fitting and reconstruct spectra/chromatograms.
deconv <- deconvolve_ms2dec(
  Y = Y,
  model_triplet = model_triplet,
  corrected_list = model_data$corrected,
  focused_peak = focused_peak
)

result <- list(
  focused_peak = focused_peak,
  candidates = model_data$candidates,
  all_peaks = model_data$all_peaks,
  sharpness_trace = sharpness_trace,
  matched_wave = matched_wave,
  model_triplet = model_triplet,
  model_matrix = deconv$model_matrix,
  coeff_matrix = deconv$coeff_matrix,
  spectra = deconv$spectra,
  reconstructed_chrom = deconv$reconstructed_chrom,
  window = deconv$window,
  mz = mz,
  rt = rt
)

cat("Focused peak RT window:", rt[focused_peak$left_idx], "to", rt[focused_peak$right_idx], "\n")
cat("Focused peak apex RT:", rt[focused_peak$apex_idx], "\n")
cat("Model peak candidates:", nrow(model_data$candidates), "\n")

triplet_summary <- do.call(
  rbind,
  lapply(names(model_triplet), function(name) {
    peak <- model_triplet[[name]]
    if (is.null(peak) || !nrow(peak)) {
      return(data.frame(model = name, mz = NA_real_, apex_rt = NA_real_,
                        ideal_slope = NA_real_, sharpness = NA_real_))
    }
    data.frame(
      model = name,
      mz = peak$mz,
      apex_rt = rt[peak$apex_idx],
      ideal_slope = peak$ideal_slope,
      sharpness = peak$sharpness
    )
  })
)
print(triplet_summary, row.names = FALSE)

for (nm in rownames(result$spectra)) {
  top_idx <- order(result$spectra[nm, ], decreasing = TRUE)[seq_len(min(10, ncol(result$spectra)))]
  cat("Top spectrum for", nm, "\n")
  print(data.frame(mz = mz[top_idx], intensity = result$spectra[nm, top_idx]), row.names = FALSE)
}

top_idx_list <- lapply(rownames(result$spectra), function(nm) {
  order(result$spectra[nm, ], decreasing = TRUE)[seq_len(min(10, ncol(result$spectra)))]
})
names(top_idx_list) <- rownames(result$spectra)

png("dev/msdial_like_result.png", width = 1600, height = 1200, res = 150)
op <- par(mfrow = c(2, 3))

plot(rt, total_chrom, type = "l", main = "Total MS2 chromatogram", xlab = "RT", ylab = "Intensity")
abline(v = rt[c(focused_peak$left_idx, focused_peak$apex_idx, focused_peak$right_idx)],
       col = c("grey40", "red", "grey40"), lty = c(2, 1, 2))
plot(rt, focused_peak$corrected, type = "l", main = "Focused peak corrected", xlab = "RT", ylab = "Intensity")
abline(v = rt[c(focused_peak$left_idx, focused_peak$apex_idx, focused_peak$right_idx)],
       col = c("grey40", "red", "grey40"), lty = c(2, 1, 2))
plot(rt, sharpness_trace, type = "h", main = "Candidate sharpness", xlab = "RT", ylab = "Sharpness")
plot(rt, matched_wave, type = "l", main = "Second Gaussian filter", xlab = "RT", ylab = "Matched wave")
matplot(rt, result$model_matrix, type = "l", lty = 1, main = "Model peaks M(n)", xlab = "RT", ylab = "Relative intensity")
if ("M2" %in% rownames(result$spectra)) {
  top_spec <- data.frame(
    mz = mz[top_idx_list[["M2"]]],
    intensity = result$spectra["M2", top_idx_list[["M2"]]]
  )
  plot(top_spec$mz, top_spec$intensity, type = "h", main = "Deconvoluted spectrum (M2)", xlab = "m/z", ylab = "Intensity")
} else {
  plot.new()
}

dev.off()
par(op)

png("dev/msdial_like_chromatograms.png", width = 1600, height = 1000, res = 150)
op <- par(mfrow = c(1, 2))

recon_sum_list <- lapply(names(result$reconstructed_chrom), function(nm) {
  rowSums(result$reconstructed_chrom[[nm]][, top_idx_list[[nm]], drop = FALSE], na.rm = TRUE)
})
names(recon_sum_list) <- names(result$reconstructed_chrom)

matplot(
  rt,
  do.call(cbind, c(as.data.frame(result$model_matrix), lapply(recon_sum_list, identity))),
  type = "l",
  lty = 1,
  lwd = 2,
  col = seq_len(ncol(result$model_matrix) + length(recon_sum_list)),
  main = "Model vs Reconstructed Chromatograms",
  xlab = "RT",
  ylab = "Intensity / Relative intensity"
)
legend(
  "topright",
  legend = c(colnames(result$model_matrix), paste0(names(recon_sum_list), "_recon_top10")),
  col = seq_len(ncol(result$model_matrix) + length(recon_sum_list)),
  lty = 1,
  lwd = 2,
  bty = "n"
)

plot.new()
title("See paired figure for component-wise outputs")

dev.off()
par(op)
png("dev/msdial_like_components.png", width = 1800, height = 900, res = 150)
comp_names <- rownames(result$spectra)
par(mfrow = c(length(comp_names), 2))

for (nm in comp_names) {
  recon_sum <- rowSums(
    result$reconstructed_chrom[[nm]][, top_idx_list[[nm]], drop = FALSE],
    na.rm = TRUE
  )
  plot(
    rt,
    recon_sum,
    type = "l",
    lwd = 2,
    main = paste("Reconstructed chromatogram", nm),
    xlab = "RT",
    ylab = "Intensity"
  )
  lines(rt, result$model_matrix[, nm] * max(recon_sum, na.rm = TRUE), col = "red", lty = 2, lwd = 2)
  legend(
    "topright",
    legend = c(paste0(nm, " recon top10"), paste0(nm, " model shape")),
    col = c("black", "red"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n",
    cex = 0.8
  )

  top_spec <- data.frame(
    mz = mz[top_idx_list[[nm]]],
    intensity = result$spectra[nm, top_idx_list[[nm]]]
  )
  plot(
    top_spec$mz,
    top_spec$intensity,
    type = "h",
    lwd = 2,
    main = paste("Deconvoluted spectrum", nm),
    xlab = "m/z",
    ylab = "Intensity"
  )
}

dev.off()
cat("Saved plots: dev/msdial_like_result.png, dev/msdial_like_chromatograms.png, dev/msdial_like_components.png\n")
