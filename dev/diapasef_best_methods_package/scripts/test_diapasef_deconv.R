options(stringsAsFactors = FALSE)

script_path <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
base_dir <- if (length(script_path)) {
  dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
repo_dir <- normalizePath(file.path(base_dir, "..", ".."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(base_dir, "diapasef_deconv_test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_dir, "MS2DecR-main/R/matchedfilter.R"))
source(file.path(repo_dir, "MS2DecR-main/R/msdial_model.R"))
source(file.path(repo_dir, "MS2DecR-main/R/msdial_deconv.R"))

chrom_file <- file.path(
  repo_dir,
  "data/fig6_ms2_matrix_plot_package/fig6_ms2_matrix_plot_package/data",
  "fig6_ms2_chromatogram_matrix_rt_by_fragment.tsv"
)
mobil_file <- file.path(
  repo_dir,
  "data/fig6_ms2_matrix_plot_package/fig6_ms2_matrix_plot_package/data",
  "fig6_ms2_mobilogram_matrix_im_by_fragment.tsv"
)

read_matrix <- function(path, axis_name) {
  x <- read.delim(path, check.names = FALSE)
  axis <- x[[1]]
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- 0
  list(axis = axis, matrix = mat, axis_name = axis_name)
}

normalize_trace <- function(x) {
  x <- pmax(as.numeric(x), 0)
  s <- sqrt(sum(x^2))
  if (!is.finite(s) || s == 0) return(x)
  x / s
}

cosine <- function(x, y) {
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(den) || den == 0) return(NA_real_)
  sum(x * y) / den
}

site_from_col <- function(x) sub("__.*$", "", x)

make_site_models <- function(mat, sites, top_n = 6) {
  models <- matrix(0, nrow(mat), length(sites), dimnames = list(NULL, sites))
  for (site in sites) {
    idx <- which(site_from_col(colnames(mat)) == site)
    totals <- colSums(mat[, idx, drop = FALSE], na.rm = TRUE)
    idx <- idx[order(totals, decreasing = TRUE)[seq_len(min(top_n, length(idx)))]]
    traces <- apply(mat[, idx, drop = FALSE], 2, normalize_trace)
    models[, site] <- normalize_trace(rowMeans(traces, na.rm = TRUE))
  }
  models
}

smooth_trace <- function(x, width = 5) {
  if (length(x) < width) return(as.numeric(x))
  y <- stats::filter(as.numeric(x), rep(1 / width, width), sides = 2)
  y[is.na(y)] <- x[is.na(y)]
  as.numeric(y)
}

peak_bounds <- function(x, apex, frac = 0.2) {
  threshold <- max(x, na.rm = TRUE) * frac
  left <- apex
  while (left > 1L && x[left - 1L] > threshold) left <- left - 1L
  right <- apex
  while (right < length(x) && x[right + 1L] > threshold) right <- right + 1L
  left:right
}

local_peak_count <- function(x, frac = 0.25) {
  if (length(x) < 3L) return(0L)
  threshold <- max(x, na.rm = TRUE) * frac
  sum(x[-c(1L, length(x))] > x[-c(length(x) - 1L, length(x))] &
        x[-c(1L, 2L)] < x[-c(1L, length(x))] &
        x[-c(1L, length(x))] > threshold)
}

make_single_peak_site_models <- function(mat, sites, window_frac = 0.18) {
  models <- matrix(0, nrow(mat), length(sites), dimnames = list(NULL, sites))
  chosen <- data.frame()
  for (site in sites) {
    idx <- which(site_from_col(colnames(mat)) == site)
    scores <- vapply(idx, function(j) {
      y <- pmax(smooth_trace(mat[, j]), 0)
      total <- sum(y)
      if (!is.finite(total) || total <= 0) return(-Inf)
      apex_height <- max(y)
      n_peaks <- local_peak_count(y)
      apex_height / sqrt(total) / (1 + n_peaks)
    }, numeric(1))
    best <- idx[which.max(scores)]
    y <- pmax(smooth_trace(mat[, best]), 0)
    apex <- which.max(y)
    keep <- peak_bounds(y, apex, frac = window_frac)
    model <- numeric(length(y))
    model[keep] <- y[keep]
    models[, site] <- normalize_trace(model)
    chosen <- rbind(
      chosen,
      data.frame(
        site = site,
        fragment = colnames(mat)[best],
        apex_index = apex,
        score = max(scores),
        stringsAsFactors = FALSE
      )
    )
  }
  list(models = models, chosen = chosen)
}

solve_common_spectrum <- function(X, Y, C, M, lambda = 1e-8) {
  k <- ncol(C)
  lhs <- crossprod(C) + crossprod(M) + lambda * diag(k)
  rhs <- crossprod(C, X) + crossprod(M, Y)
  A <- solve(lhs, rhs)
  A[!is.finite(A)] <- 0
  pmax(A, 0)
}

solve_common_spectrum_window_baseline <- function(X, Y, C, M, lambda = 1e-8) {
  k <- ncol(C)
  rt_window <- which(rowSums(C != 0, na.rm = TRUE) > 0)
  im_window <- which(rowSums(M != 0, na.rm = TRUE) > 0)
  if (!length(rt_window)) rt_window <- seq_len(nrow(C))
  if (!length(im_window)) im_window <- seq_len(nrow(M))

  rt_t <- seq_along(rt_window)
  im_t <- seq_along(im_window)
  rt_t <- as.numeric(scale(rt_t, center = TRUE, scale = TRUE))
  im_t <- as.numeric(scale(im_t, center = TRUE, scale = TRUE))
  rt_t[!is.finite(rt_t)] <- 0
  im_t[!is.finite(im_t)] <- 0

  rt_design <- cbind(
    C[rt_window, , drop = FALSE],
    rt_slope = rt_t,
    rt_intercept = 1,
    im_slope = 0,
    im_intercept = 0
  )
  im_design <- cbind(
    M[im_window, , drop = FALSE],
    rt_slope = 0,
    rt_intercept = 0,
    im_slope = im_t,
    im_intercept = 1
  )
  design <- rbind(rt_design, im_design)

  spectra <- matrix(0, nrow = k, ncol = ncol(X), dimnames = list(colnames(C), colnames(X)))
  nuisance <- matrix(
    0,
    nrow = 4,
    ncol = ncol(X),
    dimnames = list(c("rt_slope", "rt_intercept", "im_slope", "im_intercept"), colnames(X))
  )

  lhs <- crossprod(design) + lambda * diag(ncol(design))
  for (j in seq_len(ncol(X))) {
    y <- c(X[rt_window, j], Y[im_window, j])
    beta <- solve(lhs, crossprod(design, y))
    beta[!is.finite(beta)] <- 0
    spectra[, j] <- pmax(beta[seq_len(k)], 0)
    nuisance[, j] <- beta[k + seq_len(4)]
  }

  list(
    spectra = spectra,
    nuisance = nuisance,
    rt_window = rt_window,
    im_window = im_window,
    design = design
  )
}

joint_error <- function(X, Y, C, M, A) {
  sqrt(sum((X - C %*% A)^2) + sum((Y - M %*% A)^2))
}

site_purity <- function(A, labels) {
  res <- lapply(seq_len(nrow(A)), function(i) {
    total <- sum(A[i, ])
    by_site <- tapply(A[i, ], labels, sum)
    purity <- if (total > 0) max(by_site) / total else NA_real_
    data.frame(
      component = rownames(A)[i],
      assigned_site = names(which.max(by_site)),
      purity = purity,
      S1_fraction = unname(by_site["S1_phospho"] / total),
      S8_fraction = unname(by_site["S8_phospho"] / total),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

top_fragments <- function(A, n = 8) {
  out <- lapply(seq_len(nrow(A)), function(i) {
    idx <- order(A[i, ], decreasing = TRUE)[seq_len(min(n, ncol(A)))]
    data.frame(
      component = rownames(A)[i],
      rank = seq_along(idx),
      fragment = colnames(A)[idx],
      intensity = A[i, idx],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

ica_profiles <- function(X, k) {
  ica_fun <- NULL
  if (requireNamespace("fastICA", quietly = TRUE)) {
    ica_fun <- function(z, n.comp) {
      fit <- fastICA::fastICA(z, n.comp = n.comp)
      list(M = t(fit$A), S = fit$S)
    }
  } else if (exists("ica", mode = "function")) {
    ica_fun <- get("ica", mode = "function")
  }

  if (is.null(ica_fun)) {
    sv <- svd(scale(X, center = TRUE, scale = FALSE), nu = k, nv = k)
    C <- abs(sv$u[, seq_len(k), drop = FALSE])
    for (i in seq_len(k)) C[, i] <- normalize_trace(C[, i])
    return(C)
  }

  ic <- ica_fun(t(X), k)
  C <- matrix(0, nrow(X), k)
  keep <- apply(X, 2, sd) > 0
  for (i in seq_len(k)) {
    s <- ic$M[, i]
    r <- suppressWarnings(cor(X[, keep, drop = FALSE], s))
    r <- r[is.finite(r)]
    if (length(r) && r[which.max(abs(r))] < 0) s <- -s
    C[, i] <- normalize_trace(pmax(s, 0))
  }
  C
}

pair_by_spectra <- function(Ax, Ay) {
  k <- nrow(Ax)
  if (k != 2L) stop("This demo pairing helper currently expects two components.")
  score_same <- cosine(Ax[1, ], Ay[1, ]) + cosine(Ax[2, ], Ay[2, ])
  score_swap <- cosine(Ax[1, ], Ay[2, ]) + cosine(Ax[2, ], Ay[1, ])
  if (is.na(score_swap) || (!is.na(score_same) && score_same >= score_swap)) c(1, 2) else c(2, 1)
}

pair_components_by_spectra <- function(C, M, X, Y, lambda = 1e-6) {
  k <- ncol(C)
  Ax <- pmax(solve(crossprod(C) + lambda * diag(k), crossprod(C, X)), 0)
  Ay <- pmax(solve(crossprod(M) + lambda * diag(k), crossprod(M, Y)), 0)
  score <- outer(seq_len(k), seq_len(k), Vectorize(function(i, j) cosine(Ax[i, ], Ay[j, ])))
  score[!is.finite(score)] <- -Inf

  order_m <- integer(k)
  used <- logical(k)
  for (i in seq_len(k)) {
    available <- which(!used)
    best <- available[which.max(score[i, available])]
    order_m[i] <- best
    used[best] <- TRUE
  }
  order_m
}

joint_mcr_als <- function(X, Y, k = 2, lambda = 1e-6, maxiter = 500, C_init = NULL, M_init = NULL) {
  C <- if (is.null(C_init)) ica_profiles(X, k) else C_init
  M <- if (is.null(M_init)) ica_profiles(Y, k) else M_init
  M <- M[, pair_components_by_spectra(C, M, X, Y, lambda), drop = FALSE]

  Alist <- list()
  Clist <- list()
  Mlist <- list()
  err <- numeric(maxiter)

  for (iter in seq_len(maxiter)) {
    A <- solve_common_spectrum(X, Y, C, M, lambda)
    A[!is.finite(A)] <- 0
    A <- pmax(A, 0)
    for (i in seq_len(k)) {
      s <- sqrt(sum(A[i, ]^2))
      if (is.finite(s) && s > 0) A[i, ] <- A[i, ] / s
    }
    gram <- A %*% t(A) + lambda * diag(k)
    C <- X %*% t(A) %*% solve(gram)
    M <- Y %*% t(A) %*% solve(gram)
    C <- pmax(C, 0)
    M <- pmax(M, 0)

    err[iter] <- joint_error(X, Y, C, M, A)
    Alist[[iter]] <- A
    Clist[[iter]] <- C
    Mlist[[iter]] <- M
  }

  best <- which.min(err)
  A <- Alist[[best]]
  C <- Clist[[best]]
  M <- Mlist[[best]]
  colnames(C) <- colnames(M) <- rownames(A) <- paste0("ALS", seq_len(k))
  list(C = C, M = M, A = A, error = err)
}

profile_peak_score <- function(x,
                               wid = 3,
                               gap = 4,
                               max_local_peaks = 2L,
                               smooth_width = 5,
                               local_peak_frac = 0.25,
                               require_filtered_close = TRUE) {
  x <- pmax(as.numeric(x), 0)
  if (!any(x > 0) || length(x) < 7L) {
    return(list(
      pass = FALSE,
      score = 0,
      apex = NA_integer_,
      filtered_apex = NA_integer_,
      close = FALSE,
      n_peaks = NA_integer_
    ))
  }
  xs <- pmax(smooth_trace(x, width = smooth_width), 0)
  filtered <- matchedfilter(xs, wid)$filtered_signal
  filtered[!is.finite(filtered)] <- 0
  apex <- which.max(xs)
  filtered_apex <- which.max(filtered)
  close <- abs(filtered_apex - apex) <= gap
  n_peaks <- local_peak_count(xs, frac = local_peak_frac)
  sharp <- max(xs, na.rm = TRUE) / max(sum(xs > 0), 1)
  score <- max(xs, na.rm = TRUE) * sharp / (1 + n_peaks)
  list(
    pass = (!require_filtered_close || close) && n_peaks <= max_local_peaks,
    score = score,
    apex = apex,
    filtered_apex = filtered_apex,
    close = close,
    n_peaks = n_peaks
  )
}

select_joint_mcr_peak_components <- function(fit,
                                             min_keep = 2L,
                                             max_keep = 5L,
                                             c_wid = 5,
                                             c_gap = 5,
                                             c_max_local_peaks = 2L,
                                             m_wid = 3,
                                             m_gap = 5,
                                             m_max_local_peaks = 2L,
                                             m_smooth_width = 5,
                                             m_local_peak_frac = 0.25,
                                             m_require_filtered_close = TRUE) {
  k <- nrow(fit$A)
  rows <- lapply(seq_len(k), function(i) {
    c_peak <- profile_peak_score(fit$C[, i], wid = c_wid, gap = c_gap, max_local_peaks = c_max_local_peaks)
    m_peak <- profile_peak_score(
      fit$M[, i],
      wid = m_wid,
      gap = m_gap,
      max_local_peaks = m_max_local_peaks,
      smooth_width = m_smooth_width,
      local_peak_frac = m_local_peak_frac,
      require_filtered_close = m_require_filtered_close
    )
    spec_sum <- sum(fit$A[i, ], na.rm = TRUE)
    data.frame(
      component = rownames(fit$A)[i],
      idx = i,
      c_pass = c_peak$pass,
      m_pass = m_peak$pass,
      c_apex = c_peak$apex,
      m_apex = m_peak$apex,
      c_filtered_apex = c_peak$filtered_apex,
      m_filtered_apex = m_peak$filtered_apex,
      c_close = c_peak$close,
      m_close = m_peak$close,
      c_n_peaks = c_peak$n_peaks,
      m_n_peaks = m_peak$n_peaks,
      score = c_peak$score * m_peak$score * log1p(spec_sum),
      spec_sum = spec_sum,
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  keep <- which(table$c_pass & table$m_pass & table$spec_sum > 0)
  if (length(keep) < min_keep) {
    keep <- order(table$score, decreasing = TRUE)[seq_len(min(min_keep, nrow(table)))]
  }
  keep <- keep[order(table$score[keep], decreasing = TRUE)]
  keep <- keep[seq_len(min(length(keep), max_keep))]

  selected <- list(
    C = fit$C[, keep, drop = FALSE],
    M = fit$M[, keep, drop = FALSE],
    A = fit$A[keep, , drop = FALSE],
    error = fit$error,
    component_table = table,
    keep = keep
  )
  colnames(selected$C) <- colnames(selected$M) <- rownames(selected$A) <- table$component[keep]
  selected
}

summarize_site_purity <- function(method, purity_table, error = NA_real_) {
  site_best <- tapply(purity_table$purity, purity_table$assigned_site, max, na.rm = TRUE)
  data.frame(
    method = method,
    S1_best_purity = unname(if ("S1_phospho" %in% names(site_best)) site_best["S1_phospho"] else NA_real_),
    S8_best_purity = unname(if ("S8_phospho" %in% names(site_best)) site_best["S8_phospho"] else NA_real_),
    mean_best_purity = mean(c(
      if ("S1_phospho" %in% names(site_best)) site_best["S1_phospho"] else NA_real_,
      if ("S8_phospho" %in% names(site_best)) site_best["S8_phospho"] else NA_real_
    ), na.rm = TRUE),
    error = error,
    stringsAsFactors = FALSE
  )
}

plot_profiles <- function(axis, profiles, title, xlab, file) {
  png(file, width = 1200, height = 800, res = 150)
  matplot(axis, profiles, type = "l", lty = 1, lwd = 2, xlab = xlab, ylab = "normalized intensity", main = title)
  legend("topright", legend = colnames(profiles), col = seq_len(ncol(profiles)), lty = 1, lwd = 2, bty = "n")
  invisible(dev.off())
}

fragment_mz <- function(x) {
  as.numeric(sub(".*__mz", "", x))
}

select_focused_peak_simple <- function(x, frac = 0.2, min_half_width = 4L) {
  smoothed <- smooth_trace(x, width = 5)
  apex <- which.max(smoothed)
  threshold <- max(smoothed, na.rm = TRUE) * frac
  left <- apex
  while (left > 1L && smoothed[left - 1L] > threshold) left <- left - 1L
  right <- apex
  while (right < length(smoothed) && smoothed[right + 1L] > threshold) right <- right + 1L
  left <- max(1L, min(left, apex - min_half_width))
  right <- min(length(x), max(right, apex + min_half_width))
  list(left_idx = left, apex_idx = apex, right_idx = right)
}

msdial_triplet_matrix <- function(model_triplet) {
  models <- list(model_triplet$M1, model_triplet$M2, model_triplet$M3)
  names(models) <- c("M1", "M2", "M3")
  models <- models[!vapply(models, is.null, logical(1))]
  mat <- do.call(cbind, lapply(models, `[[`, "vector"))
  colnames(mat) <- names(models)
  mat
}

make_mobility_models_from_triplet <- function(Y, mz, model_triplet, n_mobility) {
  models <- list(model_triplet$M1, model_triplet$M2, model_triplet$M3)
  names(models) <- c("M1", "M2", "M3")
  models <- models[!vapply(models, is.null, logical(1))]
  mat <- matrix(0, n_mobility, length(models), dimnames = list(NULL, names(models)))
  for (nm in names(models)) {
    mz_id <- which.min(abs(mz - models[[nm]]$representative_mz))
    y <- smooth_trace(Y[, mz_id])
    apex <- which.max(y)
    keep <- peak_bounds(y, apex, frac = 0.15)
    model <- numeric(length(y))
    model[keep] <- y[keep]
    mat[, nm] <- normalize_trace(model)
  }
  mat
}

plot_deconv_summary <- function(rt, im, C, M, A, labels, title, file, top_n = 12) {
  mz <- fragment_mz(colnames(A))
  site_cols <- ifelse(labels == "S1_phospho", "#2878b5", "#c23b22")
  comp_cols <- c("#2878b5", "#c23b22", "#4d8f31", "#7a5195", "#b7791f")
  comp_cols <- rep(comp_cols, length.out = nrow(A))

  png(file, width = 1800, height = 900, res = 150)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    invisible(dev.off())
  }, add = TRUE)

  layout(matrix(seq_len(nrow(A) * 3), nrow = nrow(A), byrow = TRUE), widths = c(1, 1, 1.35))
  par(mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  for (i in seq_len(nrow(A))) {
    comp <- rownames(A)[i]
    c_scale <- sum(A[i, ])
    plot(rt, C[, i] * c_scale, type = "l", lwd = 2, col = comp_cols[i],
         xlab = "RT (min)", ylab = "scaled intensity",
         main = paste(comp, "chromatogram"))
    grid(col = "gray90")

    plot(im, M[, i] * c_scale, type = "l", lwd = 2, col = comp_cols[i],
         xlab = "1/K0", ylab = "scaled intensity",
         main = paste(comp, "mobilogram"))
    grid(col = "gray90")

    top <- order(A[i, ], decreasing = TRUE)[seq_len(min(top_n, ncol(A)))]
    bar_cols <- site_cols[top]
    bp <- barplot(A[i, top], names.arg = sprintf("%.1f", mz[top]), las = 2,
                  col = bar_cols, border = NA, cex.names = 0.75,
                  xlab = "fragment m/z", ylab = "estimated intensity",
                  main = paste(comp, "spectrum"))
    legend("topright", legend = c("S1", "S8"), fill = c("#2878b5", "#c23b22"),
           bty = "n", cex = 0.8)
    text(bp, A[i, top], labels = sub("_phospho.*", "", labels[top]),
         pos = 3, cex = 0.65, xpd = NA)
  }
  mtext(title, outer = TRUE, cex = 1.3, font = 2)
}

run_ms2decr_msdial_case <- function(X, Y, mz, rt, focused_peak, ideal_min, labels) {
  obj <- list(Y = X, mz = mz, rt = rt, focused_peak = focused_peak, result = list())
  obj <- msdial_model(
    obj,
    int_min = 300,
    detect_fwhm = 5,
    ideal_min = ideal_min,
    sharp_min = 0.02,
    scan_gap = 1L
  )
  obj <- msdial_deconv(obj)
  C <- msdial_triplet_matrix(obj$result$model$model_triplet)
  M <- make_mobility_models_from_triplet(Y, mz, obj$result$model$model_triplet, nrow(Y))
  A <- solve_common_spectrum(X, Y, C, M)
  rownames(A) <- colnames(C)

  list(
    obj = obj,
    C = C,
    M = M,
    A = A,
    direct_purity = site_purity(obj$result$deconv$spectra, labels),
    direct_top = top_fragments(obj$result$deconv$spectra),
    commonA_purity = site_purity(A, labels),
    commonA_top = top_fragments(A)
  )
}

run_axis_msdial_model <- function(mat, mz, axis, focused_peak, ideal_min, detect_fwhm) {
  obj <- list(Y = mat, mz = mz, rt = axis, focused_peak = focused_peak, result = list())
  msdial_model(
    obj,
    int_min = 300,
    detect_fwhm = detect_fwhm,
    ideal_min = ideal_min,
    sharp_min = 0.02,
    scan_gap = 1L
  )
}

triplet_model_sites <- function(model_triplet, mz, labels) {
  models <- list(model_triplet$M1, model_triplet$M2, model_triplet$M3)
  names(models) <- c("M1", "M2", "M3")
  models <- models[!vapply(models, is.null, logical(1))]
  out <- vapply(models, function(model) {
    mz_id <- which.min(abs(mz - model$representative_mz))
    labels[mz_id]
  }, character(1))
  out
}

align_models_by_site <- function(target_sites, model_matrix, model_sites) {
  aligned <- matrix(
    0,
    nrow = nrow(model_matrix),
    ncol = length(target_sites),
    dimnames = list(NULL, names(target_sites))
  )
  for (i in seq_along(target_sites)) {
    hit <- which(model_sites == target_sites[i])
    if (length(hit)) {
      aligned[, i] <- model_matrix[, hit[1]]
    }
  }
  aligned
}

chrom <- read_matrix(chrom_file, "rt_min")
mobil <- read_matrix(mobil_file, "im_bin")
stopifnot(identical(colnames(chrom$matrix), colnames(mobil$matrix)))

X <- chrom$matrix
Y <- mobil$matrix
labels <- site_from_col(colnames(X))
sites <- c("S1_phospho", "S8_phospho")
mz <- fragment_mz(colnames(X))

C_model <- make_site_models(X, sites)
M_model <- make_site_models(Y, sites)
A_model <- solve_common_spectrum(X, Y, C_model, M_model)
rownames(A_model) <- sites

single_chrom <- make_single_peak_site_models(X, sites)
C_single <- single_chrom$models
A_single_chrom <- solve_common_spectrum(X, Y, C_single, M_model)
rownames(A_single_chrom) <- sites

focused_peak <- select_focused_peak_simple(rowSums(X, na.rm = TRUE))
msdial_case_090 <- run_ms2decr_msdial_case(X, Y, mz, chrom$axis, focused_peak, 0.90, labels)
msdial_obj <- msdial_case_090$obj
C_msdial <- msdial_case_090$C
M_msdial <- msdial_case_090$M
A_msdial <- msdial_case_090$A
msdial_case_070 <- run_ms2decr_msdial_case(X, Y, mz, chrom$axis, focused_peak, 0.70, labels)

focused_mobility_peak <- select_focused_peak_simple(rowSums(Y, na.rm = TRUE), frac = 0.2, min_half_width = 4L)
mobility_model_obj <- run_axis_msdial_model(
  mat = Y,
  mz = mz,
  axis = mobil$axis,
  focused_peak = focused_mobility_peak,
  ideal_min = 0.60,
  detect_fwhm = 3
)
M_axis <- msdial_triplet_matrix(mobility_model_obj$result$model$model_triplet)
rt_model_sites <- triplet_model_sites(msdial_case_070$obj$result$model$model_triplet, mz, labels)
im_model_sites <- triplet_model_sites(mobility_model_obj$result$model$model_triplet, mz, labels)
C_rt070 <- msdial_case_070$C
M_im060 <- align_models_by_site(rt_model_sites, M_axis, im_model_sites)
A_rt070_im060 <- solve_common_spectrum(X, Y, C_rt070, M_im060)
rownames(A_rt070_im060) <- names(rt_model_sites)
fit_rt070_im060_bl <- solve_common_spectrum_window_baseline(X, Y, C_rt070, M_im060)
A_rt070_im060_bl <- fit_rt070_im060_bl$spectra

als <- joint_mcr_als(X, Y, k = 2)
als_model_init <- joint_mcr_als(X, Y, k = 2, C_init = C_model, M_init = M_model)
als_over <- joint_mcr_als(X, Y, k = 5)
als_over_selected <- select_joint_mcr_peak_components(als_over, min_keep = 2L, max_keep = 5L)
als_over_selected_mloose <- select_joint_mcr_peak_components(
  als_over,
  min_keep = 2L,
  max_keep = 5L,
  m_gap = 8,
  m_max_local_peaks = 4L
)
als_over_selected_mshape <- select_joint_mcr_peak_components(
  als_over,
  min_keep = 2L,
  max_keep = 5L,
  m_gap = 8,
  m_max_local_peaks = 6L,
  m_smooth_width = 7,
  m_local_peak_frac = 0.50,
  m_require_filtered_close = FALSE
)
als_over_selected_Aonly <- als_over_selected
als_over_selected_Aonly$A <- solve_common_spectrum_window_baseline(
  X,
  Y,
  als_over_selected$C,
  als_over_selected$M
)$spectra
als_over_selected_mloose_Aonly <- als_over_selected_mloose
als_over_selected_mloose_Aonly$A <- solve_common_spectrum_window_baseline(
  X,
  Y,
  als_over_selected_mloose$C,
  als_over_selected_mloose$M
)$spectra
als_over_selected_mshape_Aonly <- als_over_selected_mshape
als_over_selected_mshape_Aonly$A <- solve_common_spectrum_window_baseline(
  X,
  Y,
  als_over_selected_mshape$C,
  als_over_selected_mshape$M
)$spectra
als_over_selected_reals <- joint_mcr_als(
  X,
  Y,
  k = ncol(als_over_selected$C),
  maxiter = 50,
  C_init = als_over_selected$C,
  M_init = als_over_selected$M
)

model_purity <- site_purity(A_model, labels)
single_chrom_purity <- site_purity(A_single_chrom, labels)
msdial_purity <- site_purity(A_msdial, labels)
msdial_direct_purity <- msdial_case_090$direct_purity
msdial_070_purity <- msdial_case_070$commonA_purity
msdial_070_direct_purity <- msdial_case_070$direct_purity
rt070_im060_purity <- site_purity(A_rt070_im060, labels)
rt070_im060_bl_purity <- site_purity(A_rt070_im060_bl, labels)
als_purity <- site_purity(als$A, labels)
als_model_init_purity <- site_purity(als_model_init$A, labels)
als_over_purity <- site_purity(als_over$A, labels)
als_over_selected_purity <- site_purity(als_over_selected$A, labels)
als_over_selected_Aonly_purity <- site_purity(als_over_selected_Aonly$A, labels)
als_over_selected_mloose_purity <- site_purity(als_over_selected_mloose$A, labels)
als_over_selected_mloose_Aonly_purity <- site_purity(als_over_selected_mloose_Aonly$A, labels)
als_over_selected_mshape_purity <- site_purity(als_over_selected_mshape$A, labels)
als_over_selected_mshape_Aonly_purity <- site_purity(als_over_selected_mshape_Aonly$A, labels)
als_over_selected_reals_purity <- site_purity(als_over_selected_reals$A, labels)
als_strategy_comparison <- do.call(rbind, list(
  summarize_site_purity("A_peak_select_only", als_over_selected_purity, tail(als_over$error, 1)),
  summarize_site_purity("B_fixed_CM_window_baseline_A_only", als_over_selected_Aonly_purity, NA_real_),
  summarize_site_purity("B2_mobility_loose_peak_select_only", als_over_selected_mloose_purity, tail(als_over$error, 1)),
  summarize_site_purity("B3_mobility_loose_fixed_CM_A_only", als_over_selected_mloose_Aonly_purity, NA_real_),
  summarize_site_purity("B4_mobility_shape_peak_select_only", als_over_selected_mshape_purity, tail(als_over$error, 1)),
  summarize_site_purity("B5_mobility_shape_fixed_CM_A_only", als_over_selected_mshape_Aonly_purity, NA_real_),
  summarize_site_purity("C_reALS_after_peak_select", als_over_selected_reals_purity, tail(als_over_selected_reals$error, 1)),
  summarize_site_purity("MSDIAL_like_window_baseline_reference", rt070_im060_bl_purity, NA_real_)
))
model_top <- top_fragments(A_model)
single_chrom_top <- top_fragments(A_single_chrom)
msdial_top <- top_fragments(A_msdial)
msdial_direct_top <- msdial_case_090$direct_top
msdial_070_top <- msdial_case_070$commonA_top
msdial_070_direct_top <- msdial_case_070$direct_top
rt070_im060_top <- top_fragments(A_rt070_im060)
rt070_im060_bl_top <- top_fragments(A_rt070_im060_bl)
als_top <- top_fragments(als$A)
als_model_init_top <- top_fragments(als_model_init$A)
als_over_top <- top_fragments(als_over$A)
als_over_selected_top <- top_fragments(als_over_selected$A)
als_over_selected_Aonly_top <- top_fragments(als_over_selected_Aonly$A)
als_over_selected_mloose_top <- top_fragments(als_over_selected_mloose$A)
als_over_selected_mloose_Aonly_top <- top_fragments(als_over_selected_mloose_Aonly$A)
als_over_selected_mshape_top <- top_fragments(als_over_selected_mshape$A)
als_over_selected_mshape_Aonly_top <- top_fragments(als_over_selected_mshape_Aonly$A)
als_over_selected_reals_top <- top_fragments(als_over_selected_reals$A)

write.csv(model_purity, file.path(out_dir, "msdial_like_commonA_purity.csv"), row.names = FALSE)
write.csv(single_chrom$chosen, file.path(out_dir, "msdial_like_single_peak_chrom_model_fragments.csv"), row.names = FALSE)
write.csv(single_chrom_purity, file.path(out_dir, "msdial_like_single_peak_chrom_commonA_purity.csv"), row.names = FALSE)
write.csv(msdial_obj$result$model$candidates, file.path(out_dir, "ms2decr_msdial_model_candidates.csv"), row.names = FALSE)
write.csv(msdial_obj$result$model$triplet_summary, file.path(out_dir, "ms2decr_msdial_model_triplet_summary.csv"), row.names = FALSE)
write.csv(msdial_direct_purity, file.path(out_dir, "ms2decr_msdial_deconv_direct_purity.csv"), row.names = FALSE)
write.csv(msdial_purity, file.path(out_dir, "ms2decr_msdial_model_commonA_purity.csv"), row.names = FALSE)
write.csv(msdial_case_070$obj$result$model$candidates, file.path(out_dir, "ms2decr_msdial_ideal070_candidates.csv"), row.names = FALSE)
write.csv(msdial_case_070$obj$result$model$triplet_summary, file.path(out_dir, "ms2decr_msdial_ideal070_triplet_summary.csv"), row.names = FALSE)
write.csv(msdial_070_direct_purity, file.path(out_dir, "ms2decr_msdial_ideal070_deconv_direct_purity.csv"), row.names = FALSE)
write.csv(msdial_070_purity, file.path(out_dir, "ms2decr_msdial_ideal070_commonA_purity.csv"), row.names = FALSE)
write.csv(mobility_model_obj$result$model$candidates, file.path(out_dir, "ms2decr_mobilogram_ideal060_fwhm3_candidates.csv"), row.names = FALSE)
write.csv(mobility_model_obj$result$model$triplet_summary, file.path(out_dir, "ms2decr_mobilogram_ideal060_fwhm3_triplet_summary.csv"), row.names = FALSE)
write.csv(data.frame(model = names(rt_model_sites), rt_site = unname(rt_model_sites), im_site = unname(im_model_sites[match(rt_model_sites, im_model_sites)])), file.path(out_dir, "ms2decr_rt070_im060_model_alignment.csv"), row.names = FALSE)
write.csv(rt070_im060_purity, file.path(out_dir, "ms2decr_rt070_im060_commonA_purity.csv"), row.names = FALSE)
write.csv(rt070_im060_bl_purity, file.path(out_dir, "ms2decr_rt070_im060_window_baseline_commonA_purity.csv"), row.names = FALSE)
write.csv(als_purity, file.path(out_dir, "mcr_als_commonA_purity.csv"), row.names = FALSE)
write.csv(als_model_init_purity, file.path(out_dir, "mcr_als_model_init_commonA_purity.csv"), row.names = FALSE)
write.csv(als_over$A, file.path(out_dir, "mcr_als_overcomplete_commonA_spectra.csv"), row.names = TRUE)
write.csv(als_over_purity, file.path(out_dir, "mcr_als_overcomplete_commonA_purity.csv"), row.names = FALSE)
write.csv(als_over_selected$component_table, file.path(out_dir, "mcr_als_overcomplete_component_table.csv"), row.names = FALSE)
write.csv(als_over_selected_purity, file.path(out_dir, "mcr_als_overcomplete_peak_selected_commonA_purity.csv"), row.names = FALSE)
write.csv(als_over_selected_Aonly_purity, file.path(out_dir, "mcr_als_overcomplete_peak_selected_Aonly_purity.csv"), row.names = FALSE)
write.csv(als_over_selected_mloose$component_table, file.path(out_dir, "mcr_als_overcomplete_mobility_loose_component_table.csv"), row.names = FALSE)
write.csv(als_over_selected_mloose_purity, file.path(out_dir, "mcr_als_overcomplete_mobility_loose_peak_selected_purity.csv"), row.names = FALSE)
write.csv(als_over_selected_mloose_Aonly_purity, file.path(out_dir, "mcr_als_overcomplete_mobility_loose_Aonly_purity.csv"), row.names = FALSE)
write.csv(als_over_selected_mshape$component_table, file.path(out_dir, "mcr_als_overcomplete_mobility_shape_component_table.csv"), row.names = FALSE)
write.csv(als_over_selected_mshape_purity, file.path(out_dir, "mcr_als_overcomplete_mobility_shape_peak_selected_purity.csv"), row.names = FALSE)
write.csv(als_over_selected_mshape_Aonly_purity, file.path(out_dir, "mcr_als_overcomplete_mobility_shape_Aonly_purity.csv"), row.names = FALSE)
write.csv(als_over_selected_reals_purity, file.path(out_dir, "mcr_als_overcomplete_peak_selected_reals_purity.csv"), row.names = FALSE)
write.csv(als_strategy_comparison, file.path(out_dir, "mcr_als_strategy_comparison.csv"), row.names = FALSE)
write.csv(model_top, file.path(out_dir, "msdial_like_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(single_chrom_top, file.path(out_dir, "msdial_like_single_peak_chrom_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(msdial_direct_top, file.path(out_dir, "ms2decr_msdial_deconv_direct_top_fragments.csv"), row.names = FALSE)
write.csv(msdial_top, file.path(out_dir, "ms2decr_msdial_model_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(msdial_070_direct_top, file.path(out_dir, "ms2decr_msdial_ideal070_deconv_direct_top_fragments.csv"), row.names = FALSE)
write.csv(msdial_070_top, file.path(out_dir, "ms2decr_msdial_ideal070_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(rt070_im060_top, file.path(out_dir, "ms2decr_rt070_im060_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(rt070_im060_bl_top, file.path(out_dir, "ms2decr_rt070_im060_window_baseline_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(als_top, file.path(out_dir, "mcr_als_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(als_model_init_top, file.path(out_dir, "mcr_als_model_init_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(als_over_top, file.path(out_dir, "mcr_als_overcomplete_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(als_over_selected_top, file.path(out_dir, "mcr_als_overcomplete_peak_selected_commonA_top_fragments.csv"), row.names = FALSE)
write.csv(als_over_selected_Aonly_top, file.path(out_dir, "mcr_als_overcomplete_peak_selected_Aonly_top_fragments.csv"), row.names = FALSE)
write.csv(als_over_selected_mloose_top, file.path(out_dir, "mcr_als_overcomplete_mobility_loose_peak_selected_top_fragments.csv"), row.names = FALSE)
write.csv(als_over_selected_mloose_Aonly_top, file.path(out_dir, "mcr_als_overcomplete_mobility_loose_Aonly_top_fragments.csv"), row.names = FALSE)
write.csv(als_over_selected_mshape_top, file.path(out_dir, "mcr_als_overcomplete_mobility_shape_peak_selected_top_fragments.csv"), row.names = FALSE)
write.csv(als_over_selected_mshape_Aonly_top, file.path(out_dir, "mcr_als_overcomplete_mobility_shape_Aonly_top_fragments.csv"), row.names = FALSE)
write.csv(als_over_selected_reals_top, file.path(out_dir, "mcr_als_overcomplete_peak_selected_reals_top_fragments.csv"), row.names = FALSE)

plot_profiles(chrom$axis, C_model, "MS-DIAL-like model chromatograms", "RT (min)", file.path(out_dir, "msdial_like_model_chromatograms.png"))
plot_profiles(mobil$axis, M_model, "MS-DIAL-like model mobilograms", "1/K0", file.path(out_dir, "msdial_like_model_mobilograms.png"))
plot_profiles(chrom$axis, C_single, "MS-DIAL-like single-peak model chromatograms", "RT (min)", file.path(out_dir, "msdial_like_single_peak_model_chromatograms.png"))
plot_profiles(chrom$axis, C_msdial, "MS2DecR msdial_model chromatograms", "RT (min)", file.path(out_dir, "ms2decr_msdial_model_chromatograms.png"))
plot_profiles(mobil$axis, M_msdial, "MS2DecR msdial_model-linked mobilograms", "1/K0", file.path(out_dir, "ms2decr_msdial_model_mobilograms.png"))
plot_profiles(chrom$axis, msdial_case_070$C, "MS2DecR msdial_model chromatograms ideal_min 0.70", "RT (min)", file.path(out_dir, "ms2decr_msdial_ideal070_chromatograms.png"))
plot_profiles(mobil$axis, msdial_case_070$M, "MS2DecR msdial_model-linked mobilograms ideal_min 0.70", "1/K0", file.path(out_dir, "ms2decr_msdial_ideal070_mobilograms.png"))
plot_profiles(mobil$axis, M_axis, "MS2DecR mobilogram models ideal_min 0.60 fwhm 3", "1/K0", file.path(out_dir, "ms2decr_mobilogram_ideal060_fwhm3_models.png"))
plot_profiles(mobil$axis, M_im060, "Aligned MS2DecR mobilogram models ideal_min 0.60 fwhm 3", "1/K0", file.path(out_dir, "ms2decr_rt070_im060_aligned_mobilograms.png"))
plot_profiles(chrom$axis, als$C, "Joint MCR-ALS chromatographic profiles", "RT (min)", file.path(out_dir, "mcr_als_chromatograms.png"))
plot_profiles(mobil$axis, als$M, "Joint MCR-ALS mobility profiles", "1/K0", file.path(out_dir, "mcr_als_mobilograms.png"))
plot_profiles(chrom$axis, als_model_init$C, "Model-initialized MCR-ALS chromatographic profiles", "RT (min)", file.path(out_dir, "mcr_als_model_init_chromatograms.png"))
plot_profiles(mobil$axis, als_model_init$M, "Model-initialized MCR-ALS mobility profiles", "1/K0", file.path(out_dir, "mcr_als_model_init_mobilograms.png"))
plot_profiles(chrom$axis, als_over$C, "Overcomplete joint MCR-ALS chromatographic profiles", "RT (min)", file.path(out_dir, "mcr_als_overcomplete_chromatograms.png"))
plot_profiles(mobil$axis, als_over$M, "Overcomplete joint MCR-ALS mobility profiles", "1/K0", file.path(out_dir, "mcr_als_overcomplete_mobilograms.png"))
plot_profiles(chrom$axis, als_over_selected$C, "Peak-selected overcomplete MCR-ALS chromatographic profiles", "RT (min)", file.path(out_dir, "mcr_als_overcomplete_peak_selected_chromatograms.png"))
plot_profiles(mobil$axis, als_over_selected$M, "Peak-selected overcomplete MCR-ALS mobility profiles", "1/K0", file.path(out_dir, "mcr_als_overcomplete_peak_selected_mobilograms.png"))
plot_profiles(chrom$axis, als_over_selected_mshape$C, "Mobility-shape selected overcomplete MCR-ALS chromatographic profiles", "RT (min)", file.path(out_dir, "mcr_als_overcomplete_mobility_shape_chromatograms.png"))
plot_profiles(mobil$axis, als_over_selected_mshape$M, "Mobility-shape selected overcomplete MCR-ALS mobility profiles", "1/K0", file.path(out_dir, "mcr_als_overcomplete_mobility_shape_mobilograms.png"))
plot_deconv_summary(chrom$axis, mobil$axis, C_model, M_model, A_model, labels,
                    "MS-DIAL-like common-A deconvolution",
                    file.path(out_dir, "msdial_like_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, C_single, M_model, A_single_chrom, labels,
                    "MS-DIAL-like single-peak chromatogram common-A deconvolution",
                    file.path(out_dir, "msdial_like_single_peak_chrom_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, C_msdial, M_msdial, A_msdial, labels,
                    "MS2DecR msdial_model RT-selected common-A deconvolution",
                    file.path(out_dir, "ms2decr_msdial_model_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, msdial_case_070$C, msdial_case_070$M, msdial_case_070$A, labels,
                    "MS2DecR msdial_model ideal_min 0.70 common-A deconvolution",
                    file.path(out_dir, "ms2decr_msdial_ideal070_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, C_rt070, M_im060, A_rt070_im060, labels,
                    "MS2DecR RT ideal 0.70 + IM ideal 0.60 fwhm 3 common-A deconvolution",
                    file.path(out_dir, "ms2decr_rt070_im060_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, C_rt070, M_im060, A_rt070_im060_bl, labels,
                    "MS2DecR RT+IM window baseline common-A deconvolution",
                    file.path(out_dir, "ms2decr_rt070_im060_window_baseline_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als$C, als$M, als$A, labels,
                    "Joint MCR-ALS common-A deconvolution",
                    file.path(out_dir, "mcr_als_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als_model_init$C, als_model_init$M, als_model_init$A, labels,
                    "Model-initialized MCR-ALS common-A deconvolution",
                    file.path(out_dir, "mcr_als_model_init_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als_over$C, als_over$M, als_over$A, labels,
                    "Overcomplete MCR-ALS all 5 components before peak selection",
                    file.path(out_dir, "mcr_als_overcomplete_all_components_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als_over_selected$C, als_over_selected$M, als_over_selected$A, labels,
                    "Overcomplete peak-selected MCR-ALS common-A deconvolution",
                    file.path(out_dir, "mcr_als_overcomplete_peak_selected_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als_over_selected_Aonly$C, als_over_selected_Aonly$M, als_over_selected_Aonly$A, labels,
                    "Overcomplete MCR-ALS selected C/M fixed, A-only refit",
                    file.path(out_dir, "mcr_als_overcomplete_peak_selected_Aonly_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als_over_selected_mloose$C, als_over_selected_mloose$M, als_over_selected_mloose$A, labels,
                    "Overcomplete MCR-ALS mobility-loose peak selection",
                    file.path(out_dir, "mcr_als_overcomplete_mobility_loose_peak_selected_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als_over_selected_mloose_Aonly$C, als_over_selected_mloose_Aonly$M, als_over_selected_mloose_Aonly$A, labels,
                    "Overcomplete MCR-ALS mobility-loose selected C/M fixed, A-only refit",
                    file.path(out_dir, "mcr_als_overcomplete_mobility_loose_Aonly_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als_over_selected_mshape$C, als_over_selected_mshape$M, als_over_selected_mshape$A, labels,
                    "Overcomplete MCR-ALS mobility-shape peak selection",
                    file.path(out_dir, "mcr_als_overcomplete_mobility_shape_peak_selected_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als_over_selected_mshape_Aonly$C, als_over_selected_mshape_Aonly$M, als_over_selected_mshape_Aonly$A, labels,
                    "Overcomplete MCR-ALS mobility-shape selected C/M fixed, A-only refit",
                    file.path(out_dir, "mcr_als_overcomplete_mobility_shape_Aonly_deconvolution_summary.png"))
plot_deconv_summary(chrom$axis, mobil$axis, als_over_selected_reals$C, als_over_selected_reals$M, als_over_selected_reals$A, labels,
                    "Overcomplete MCR-ALS re-ALS after peak selection",
                    file.path(out_dir, "mcr_als_overcomplete_peak_selected_reals_deconvolution_summary.png"))

summary_file <- file.path(out_dir, "summary.md")
cat(
  "# diaPASEF Fig.6 common-spectrum deconvolution test\n\n",
  "Input matrices:\n\n",
  sprintf("- Chromatogram matrix: `%d x %d`\n", nrow(X), ncol(X)),
  sprintf("- Mobilogram matrix: `%d x %d`\n", nrow(Y), ncol(Y)),
  "\nThe same fragment columns are used in both matrices, so the shared spectrum matrix `A` can be estimated from\n",
  "`argmin ||X - C A||_F^2 + ||Y - M A||_F^2`.\n\n",
  "## MS-DIAL-like fixed-model result\n\n",
  paste(capture.output(print(model_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(model_top, row.names = FALSE)), collapse = "\n"),
  "\n\n## MS-DIAL-like single-peak chromatogram model result\n\n",
  "Selected chromatogram model fragments:\n\n",
  paste(capture.output(print(single_chrom$chosen, row.names = FALSE)), collapse = "\n"),
  "\n\n",
  paste(capture.output(print(single_chrom_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(single_chrom_top, row.names = FALSE)), collapse = "\n"),
  "\n\n## MS2DecR msdial_model RT-selected model result\n\n",
  "Focused peak:\n\n",
  paste(capture.output(print(data.frame(
    left_rt = chrom$axis[focused_peak$left_idx],
    apex_rt = chrom$axis[focused_peak$apex_idx],
    right_rt = chrom$axis[focused_peak$right_idx]
  ), row.names = FALSE)), collapse = "\n"),
  "\n\nTriplet summary:\n\n",
  paste(capture.output(print(msdial_obj$result$model$triplet_summary, row.names = FALSE)), collapse = "\n"),
  "\n\nDirect MS2DecR `msdial_deconv()` purity:\n\n",
  paste(capture.output(print(msdial_direct_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nDirect MS2DecR top fragments:\n\n",
  paste(capture.output(print(msdial_direct_top, row.names = FALSE)), collapse = "\n"),
  "\n\nRT+IM common-A purity using the same RT models:\n\n",
  paste(capture.output(print(msdial_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(msdial_top, row.names = FALSE)), collapse = "\n"),
  "\n\n## MS2DecR msdial_model ideal_min 0.70 result\n\n",
  "Triplet summary:\n\n",
  paste(capture.output(print(msdial_case_070$obj$result$model$triplet_summary, row.names = FALSE)), collapse = "\n"),
  "\n\nDirect MS2DecR `msdial_deconv()` purity:\n\n",
  paste(capture.output(print(msdial_070_direct_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nDirect MS2DecR top fragments:\n\n",
  paste(capture.output(print(msdial_070_direct_top, row.names = FALSE)), collapse = "\n"),
  "\n\nRT+IM common-A purity using the same RT models:\n\n",
  paste(capture.output(print(msdial_070_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(msdial_070_top, row.names = FALSE)), collapse = "\n"),
  "\n\n## MS2DecR RT ideal_min 0.70 + IM ideal_min 0.60 fwhm 3 result\n\n",
  "RT triplet summary:\n\n",
  paste(capture.output(print(msdial_case_070$obj$result$model$triplet_summary, row.names = FALSE)), collapse = "\n"),
  "\n\nIM triplet summary:\n\n",
  paste(capture.output(print(mobility_model_obj$result$model$triplet_summary, row.names = FALSE)), collapse = "\n"),
  "\n\nModel site alignment:\n\n",
  paste(capture.output(print(data.frame(
    model = names(rt_model_sites),
    rt_site = unname(rt_model_sites),
    im_site = unname(im_model_sites[match(rt_model_sites, im_model_sites)])
  ), row.names = FALSE)), collapse = "\n"),
  "\n\nRT+IM common-A purity:\n\n",
  paste(capture.output(print(rt070_im060_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(rt070_im060_top, row.names = FALSE)), collapse = "\n"),
  "\n\nRT+IM window + baseline common-A purity:\n\n",
  paste(capture.output(print(rt070_im060_bl_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(rt070_im060_bl_top, row.names = FALSE)), collapse = "\n"),
  "\n\n## Joint MCR-ALS result\n\n",
  paste(capture.output(print(als_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(als_top, row.names = FALSE)), collapse = "\n"),
  "\n\n## ALS error\n\n",
  sprintf("- initial: %.4f\n- final: %.4f\n", als$error[1], tail(als$error, 1)),
  "\n## Model-initialized joint MCR-ALS result\n\n",
  paste(capture.output(print(als_model_init_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(als_model_init_top, row.names = FALSE)), collapse = "\n"),
  "\n\n## Model-initialized ALS error\n\n",
  sprintf("- initial: %.4f\n- final: %.4f\n", als_model_init$error[1], tail(als_model_init$error, 1)),
  "\n## Overcomplete peak-selected MCR-ALS result\n\n",
  "Component selection table:\n\n",
  paste(capture.output(print(als_over_selected$component_table, row.names = FALSE)), collapse = "\n"),
  "\n\nSelected component purity:\n\n",
  paste(capture.output(print(als_over_selected_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(als_over_selected_top, row.names = FALSE)), collapse = "\n"),
  "\n\n## MCR-ALS post-selection strategy comparison\n\n",
  paste(capture.output(print(als_strategy_comparison, row.names = FALSE)), collapse = "\n"),
  "\n\n### B. Fixed C/M, A-only window-baseline refit\n\n",
  paste(capture.output(print(als_over_selected_Aonly_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(als_over_selected_Aonly_top, row.names = FALSE)), collapse = "\n"),
  "\n\n### B2. Mobility-loose peak selection\n\n",
  "Component selection table:\n\n",
  paste(capture.output(print(als_over_selected_mloose$component_table, row.names = FALSE)), collapse = "\n"),
  "\n\n",
  paste(capture.output(print(als_over_selected_mloose_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(als_over_selected_mloose_top, row.names = FALSE)), collapse = "\n"),
  "\n\n### B3. Mobility-loose fixed C/M, A-only window-baseline refit\n\n",
  paste(capture.output(print(als_over_selected_mloose_Aonly_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(als_over_selected_mloose_Aonly_top, row.names = FALSE)), collapse = "\n"),
  "\n\n### B4. Mobility-shape peak selection\n\n",
  "Component selection table:\n\n",
  paste(capture.output(print(als_over_selected_mshape$component_table, row.names = FALSE)), collapse = "\n"),
  "\n\n",
  paste(capture.output(print(als_over_selected_mshape_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(als_over_selected_mshape_top, row.names = FALSE)), collapse = "\n"),
  "\n\n### B5. Mobility-shape fixed C/M, A-only window-baseline refit\n\n",
  paste(capture.output(print(als_over_selected_mshape_Aonly_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(als_over_selected_mshape_Aonly_top, row.names = FALSE)), collapse = "\n"),
  "\n\n### C. Re-ALS after peak selection\n\n",
  paste(capture.output(print(als_over_selected_reals_purity, row.names = FALSE)), collapse = "\n"),
  "\n\nTop fragments:\n\n",
  paste(capture.output(print(als_over_selected_reals_top, row.names = FALSE)), collapse = "\n"),
  "\n\n## Overcomplete ALS error\n\n",
  sprintf("- initial: %.4f\n- final: %.4f\n", als_over$error[1], tail(als_over$error, 1)),
  "\n## Summary Images\n\n",
  "- `msdial_like_deconvolution_summary.png`\n",
  "- `msdial_like_single_peak_chrom_deconvolution_summary.png`\n",
  "- `ms2decr_msdial_model_deconvolution_summary.png`\n",
  "- `ms2decr_msdial_ideal070_deconvolution_summary.png`\n",
  "- `ms2decr_rt070_im060_deconvolution_summary.png`\n",
  "- `ms2decr_rt070_im060_window_baseline_deconvolution_summary.png`\n",
  "- `mcr_als_deconvolution_summary.png`\n",
  "- `mcr_als_model_init_deconvolution_summary.png`\n",
  "- `mcr_als_overcomplete_all_components_deconvolution_summary.png`\n",
  "- `mcr_als_overcomplete_peak_selected_deconvolution_summary.png`\n",
  "- `mcr_als_overcomplete_peak_selected_Aonly_deconvolution_summary.png`\n",
  "- `mcr_als_overcomplete_mobility_loose_peak_selected_deconvolution_summary.png`\n",
  "- `mcr_als_overcomplete_mobility_loose_Aonly_deconvolution_summary.png`\n",
  "- `mcr_als_overcomplete_mobility_shape_peak_selected_deconvolution_summary.png`\n",
  "- `mcr_als_overcomplete_mobility_shape_Aonly_deconvolution_summary.png`\n",
  "- `mcr_als_overcomplete_peak_selected_reals_deconvolution_summary.png`\n",
  file = summary_file,
  sep = ""
)

cat("Wrote outputs to:", normalizePath(out_dir, winslash = "/"), "\n")
cat("\nMS-DIAL-like fixed-model purity:\n")
print(model_purity, row.names = FALSE)
cat("\nMS-DIAL-like single-peak chromatogram model fragments:\n")
print(single_chrom$chosen, row.names = FALSE)
cat("\nMS-DIAL-like single-peak chromatogram purity:\n")
print(single_chrom_purity, row.names = FALSE)
cat("\nMS2DecR msdial_model focused peak:\n")
print(data.frame(
  left_rt = chrom$axis[focused_peak$left_idx],
  apex_rt = chrom$axis[focused_peak$apex_idx],
  right_rt = chrom$axis[focused_peak$right_idx]
), row.names = FALSE)
cat("\nMS2DecR msdial_model triplet summary:\n")
print(msdial_obj$result$model$triplet_summary, row.names = FALSE)
cat("\nDirect MS2DecR msdial_deconv purity:\n")
print(msdial_direct_purity, row.names = FALSE)
cat("\nMS2DecR msdial_model RT-selected purity:\n")
print(msdial_purity, row.names = FALSE)
cat("\nMS2DecR msdial_model ideal_min 0.70 triplet summary:\n")
print(msdial_case_070$obj$result$model$triplet_summary, row.names = FALSE)
cat("\nDirect MS2DecR msdial_deconv ideal_min 0.70 purity:\n")
print(msdial_070_direct_purity, row.names = FALSE)
cat("\nMS2DecR msdial_model ideal_min 0.70 common-A purity:\n")
print(msdial_070_purity, row.names = FALSE)
cat("\nMS2DecR mobilogram ideal_min 0.60 fwhm 3 triplet summary:\n")
print(mobility_model_obj$result$model$triplet_summary, row.names = FALSE)
cat("\nMS2DecR RT ideal_min 0.70 + IM ideal_min 0.60 fwhm 3 common-A purity:\n")
print(rt070_im060_purity, row.names = FALSE)
cat("\nMS2DecR RT+IM window baseline common-A purity:\n")
print(rt070_im060_bl_purity, row.names = FALSE)
cat("\nJoint MCR-ALS purity:\n")
print(als_purity, row.names = FALSE)
cat(sprintf("\nJoint MCR-ALS error: %.4f -> %.4f\n", als$error[1], tail(als$error, 1)))
cat("\nModel-initialized joint MCR-ALS purity:\n")
print(als_model_init_purity, row.names = FALSE)
cat(sprintf("\nModel-initialized joint MCR-ALS error: %.4f -> %.4f\n", als_model_init$error[1], tail(als_model_init$error, 1)))
cat("\nOvercomplete peak-selected MCR-ALS component table:\n")
print(als_over_selected$component_table, row.names = FALSE)
cat("\nOvercomplete peak-selected MCR-ALS purity:\n")
print(als_over_selected_purity, row.names = FALSE)
cat(sprintf("\nOvercomplete MCR-ALS error: %.4f -> %.4f\n", als_over$error[1], tail(als_over$error, 1)))
cat("\nMCR-ALS post-selection strategy comparison:\n")
print(als_strategy_comparison, row.names = FALSE)
cat("\nFixed C/M A-only refit purity:\n")
print(als_over_selected_Aonly_purity, row.names = FALSE)
cat("\nMobility-loose peak-selected MCR-ALS component table:\n")
print(als_over_selected_mloose$component_table, row.names = FALSE)
cat("\nMobility-loose peak-selected MCR-ALS purity:\n")
print(als_over_selected_mloose_purity, row.names = FALSE)
cat("\nMobility-loose fixed C/M A-only refit purity:\n")
print(als_over_selected_mloose_Aonly_purity, row.names = FALSE)
cat("\nMobility-shape peak-selected MCR-ALS component table:\n")
print(als_over_selected_mshape$component_table, row.names = FALSE)
cat("\nMobility-shape peak-selected MCR-ALS purity:\n")
print(als_over_selected_mshape_purity, row.names = FALSE)
cat("\nMobility-shape fixed C/M A-only refit purity:\n")
print(als_over_selected_mshape_Aonly_purity, row.names = FALSE)
cat("\nRe-ALS after peak selection purity:\n")
print(als_over_selected_reals_purity, row.names = FALSE)
