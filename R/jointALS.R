.jointALS_normalize_trace <- function(x) {
  x <- pmax(as.numeric(x), 0)
  s <- sqrt(sum(x^2))
  if (!is.finite(s) || s == 0) return(x)
  x / s
}

.jointALS_cosine <- function(x, y) {
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(den) || den == 0) return(NA_real_)
  sum(x * y) / den
}

.jointALS_smooth_trace <- function(x, width = 5) {
  if (length(x) < width) return(as.numeric(x))
  y <- stats::filter(as.numeric(x), rep(1 / width, width), sides = 2)
  y[is.na(y)] <- x[is.na(y)]
  as.numeric(y)
}

.jointALS_local_peak_count <- function(x, frac = 0.25) {
  if (length(x) < 3L) return(0L)
  threshold <- max(x, na.rm = TRUE) * frac
  sum(x[-c(1L, length(x))] > x[-c(length(x) - 1L, length(x))] &
        x[-c(1L, 2L)] < x[-c(1L, length(x))] &
        x[-c(1L, length(x))] > threshold)
}

.jointALS_profile_peak_score <- function(x,
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
  xs <- pmax(.jointALS_smooth_trace(x, width = smooth_width), 0)
  filtered <- matchedfilter(xs, wid)$filtered_signal
  filtered[!is.finite(filtered)] <- 0
  apex <- which.max(xs)
  filtered_apex <- which.max(filtered)
  close <- abs(filtered_apex - apex) <= gap
  n_peaks <- .jointALS_local_peak_count(xs, frac = local_peak_frac)
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

.jointALS_ica_profiles <- function(X, com) {
  if (!requireNamespace("fastICA", quietly = TRUE)) {
    sv <- svd(scale(X, center = TRUE, scale = FALSE), nu = com, nv = com)
    C <- abs(sv$u[, seq_len(com), drop = FALSE])
    for (i in seq_len(com)) C[, i] <- .jointALS_normalize_trace(C[, i])
    return(C)
  }

  ic <- fastICA::fastICA(t(X), n.comp = com)
  C <- matrix(0, nrow(X), com)
  keep <- apply(X, 2, sd) > 0
  for (i in seq_len(com)) {
    s <- t(ic$A)[, i]
    r <- suppressWarnings(cor(X[, keep, drop = FALSE], s))
    r <- r[is.finite(r)]
    if (length(r) && r[which.max(abs(r))] < 0) s <- -s
    C[, i] <- .jointALS_normalize_trace(pmax(s, 0))
  }
  C
}

.jointALS_pair_components_by_spectra <- function(C, M, X, Y, lambda = 1e-6) {
  com <- ncol(C)
  Ax <- pmax(solve(crossprod(C) + lambda * diag(com), crossprod(C, X)), 0)
  Ay <- pmax(solve(crossprod(M) + lambda * diag(com), crossprod(M, Y)), 0)
  score <- outer(seq_len(com), seq_len(com), Vectorize(function(i, j) {
    .jointALS_cosine(Ax[i, ], Ay[j, ])
  }))
  score[!is.finite(score)] <- -Inf

  order_m <- integer(com)
  used <- logical(com)
  for (i in seq_len(com)) {
    available <- which(!used)
    best <- available[which.max(score[i, available])]
    order_m[i] <- best
    used[best] <- TRUE
  }
  order_m
}

# Detect valid joint ALS components from RT and IM profiles
detectJointALSPeaks <- function(C, M, A, c_wid = 5, m_wid = 5,
                                c_gap = 5, m_gap = 5,
                                mode = c("both", "either", "rt", "im")) {
  C <- as.matrix(C)
  M <- as.matrix(M)
  A <- as.matrix(A)
  mode <- match.arg(mode)

  if (ncol(C) != ncol(M)) stop("C and M must have the same number of components.")
  if (ncol(C) != nrow(A)) stop("ncol(C) and ncol(M) must match nrow(A).")

  component_names <- colnames(C)
  if (is.null(component_names)) component_names <- paste0("ALS", seq_len(ncol(C)))

  peak_pass <- function(x, wid, gap) {
    x <- as.numeric(x)
    if (!length(x) || !any(is.finite(x))) return(list(pass = FALSE, apex = NA_integer_, filtered_apex = NA_integer_))
    x[!is.finite(x)] <- 0
    filtered_signal <- matchedfilter(x, wid)[[1]]
    filtered_signal[!is.finite(filtered_signal)] <- 0
    apex <- which.max(x)
    filtered_apex <- which.max(filtered_signal)
    list(
      pass = abs(filtered_apex - apex) <= gap,
      apex = apex,
      filtered_apex = filtered_apex
    )
  }

  component_table <- do.call(rbind, lapply(seq_len(ncol(C)), function(i) {
    c_peak <- peak_pass(C[, i], c_wid, c_gap)
    m_peak <- peak_pass(M[, i], m_wid, m_gap)
    spectrum_nonzero <- sum(A[i, ] != 0, na.rm = TRUE) >= 1L
    keep <- switch(
      mode,
      both = c_peak$pass && m_peak$pass && spectrum_nonzero,
      either = (c_peak$pass || m_peak$pass) && spectrum_nonzero,
      rt = c_peak$pass && spectrum_nonzero,
      im = m_peak$pass && spectrum_nonzero
    )
    data.frame(
      component = component_names[i],
      index = i,
      c_pass = c_peak$pass,
      m_pass = m_peak$pass,
      c_apex = c_peak$apex,
      c_filtered_apex = c_peak$filtered_apex,
      m_apex = m_peak$apex,
      m_filtered_apex = m_peak$filtered_apex,
      spectrum_nonzero = spectrum_nonzero,
      keep = keep,
      stringsAsFactors = FALSE
    )
  }))

  keep_idx <- component_table$index[component_table$keep]
  list(
    C_fin = C[, keep_idx, drop = FALSE],
    M_fin = M[, keep_idx, drop = FALSE],
    A_fin = A[keep_idx, , drop = FALSE],
    keep_idx = keep_idx,
    component_table = component_table
  )
}

# Select components from an overcomplete joint ALS result
selectJointALSComponents <- function(fit,
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
  if (is.null(fit$C) || is.null(fit$M) || is.null(fit$A)) {
    stop("fit must contain C, M, and A.")
  }
  if (ncol(fit$C) != ncol(fit$M) || ncol(fit$C) != nrow(fit$A)) {
    stop("fit$C, fit$M, and fit$A must have the same number of components.")
  }

  com <- nrow(fit$A)
  component_names <- rownames(fit$A)
  if (is.null(component_names)) component_names <- paste0("ALS", seq_len(com))

  rows <- lapply(seq_len(com), function(i) {
    c_peak <- .jointALS_profile_peak_score(
      fit$C[, i],
      wid = c_wid,
      gap = c_gap,
      max_local_peaks = c_max_local_peaks
    )
    m_peak <- .jointALS_profile_peak_score(
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
      component = component_names[i],
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
  component_table <- do.call(rbind, rows)

  keep <- which(component_table$c_pass & component_table$m_pass & component_table$spec_sum > 0)
  if (length(keep) < min_keep) {
    keep <- order(component_table$score, decreasing = TRUE)[seq_len(min(min_keep, nrow(component_table)))]
  }
  keep <- keep[order(component_table$score[keep], decreasing = TRUE)]
  keep <- keep[seq_len(min(length(keep), max_keep))]

  selected <- list(
    C = fit$C[, keep, drop = FALSE],
    M = fit$M[, keep, drop = FALSE],
    A = fit$A[keep, , drop = FALSE],
    E = fit$E,
    error = if (!is.null(fit$error)) fit$error else fit$E,
    component_table = component_table,
    keep = keep
  )
  colnames(selected$C) <- colnames(selected$M) <- rownames(selected$A) <-
    component_table$component[keep]
  selected
}

# Refit shared MS/MS spectra from fixed joint ALS profiles
refitJointALSSpectra <- function(X, Y, C, M, lambda = 1e-8,
                                 baseline = c("window", "none")) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  C <- as.matrix(C)
  M <- as.matrix(M)
  baseline <- match.arg(baseline)

  if (nrow(C) != nrow(X)) stop("nrow(C) must match nrow(X).")
  if (nrow(M) != nrow(Y)) stop("nrow(M) must match nrow(Y).")
  if (ncol(C) != ncol(M)) stop("C and M must have the same number of components.")
  if (ncol(X) != ncol(Y)) stop("X and Y must have the same fragment columns.")

  com <- ncol(C)
  if (baseline == "none") {
    lhs <- crossprod(C) + crossprod(M) + lambda * diag(com)
    rhs <- crossprod(C, X) + crossprod(M, Y)
    spectra <- solve(lhs, rhs)
    spectra[!is.finite(spectra)] <- 0
    spectra <- pmax(spectra, 0)
    rownames(spectra) <- colnames(C)
    colnames(spectra) <- colnames(X)
    return(list(
      spectra = spectra,
      A = spectra,
      fitted_X = C %*% spectra,
      fitted_Y = M %*% spectra,
      baseline = NULL
    ))
  }

  rt_window <- which(rowSums(C != 0, na.rm = TRUE) > 0)
  im_window <- which(rowSums(M != 0, na.rm = TRUE) > 0)
  if (!length(rt_window)) rt_window <- seq_len(nrow(C))
  if (!length(im_window)) im_window <- seq_len(nrow(M))

  rt_t <- as.numeric(scale(seq_along(rt_window), center = TRUE, scale = TRUE))
  im_t <- as.numeric(scale(seq_along(im_window), center = TRUE, scale = TRUE))
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

  spectra <- matrix(0, nrow = com, ncol = ncol(X),
                    dimnames = list(colnames(C), colnames(X)))
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
    spectra[, j] <- pmax(beta[seq_len(com)], 0)
    nuisance[, j] <- beta[com + seq_len(4)]
  }

  list(
    spectra = spectra,
    A = spectra,
    nuisance = nuisance,
    rt_window = rt_window,
    im_window = im_window,
    design = design,
    fitted_X = C %*% spectra,
    fitted_Y = M %*% spectra
  )
}

# Perform joint Alternating Least Squares optimization for RT and IM matrices
optimizeJointALS <- function(X, Y, C, M, com = ncol(C), lambda = 1e-6, maxiter = 500) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  C <- as.matrix(C)
  M <- as.matrix(M)

  if (nrow(C) != nrow(X)) stop("nrow(C) must match nrow(X).")
  if (nrow(M) != nrow(Y)) stop("nrow(M) must match nrow(Y).")
  if (ncol(C) != com || ncol(M) != com) stop("C and M must have com columns.")
  if (ncol(X) != ncol(Y)) stop("X and Y must have the same fragment columns.")
  if (com < 1L) stop("com must be at least 1.")
  if (lambda < 0) stop("lambda must be non-negative.")
  if (maxiter < 1L) stop("maxiter must be at least 1.")

  Alist <- vector("list", maxiter)
  Clist <- vector("list", maxiter)
  Mlist <- vector("list", maxiter)
  E <- numeric(maxiter)

  for (k in seq_len(maxiter)) {
    A <- solve(crossprod(C) + crossprod(M) + lambda * diag(1, com)) %*%
      (crossprod(C, X) + crossprod(M, Y))
    A[!is.finite(A)] <- 0
    A[A < 0] <- 0

    for (i in seq_len(com)) {
      denom <- sqrt(sum(A[i, ]^2))
      if (is.finite(denom) && denom > 0) A[i, ] <- A[i, ] / denom
    }

    gram <- A %*% t(A) + lambda * diag(1, com)
    C <- X %*% t(A) %*% solve(gram)
    M <- Y %*% t(A) %*% solve(gram)
    C[!is.finite(C)] <- 0
    M[!is.finite(M)] <- 0
    C[C < 0] <- 0
    M[M < 0] <- 0

    E[k] <- sqrt(sum((X - C %*% A)^2) + sum((Y - M %*% A)^2))
    Alist[[k]] <- A
    Clist[[k]] <- C
    Mlist[[k]] <- M
  }

  min_index <- which.min(E)
  A <- Alist[[min_index]]
  C <- Clist[[min_index]]
  M <- Mlist[[min_index]]

  rownames(A) <- colnames(C)
  colnames(A) <- colnames(X)

  list(
    A = A,
    C = C,
    M = M,
    E = E,
    error = E,
    best_iteration = min_index,
    fitted_X = C %*% A,
    fitted_Y = M %*% A
  )
}

# Run joint ALS with optional RT and IM initial component profiles
jointALS <- function(X, Y, com = 2L, C_init = NULL, M_init = NULL,
                     lambda = 1e-6, maxiter = 500, align_m = TRUE,
                     detect_peaks = TRUE, c_wid = 5, m_wid = 5,
                     c_gap = 5, m_gap = 5,
                     peak_mode = c("both", "either", "rt", "im")) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (ncol(X) != ncol(Y)) stop("X and Y must have the same fragment columns.")
  if (com < 1L) stop("com must be at least 1.")

  if (is.null(C_init)) {
    C_init <- .jointALS_ica_profiles(X, com)
  }
  if (is.null(M_init)) {
    M_init <- .jointALS_ica_profiles(Y, com)
  }

  if (isTRUE(align_m)) {
    M_init <- M_init[, .jointALS_pair_components_by_spectra(C_init, M_init, X, Y, lambda),
                     drop = FALSE]
  }

  component_names <- paste0("ALS", seq_len(com))
  if (is.null(colnames(C_init))) colnames(C_init) <- component_names
  if (is.null(colnames(M_init))) colnames(M_init) <- component_names

  fit <- optimizeJointALS(
    X = X,
    Y = Y,
    C = C_init,
    M = M_init,
    com = com,
    lambda = lambda,
    maxiter = maxiter
  )
  colnames(fit$C) <- component_names
  colnames(fit$M) <- component_names
  rownames(fit$A) <- component_names
  fit$com <- com

  if (isTRUE(detect_peaks)) {
    fit$peaks <- detectJointALSPeaks(
      C = fit$C,
      M = fit$M,
      A = fit$A,
      c_wid = c_wid,
      m_wid = m_wid,
      c_gap = c_gap,
      m_gap = m_gap,
      mode = match.arg(peak_mode)
    )
    fit$selected <- list(
      C = fit$peaks$C_fin,
      M = fit$peaks$M_fin,
      A = fit$peaks$A_fin,
      component_table = fit$peaks$component_table
    )
  }

  fit
}
