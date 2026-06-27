# Simulate MS1-guided DDA-MS2 spectral unmixing for leucine/isoleucine-like isomers.
#
# This is a toy simulation, not a real Leu/Ile spectral library example.
# It demonstrates:
#   1. same precursor m/z can be represented by two RT-resolved MS1 Gaussian peaks;
#   2. the Gaussian peak heights at DDA-MS2 scan times form the mixing matrix A;
#   3. chimeric DDA-MS2 spectra can be unmixed by A-fixed non-negative least squares;
#   4. if MS1 is treated as one merged same-m/z EIC, only a mixed MS2 spectrum is obtained.

set.seed(7)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
  script_dir <- dirname(script_path)
} else {
  script_dir <- getwd()
}

out_dir <- file.path(script_dir, "output")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

gaussian <- function(x, apex, sigma, height) {
  height * exp(-0.5 * ((x - apex) / sigma)^2)
}

normalize_rows <- function(x) {
  out <- x
  for (i in seq_len(nrow(out))) {
    mx <- max(out[i, ])
    if (mx > 0) {
      out[i, ] <- out[i, ] / mx
    }
  }
  out
}

# A compact NNLS solver for this simulation.
# It solves min ||A x - y||^2 subject to x >= 0.
# For general production use, use a standard NNLS implementation.
nnls_small <- function(A, y) {
  p <- ncol(A)
  best_x <- rep(0, p)
  best_rss <- sum(y^2)

  for (mask in seq_len(2^p - 1)) {
    active <- as.logical(intToBits(mask)[seq_len(p)])
    Aa <- A[, active, drop = FALSE]
    x_active <- as.vector(qr.solve(Aa, y))

    if (all(is.finite(x_active)) && all(x_active >= 0)) {
      x <- rep(0, p)
      x[active] <- x_active
      rss <- sum((y - A %*% x)^2)
      if (rss < best_rss) {
        best_rss <- rss
        best_x <- x
      }
    }
  }

  list(x = best_x, rss = best_rss)
}

write_plot <- function(filename, expr, width = 900, height = 560) {
  png(file.path(out_dir, filename), width = width, height = height, res = 120)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

# MS1 EIC: two isomer-like components with the same precursor m/z but shifted RT.
rt <- seq(4.6, 5.8, length.out = 241)
leu_eic <- gaussian(rt, apex = 5.08, sigma = 0.13, height = 120000)
ile_eic <- gaussian(rt, apex = 5.32, sigma = 0.14, height = 100000)
merged_eic <- leu_eic + ile_eic

# Sparse DDA-MS2 scan times.
ms2_rt <- seq(4.85, 5.58, length.out = 44) + rnorm(44, mean = 0, sd = 0.006)
ms2_rt <- sort(ms2_rt)

leu_at_ms2 <- approx(rt, leu_eic, xout = ms2_rt, rule = 2)$y
ile_at_ms2 <- approx(rt, ile_eic, xout = ms2_rt, rule = 2)$y

# MS1-derived mixing matrix.
A <- cbind(Leu = leu_at_ms2, Ile = ile_at_ms2)

# Toy fragment bins and toy spectra.
fragment_mz <- c(30, 41, 44, 69, 86, 132)
S_true <- rbind(
  Leu = c(0.12, 0.20, 0.75, 0.08, 1.00, 0.18),
  Ile = c(0.10, 0.62, 0.22, 0.80, 1.00, 0.15)
)
colnames(S_true) <- paste0("mz_", fragment_mz)

# Chimeric DDA-MS2 spectra.
Y_clean <- A %*% S_true
noise <- matrix(
  rnorm(length(Y_clean), mean = 0, sd = 0.025 * max(Y_clean)),
  nrow = nrow(Y_clean),
  ncol = ncol(Y_clean)
)
Y <- pmax(Y_clean + noise, 0)

# Estimate Leu/Ile MS2 spectra by NNLS, one fragment bin at a time.
S_hat <- matrix(0, nrow = 2, ncol = ncol(Y))
rownames(S_hat) <- c("Leu", "Ile")
colnames(S_hat) <- colnames(S_true)
rss <- numeric(ncol(Y))

for (j in seq_len(ncol(Y))) {
  fit <- nnls_small(A, Y[, j])
  S_hat[, j] <- fit$x
  rss[j] <- fit$rss
}

S_true_norm <- normalize_rows(S_true)
S_hat_norm <- normalize_rows(S_hat)

# Control analysis: one merged same-m/z MS1 EIC gives only one mixed spectrum.
A_merged <- matrix(approx(rt, merged_eic, xout = ms2_rt, rule = 2)$y, ncol = 1)
S_merged <- numeric(ncol(Y))
for (j in seq_len(ncol(Y))) {
  fit <- nnls_small(A_merged, Y[, j])
  S_merged[j] <- fit$x[1]
}
S_merged_norm <- S_merged / max(S_merged)

# Diagnostics.
corr_leu <- cor(S_true_norm["Leu", ], S_hat_norm["Leu", ])
corr_ile <- cor(S_true_norm["Ile", ], S_hat_norm["Ile", ])
condition_number <- kappa(A)

# Write outputs.
write.csv(
  data.frame(
    rt_min = rt,
    leucine_ms1_eic = leu_eic,
    isoleucine_ms1_eic = ile_eic,
    merged_same_mz_eic = merged_eic
  ),
  file.path(out_dir, "ms1_eic.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    ms2_scan = seq_along(ms2_rt),
    rt_min = ms2_rt,
    A_leucine = A[, "Leu"],
    A_isoleucine = A[, "Ile"]
  ),
  file.path(out_dir, "design_matrix_A.csv"),
  row.names = FALSE
)

estimated_spectra <- data.frame(
  fragment_mz = fragment_mz,
  true_leucine = as.numeric(S_true_norm["Leu", ]),
  estimated_leucine = as.numeric(S_hat_norm["Leu", ]),
  true_isoleucine = as.numeric(S_true_norm["Ile", ]),
  estimated_isoleucine = as.numeric(S_hat_norm["Ile", ]),
  wrong_merged_estimate = as.numeric(S_merged_norm)
)
write.csv(
  estimated_spectra,
  file.path(out_dir, "estimated_spectra.csv"),
  row.names = FALSE
)

observed_ms2 <- as.data.frame(Y)
colnames(observed_ms2) <- paste0("frag_", fragment_mz)
observed_ms2$rt_min <- ms2_rt
write.csv(
  observed_ms2,
  file.path(out_dir, "observed_chimeric_ms2.csv"),
  row.names = FALSE
)

summary_lines <- c(
  paste("Output directory:", out_dir),
  paste("Correlation true vs estimated Leu:", round(corr_leu, 4)),
  paste("Correlation true vs estimated Ile:", round(corr_ile, 4)),
  paste("Condition number of A:", round(condition_number, 3)),
  "",
  "Estimated spectra:",
  paste(capture.output(print(estimated_spectra, row.names = FALSE)), collapse = "\n")
)
writeLines(summary_lines, file.path(out_dir, "summary.txt"))

# Plot 1: MS1 EIC and DDA-MS2 scan times.
write_plot("01_ms1_eic_and_dda_times.png", {
  plot(rt, leu_eic, type = "l", lwd = 2, col = "#1f77b4",
       xlab = "Retention time (min)", ylab = "MS1 intensity",
       main = "Same precursor m/z, separated only by RT",
       ylim = c(0, max(merged_eic) * 1.08))
  lines(rt, ile_eic, lwd = 2, col = "#d62728")
  lines(rt, merged_eic, lwd = 2, col = "gray25", lty = 2)
  points(ms2_rt, approx(rt, merged_eic, xout = ms2_rt, rule = 2)$y, pch = 16, cex = 0.7)
  legend("topright",
         legend = c("Leucine MS1 component", "Isoleucine MS1 component",
                    "Merged same-m/z EIC", "DDA-MS2 scan times"),
         col = c("#1f77b4", "#d62728", "gray25", "black"),
         lty = c(1, 1, 2, NA), pch = c(NA, NA, NA, 16), lwd = c(2, 2, 2, NA),
         bty = "n")
})

# Plot 2: design matrix A.
write_plot("02_design_matrix_A.png", {
  image(
    x = seq_len(ncol(A)),
    y = seq_len(nrow(A)),
    z = t(A / max(A)),
    axes = FALSE,
    xlab = "MS1 component",
    ylab = "DDA-MS2 scan",
    main = "MS1-derived design matrix A",
    col = hcl.colors(64, "Viridis")
  )
  axis(1, at = c(1, 2), labels = c("Leu", "Ile"))
  tick_idx <- seq(1, length(ms2_rt), by = 5)
  axis(2, at = tick_idx, labels = sprintf("%.2f", ms2_rt[tick_idx]), las = 1)
  box()
})

# Plot 3: true vs estimated spectra.
write_plot("03_true_vs_estimated_spectra.png", expr = {
  old_par <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  on.exit(par(old_par), add = TRUE)
  for (compound in c("Leu", "Ile")) {
    barplot(
      rbind(S_true_norm[compound, ], S_hat_norm[compound, ]),
      beside = TRUE,
      col = c("gray35", "#2ca02c"),
      names.arg = fragment_mz,
      ylim = c(0, 1.12),
      ylab = "Relative intensity",
      xlab = "Fragment m/z bin",
      main = compound
    )
    legend("topright", legend = c("True", "Estimated from chimeric DDA-MS2"),
           fill = c("gray35", "#2ca02c"), bty = "n")
  }
}, height = 720)

# Plot 4: merged-feature control.
write_plot("04_wrong_merged_feature.png", {
  barplot(
    rbind(S_true_norm["Leu", ], S_true_norm["Ile", ], S_merged_norm),
    beside = TRUE,
    col = c("#1f77b4", "#d62728", "gray35"),
    names.arg = fragment_mz,
    ylim = c(0, 1.12),
    ylab = "Relative intensity",
    xlab = "Fragment m/z bin",
    main = "Merged MS1 feature gives only a mixed spectrum"
  )
  legend("topright",
         legend = c("True Leu", "True Ile", "Estimated from merged same-m/z EIC"),
         fill = c("#1f77b4", "#d62728", "gray35"), bty = "n")
})

cat(paste(summary_lines, collapse = "\n"))
cat("\n")
