# Simulate MS1-guided unmixing of two chimeric DDA-MS2 spectra.
#
# This is a toy leucine/isoleucine-like example.
# It assumes only two DDA-MS2 spectra:
#   1. one acquired near the Leu-like peak top;
#   2. one acquired near the Ile-like peak top.
#
# The point is to demonstrate a simple linear unmixing/inverse-matrix idea,
# not MCR-ALS and not a many-scan DDA setting.

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

# Compact NNLS solver for small examples.
# It solves min ||A x - y||^2 subject to x >= 0.
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

# MS1 EIC: two same-m/z isomer-like components separated by RT.
rt <- seq(4.6, 5.8, length.out = 241)
leu_eic <- gaussian(rt, apex = 5.08, sigma = 0.13, height = 120000)
ile_eic <- gaussian(rt, apex = 5.32, sigma = 0.14, height = 100000)
merged_eic <- leu_eic + ile_eic

# Only two DDA-MS2 spectra are assumed.
# One is acquired near the Leu-like apex, the other near the Ile-like apex.
ms2_rt <- c(5.08, 5.32)
scan_label <- c("Leu_apex_scan", "Ile_apex_scan")

leu_at_ms2 <- approx(rt, leu_eic, xout = ms2_rt, rule = 2)$y
ile_at_ms2 <- approx(rt, ile_eic, xout = ms2_rt, rule = 2)$y

# MS1-derived 2 x 2 mixing matrix.
A <- cbind(Leu = leu_at_ms2, Ile = ile_at_ms2)
rownames(A) <- scan_label

# Toy fragment bins and toy spectra.
# These are illustrative values, not real Leu/Ile library spectra.
fragment_mz <- c(30, 41, 44, 69, 86, 132)
S_true <- rbind(
  Leu = c(0.12, 0.20, 0.75, 0.08, 1.00, 0.18),
  Ile = c(0.10, 0.62, 0.22, 0.80, 1.00, 0.15)
)
colnames(S_true) <- paste0("mz_", fragment_mz)

# Two observed chimeric DDA-MS2 spectra.
Y_clean <- A %*% S_true
noise <- matrix(
  rnorm(length(Y_clean), mean = 0, sd = 0.01 * max(Y_clean)),
  nrow = nrow(Y_clean),
  ncol = ncol(Y_clean)
)
Y <- pmax(Y_clean + noise, 0)
rownames(Y) <- scan_label
colnames(Y) <- paste0("frag_", fragment_mz)

# Direct inverse solution. This is the conceptual core for the 2-scan case.
S_inverse <- solve(A, Y)
rownames(S_inverse) <- c("Leu", "Ile")
colnames(S_inverse) <- colnames(S_true)

# NNLS solution, one fragment bin at a time.
S_nnls <- matrix(0, nrow = 2, ncol = ncol(Y))
rownames(S_nnls) <- c("Leu", "Ile")
colnames(S_nnls) <- colnames(S_true)
rss <- numeric(ncol(Y))

for (j in seq_len(ncol(Y))) {
  fit <- nnls_small(A, Y[, j])
  S_nnls[, j] <- fit$x
  rss[j] <- fit$rss
}

S_true_norm <- normalize_rows(S_true)
S_inverse_norm <- normalize_rows(pmax(S_inverse, 0))
S_nnls_norm <- normalize_rows(S_nnls)

# Control: if MS1 is treated as one merged same-m/z EIC, separation is impossible.
A_merged <- matrix(approx(rt, merged_eic, xout = ms2_rt, rule = 2)$y, ncol = 1)
rownames(A_merged) <- scan_label
colnames(A_merged) <- "Merged"

S_merged <- numeric(ncol(Y))
for (j in seq_len(ncol(Y))) {
  fit <- nnls_small(A_merged, Y[, j])
  S_merged[j] <- fit$x[1]
}
S_merged_norm <- S_merged / max(S_merged)

corr_inverse_leu <- cor(S_true_norm["Leu", ], S_inverse_norm["Leu", ])
corr_inverse_ile <- cor(S_true_norm["Ile", ], S_inverse_norm["Ile", ])
corr_nnls_leu <- cor(S_true_norm["Leu", ], S_nnls_norm["Leu", ])
corr_nnls_ile <- cor(S_true_norm["Ile", ], S_nnls_norm["Ile", ])
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
    scan = scan_label,
    rt_min = ms2_rt,
    A_leucine = A[, "Leu"],
    A_isoleucine = A[, "Ile"],
    leu_fraction = A[, "Leu"] / rowSums(A),
    ile_fraction = A[, "Ile"] / rowSums(A)
  ),
  file.path(out_dir, "design_matrix_A.csv"),
  row.names = FALSE
)

estimated_spectra <- data.frame(
  fragment_mz = fragment_mz,
  true_leucine = as.numeric(S_true_norm["Leu", ]),
  inverse_leucine = as.numeric(S_inverse_norm["Leu", ]),
  nnls_leucine = as.numeric(S_nnls_norm["Leu", ]),
  true_isoleucine = as.numeric(S_true_norm["Ile", ]),
  inverse_isoleucine = as.numeric(S_inverse_norm["Ile", ]),
  nnls_isoleucine = as.numeric(S_nnls_norm["Ile", ]),
  wrong_merged_estimate = as.numeric(S_merged_norm)
)
write.csv(
  estimated_spectra,
  file.path(out_dir, "estimated_spectra.csv"),
  row.names = FALSE
)

observed_ms2 <- as.data.frame(Y)
observed_ms2$scan <- scan_label
observed_ms2$rt_min <- ms2_rt
write.csv(
  observed_ms2,
  file.path(out_dir, "observed_chimeric_ms2.csv"),
  row.names = FALSE
)

summary_lines <- c(
  paste("Output directory:", out_dir),
  "DDA-MS2 scan count: 2",
  paste("Condition number of A:", round(condition_number, 3)),
  paste("Correlation true vs inverse Leu:", round(corr_inverse_leu, 4)),
  paste("Correlation true vs inverse Ile:", round(corr_inverse_ile, 4)),
  paste("Correlation true vs NNLS Leu:", round(corr_nnls_leu, 4)),
  paste("Correlation true vs NNLS Ile:", round(corr_nnls_ile, 4)),
  "",
  "Design matrix A:",
  paste(capture.output(print(A)), collapse = "\n"),
  "",
  "Estimated spectra:",
  paste(capture.output(print(estimated_spectra, row.names = FALSE)), collapse = "\n")
)
writeLines(summary_lines, file.path(out_dir, "summary.txt"))

# Plot 1: MS1 EIC and two DDA-MS2 scan times.
write_plot("01_ms1_eic_and_dda_times.png", {
  plot(rt, leu_eic, type = "l", lwd = 2, col = "#1f77b4",
       xlab = "Retention time (min)", ylab = "MS1 intensity",
       main = "Two DDA-MS2 spectra at Leu/Ile peak tops",
       ylim = c(0, max(merged_eic) * 1.08))
  lines(rt, ile_eic, lwd = 2, col = "#d62728")
  lines(rt, merged_eic, lwd = 2, col = "gray25", lty = 2)
  points(ms2_rt, approx(rt, merged_eic, xout = ms2_rt, rule = 2)$y, pch = 16, cex = 1.0)
  text(ms2_rt, approx(rt, merged_eic, xout = ms2_rt, rule = 2)$y,
       labels = c("MS2-1", "MS2-2"), pos = 3, cex = 0.8)
  legend("topright",
         legend = c("Leu-like MS1 component", "Ile-like MS1 component",
                    "Merged same-m/z EIC", "DDA-MS2 scan times"),
         col = c("#1f77b4", "#d62728", "gray25", "black"),
         lty = c(1, 1, 2, NA), pch = c(NA, NA, NA, 16), lwd = c(2, 2, 2, NA),
         bty = "n")
})

# Plot 2: 2 x 2 design matrix A.
write_plot("02_design_matrix_A.png", {
  image(
    x = seq_len(ncol(A)),
    y = seq_len(nrow(A)),
    z = t(A / max(A)),
    axes = FALSE,
    xlab = "MS1 component",
    ylab = "DDA-MS2 scan",
    main = "2 x 2 MS1-derived mixing matrix A",
    col = hcl.colors(64, "Viridis")
  )
  axis(1, at = c(1, 2), labels = c("Leu", "Ile"))
  axis(2, at = c(1, 2), labels = c("Leu apex", "Ile apex"), las = 1)
  text(rep(1:2, each = 2), rep(1:2, times = 2),
       labels = sprintf("%.2f", as.vector(A / max(A))), col = "white")
  box()
})

# Plot 3: true vs NNLS-estimated spectra.
write_plot("03_true_vs_estimated_spectra.png", expr = {
  old_par <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  on.exit(par(old_par), add = TRUE)
  for (compound in c("Leu", "Ile")) {
    barplot(
      rbind(S_true_norm[compound, ], S_nnls_norm[compound, ]),
      beside = TRUE,
      col = c("gray35", "#2ca02c"),
      names.arg = fragment_mz,
      ylim = c(0, 1.12),
      ylab = "Relative intensity",
      xlab = "Fragment m/z bin",
      main = compound
    )
    legend("topright", legend = c("True", "Estimated by NNLS"),
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
    main = "One merged MS1 feature gives only a mixed spectrum"
  )
  legend("topright",
         legend = c("True Leu", "True Ile", "Estimated from merged same-m/z EIC"),
         fill = c("#1f77b4", "#d62728", "gray35"), bty = "n")
})

cat(paste(summary_lines, collapse = "\n"))
cat("\n")
