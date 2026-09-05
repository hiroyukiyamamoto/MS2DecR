repo_dir <- normalizePath(getwd(), mustWork = TRUE)

source(file.path(repo_dir, "R", "matchedfilter.R"))
source(file.path(repo_dir, "R", "msdial_model.R"))
source(file.path(repo_dir, "R", "msdial_deconv.R"))
source(file.path(repo_dir, "R", "msdial_joint.R"))
source(file.path(repo_dir, "R", "jointALS.R"))

data_file <- file.path(repo_dir, "data", "diapasef_phospho_fig6.rda")
load(data_file)

assert_true <- function(value, label) {
  if (!isTRUE(value)) stop(label)
}

assert_close <- function(value, expected, tol, label) {
  if (abs(value - expected) > tol) {
    stop(sprintf("%s: got %.6f, expected %.6f +/- %.6f", label, value, expected, tol))
  }
}

model <- msdial_joint_model(
  X = diapasef_phospho_fig6$X,
  Y = diapasef_phospho_fig6$Y,
  mz = diapasef_phospho_fig6$mz,
  rt = diapasef_phospho_fig6$rt,
  im = diapasef_phospho_fig6$im,
  rt_range = c(17.6, 18.6),
  im_range = range(diapasef_phospho_fig6$im),
  rt_max_peaks = 2L,
  im_max_peaks = 2L,
  rt_min_rel_height = 0.02,
  im_min_rel_height = 0.2
)
joint <- msdial_joint(model, rt_models = "center", im_models = "candidates")

fit <- jointALS(
  X = diapasef_phospho_fig6$X,
  Y = diapasef_phospho_fig6$Y,
  com = ncol(joint$C),
  C_init = joint$C,
  M_init = joint$M,
  lambda = 1e-6,
  maxiter = 50,
  align_m = FALSE,
  detect_peaks = TRUE
)

assert_true(ncol(fit$C) == 2L, "jointALS should keep two RT components.")
assert_true(ncol(fit$M) == 2L, "jointALS should keep two IM components.")
assert_true(nrow(fit$A) == 2L, "jointALS should estimate two shared MS/MS spectra.")
assert_true(ncol(fit$A) == length(diapasef_phospho_fig6$mz),
            "jointALS spectra should have one column per fragment.")
assert_true(all(fit$A >= 0) && all(fit$C >= 0) && all(fit$M >= 0),
            "jointALS should keep non-negative A, C, and M.")
assert_true(tail(fit$E, 1) <= fit$E[1],
            "jointALS should not increase the reconstruction error.")
assert_true(!is.null(fit$peaks), "jointALS should return peak-detection results.")
assert_true(!is.null(fit$selected), "jointALS should return selected components.")
assert_true(all(c("c_pass", "m_pass", "keep") %in% colnames(fit$peaks$component_table)),
            "jointALS peak detection should report RT and IM pass flags.")

top_mz <- diapasef_phospho_fig6$mz[max.col(fit$A, ties.method = "first")]
assert_true(all(is.finite(top_mz)), "jointALS top m/z values should be finite.")

fit_k5 <- jointALS(
  X = diapasef_phospho_fig6$X,
  Y = diapasef_phospho_fig6$Y,
  com = 5L,
  lambda = 1e-6,
  maxiter = 50,
  detect_peaks = TRUE
)
assert_true(ncol(fit_k5$C) == 5L, "jointALS k=5 should keep raw five components.")
assert_true(nrow(fit_k5$peaks$component_table) == 5L,
            "jointALS k=5 should evaluate five components.")
assert_true(ncol(fit_k5$selected$C) <= 5L,
            "jointALS selected k=5 components should not exceed five.")

selected_k5 <- selectJointALSComponents(
  fit_k5,
  min_keep = 2L,
  max_keep = 5L,
  m_gap = 8,
  m_max_local_peaks = 6L,
  m_smooth_width = 7,
  m_local_peak_frac = 0.50,
  m_require_filtered_close = FALSE
)
refit_k5 <- refitJointALSSpectra(
  X = diapasef_phospho_fig6$X,
  Y = diapasef_phospho_fig6$Y,
  C = selected_k5$C,
  M = selected_k5$M,
  baseline = "window"
)
assert_true(ncol(selected_k5$C) >= 2L, "selectJointALSComponents should keep at least two components.")
assert_true(nrow(refit_k5$A) == ncol(selected_k5$C),
            "refitJointALSSpectra should estimate one spectrum per selected component.")
assert_true(ncol(refit_k5$A) == length(diapasef_phospho_fig6$mz),
            "refitJointALSSpectra spectra should have one column per fragment.")

out_dir <- file.path(repo_dir, "dev/diapasef_best_methods_package/results/images")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  data.frame(iteration = seq_along(fit$E), error = fit$E),
  file.path(out_dir, "package_joint_als_error.csv"),
  row.names = FALSE
)
write.csv(
  fit$A,
  file.path(out_dir, "package_joint_als_spectra.csv")
)
write.csv(
  fit_k5$peaks$component_table,
  file.path(out_dir, "package_joint_als_k5_component_table.csv"),
  row.names = FALSE
)
write.csv(
  selected_k5$component_table,
  file.path(out_dir, "package_joint_als_k5_selection_table.csv"),
  row.names = FALSE
)
write.csv(
  refit_k5$A,
  file.path(out_dir, "package_joint_als_k5_selected_refit_spectra.csv")
)

cat("Package jointALS demo-data test: PASS\n")
cat(sprintf("Error: %.6f -> %.6f\n", fit$E[1], tail(fit$E, 1)))
cat("Top m/z:", paste(top_mz, collapse = ", "), "\n")
cat("k=5 selected components:", ncol(fit_k5$selected$C), "\n")
cat("k=5 advanced selected components:", ncol(selected_k5$C), "\n")
