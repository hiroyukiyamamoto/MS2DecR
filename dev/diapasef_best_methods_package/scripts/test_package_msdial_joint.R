repo_dir <- normalizePath(getwd(), mustWork = TRUE)

source(file.path(repo_dir, "R", "matchedfilter.R"))
source(file.path(repo_dir, "R", "msdial_model.R"))
source(file.path(repo_dir, "R", "msdial_deconv.R"))
source(file.path(repo_dir, "R", "msdial_joint.R"))

data_file <- file.path(repo_dir, "data", "diapasef_phospho_fig6.rda")
load(data_file)

assert_true <- function(value, label) {
  if (!isTRUE(value)) stop(label)
}

assert_error <- function(expr, pattern, label) {
  err <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e
  )
  if (is.null(err)) stop(sprintf("%s: expected an error.", label))
  if (!grepl(pattern, conditionMessage(err))) {
    stop(sprintf("%s: unexpected error: %s", label, conditionMessage(err)))
  }
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

assert_true(ncol(model$C) == 2L, "msdial_joint_model should return two RT center models.")
assert_true(ncol(model$M) == 1L, "msdial_joint_model should return one IM center model.")
assert_true(ncol(model$M_candidates) == 2L, "msdial_joint_model should keep two IM candidates.")
assert_close(model$diagnostics$rt_peaks$apex_value[1], 18.43304, 0.08, "RT peak 1")
assert_close(model$diagnostics$rt_peaks$apex_value[2], 17.81283, 0.08, "RT peak 2")
assert_close(model$diagnostics$im_peaks$apex_value[1], 0.9410, 0.006, "IM peak")

assert_error(
  msdial_joint(model, rt_models = "center", im_models = "center"),
  "must match for msdial_joint",
  "center RT/IM model counts should fail"
)

joint <- msdial_joint(model, rt_models = "center", im_models = "candidates")
fit <- msdial_joint_deconv(joint, diapasef_phospho_fig6$X, diapasef_phospho_fig6$Y)

top_mz <- diapasef_phospho_fig6$mz[max.col(fit$A, ties.method = "first")]
assert_close(top_mz[1], 654.3333, 0.02, "component 1 top m/z")
assert_close(top_mz[2], 694.3114, 0.02, "component 2 top m/z")
assert_true(any(grepl("IM peak 1 M2 0.9410", joint$alignment$im_model)),
            "joint alignment should include IM M2.")
assert_true(any(grepl("IM peak 1 M3 0.9635", joint$alignment$im_model)),
            "joint alignment should include IM M3.")

out_dir <- file.path(repo_dir, "dev/diapasef_best_methods_package/results/images")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(joint$alignment,
          file.path(out_dir, "package_msdial_joint_alignment.csv"),
          row.names = FALSE)
write.csv(fit$A,
          file.path(out_dir, "package_msdial_joint_spectra.csv"))

cat("Package msdial_joint demo-data test: PASS\n")
cat("Alignment:\n")
print(joint$alignment, row.names = FALSE)
cat("Top m/z:", paste(top_mz, collapse = ", "), "\n")
