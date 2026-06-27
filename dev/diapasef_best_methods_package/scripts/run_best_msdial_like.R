options(stringsAsFactors = FALSE)

script_path <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(script_path)) {
  dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

source(normalizePath(file.path(script_dir, "..", "..", "test_diapasef_deconv.R"), winslash = "/", mustWork = TRUE))

cat("\nBest MS-DIAL-like result:\n")
cat(normalizePath(file.path(script_dir, "..", "results", "images", "best_msdial_like_window_baseline_summary.png"), winslash = "/", mustWork = FALSE), "\n")
cat("Purity table:\n")
cat(normalizePath(file.path(script_dir, "..", "results", "tables", "ms2decr_rt070_im060_window_baseline_commonA_purity.csv"), winslash = "/", mustWork = FALSE), "\n")
