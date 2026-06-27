options(stringsAsFactors = FALSE)

script_path <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(script_path)) {
  dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

source(normalizePath(file.path(script_dir, "..", "..", "test_diapasef_deconv.R"), winslash = "/", mustWork = TRUE))

cat("\nBest MCR-ALS result:\n")
cat(normalizePath(file.path(script_dir, "..", "results", "images", "best_mcr_als_peak_selected_Aonly_summary.png"), winslash = "/", mustWork = FALSE), "\n")
cat("Strategy comparison table:\n")
cat(normalizePath(file.path(script_dir, "..", "results", "tables", "mcr_als_strategy_comparison.csv"), winslash = "/", mustWork = FALSE), "\n")
