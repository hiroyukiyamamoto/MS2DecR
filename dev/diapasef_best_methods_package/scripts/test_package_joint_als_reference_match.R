repo_dir <- normalizePath(getwd(), mustWork = TRUE)

source(file.path(repo_dir, "R", "jointALS.R"))
load(file.path(repo_dir, "data", "diapasef_phospho_fig6.rda"))

X <- diapasef_phospho_fig6$X
Y <- diapasef_phospho_fig6$Y

assert_true <- function(value, label) {
  if (!isTRUE(value)) stop(label)
}

top_fragment_names <- function(A) {
  colnames(A)[max.col(A, ties.method = "first")]
}

cosine <- function(x, y) {
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(den) || den == 0) return(NA_real_)
  sum(x * y) / den
}

greedy_match <- function(score) {
  out <- integer(nrow(score))
  used <- logical(ncol(score))
  for (i in seq_len(nrow(score))) {
    available <- which(!used)
    best <- available[which.max(score[i, available])]
    out[i] <- best
    used[best] <- TRUE
  }
  out
}

read_reference_top1 <- function(path) {
  x <- read.csv(path, check.names = FALSE)
  x <- x[x$rank == 1L, , drop = FALSE]
  setNames(x$fragment, x$component)
}

out_dir <- file.path(repo_dir, "dev/diapasef_best_methods_package/results/images")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reference_dir <- file.path(repo_dir, "dev/diapasef_best_methods_package/diapasef_deconv_test")
ref_k2_file <- file.path(reference_dir, "mcr_als_commonA_top_fragments.csv")
ref_k5_file <- file.path(reference_dir, "mcr_als_overcomplete_commonA_spectra.csv")
if (!file.exists(ref_k2_file) || !file.exists(ref_k5_file)) {
  cat("Package jointALS reference-match test: SKIP\n")
  cat("Reference result files are not available in this checkout.\n")
  quit(status = 0)
}

ref_k2_top1 <- read_reference_top1(ref_k2_file)

fit_k2 <- jointALS(X, Y, com = 2L, lambda = 1e-6, maxiter = 500,
                   detect_peaks = FALSE)
new_k2_top1 <- top_fragment_names(fit_k2$A)
names(new_k2_top1) <- rownames(fit_k2$A)

ref_over_A <- as.matrix(read.csv(ref_k5_file, row.names = 1, check.names = FALSE))
storage.mode(ref_over_A) <- "double"

fit_k5 <- jointALS(X, Y, com = 5L, lambda = 1e-6, maxiter = 500,
                   detect_peaks = FALSE)
new_k5_top1 <- top_fragment_names(fit_k5$A)
ref_k5_top1 <- top_fragment_names(ref_over_A)

k5_cosine <- outer(seq_len(nrow(ref_over_A)), seq_len(nrow(fit_k5$A)),
                   Vectorize(function(i, j) cosine(ref_over_A[i, ], fit_k5$A[j, ])))
k5_match <- greedy_match(k5_cosine)
k5_spectrum_match <- data.frame(
  previous_component = rownames(ref_over_A),
  new_component = rownames(fit_k5$A)[k5_match],
  cosine_similarity = k5_cosine[cbind(seq_len(nrow(k5_cosine)), k5_match)],
  stringsAsFactors = FALSE
)

comparison <- rbind(
  data.frame(
    case = "k2_ICA_initialized",
    component = names(new_k2_top1),
    previous_top_fragment = unname(ref_k2_top1),
    new_top_fragment = unname(new_k2_top1),
    match = unname(ref_k2_top1) == unname(new_k2_top1),
    unordered_case_match = identical(sort(unname(new_k2_top1)), sort(unname(ref_k2_top1))),
    stringsAsFactors = FALSE
  ),
  data.frame(
    case = "k5_ICA_initialized",
    component = rownames(fit_k5$A),
    previous_top_fragment = unname(ref_k5_top1),
    new_top_fragment = unname(new_k5_top1),
    match = unname(ref_k5_top1) == unname(new_k5_top1),
    unordered_case_match = identical(sort(unname(new_k5_top1)), sort(unname(ref_k5_top1))),
    stringsAsFactors = FALSE
  )
)

write.csv(
  comparison,
  file.path(out_dir, "package_joint_als_reference_match.csv"),
  row.names = FALSE
)
write.csv(
  k5_spectrum_match,
  file.path(out_dir, "package_joint_als_reference_k5_spectrum_similarity.csv"),
  row.names = FALSE
)

assert_true(identical(sort(unname(new_k2_top1)), sort(unname(ref_k2_top1))),
            "jointALS k=2 top fragments should match the previous joint_mcr_als result up to component order.")
assert_true(identical(sort(unname(new_k5_top1)), sort(unname(ref_k5_top1))),
            "jointALS k=5 top fragments should match the previous overcomplete joint_mcr_als result up to component order.")
assert_true(min(k5_spectrum_match$cosine_similarity, na.rm = TRUE) > 0.97,
            "jointALS k=5 spectra should be highly similar to previous spectra after component matching.")

cat("Package jointALS reference-match test: PASS\n")
print(comparison, row.names = FALSE)
cat("k=5 spectrum cosine similarity after component matching:\n")
print(k5_spectrum_match, row.names = FALSE)
