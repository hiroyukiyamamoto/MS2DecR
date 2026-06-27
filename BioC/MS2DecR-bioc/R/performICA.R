# Perform Independent Component Analysis (ICA) and preprocess the components
performICA <- function(Y0, com) {
  if (!is.matrix(Y0)) {
    Y0 <- as.matrix(Y0)
  }
  if (com < 1L) {
    stop("com must be a positive integer.")
  }
  if (com > min(dim(Y0))) {
    stop("com must be no larger than the smaller dimension of Y0.")
  }

  if (!requireNamespace("fastICA", quietly = TRUE)) {
    stop("The fastICA package is required for performICA(). Please install fastICA.")
  }

  # Estimate time-domain component profiles by fastICA.
  icacom <- fastICA::fastICA(t(Y0), n.comp = com)
  component_profiles <- t(icacom$A)
  source_scores <- icacom$S
  
  # Initialize the component matrix
  C <- matrix(NA, com, nrow(Y0))
  keep <- apply(Y0, 2, sd) > 0
  
  # Iterate through each component to align polarity and apply non-negativity constraint
  for (k in 1:com) {
    S <- component_profiles[, k]  # Extract the k-th component profile
    R <- cor(Y0[, keep, drop = FALSE], S)  # Correlate against varying m/z traces only
    R <- R[is.finite(R)]
    
    # Adjust polarity if the maximum correlation is negative
    if (length(R) && R[which.max(abs(R))] < 0) {
      S <- -S
    }
    # Apply non-negativity constraint
    S[S < 0] <- 0
    C[k, ] <- S  # Store the component
  }
  
  # Return the processed components and ICA results
  return(list(C = t(C), S = source_scores, M = component_profiles))
}
