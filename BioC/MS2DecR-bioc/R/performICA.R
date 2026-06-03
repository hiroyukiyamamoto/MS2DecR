# Perform Independent Component Analysis (ICA) and preprocess the components
performICA <- function(Y0, com) {
  # Apply ICA to the transposed data
  icacom <- .ms2decr_ica(t(Y0), com)
  
  # Initialize the component matrix
  C <- matrix(NA, com, nrow(Y0))
  keep <- apply(Y0, 2, sd) > 0
  
  # Iterate through each component to align polarity and apply non-negativity constraint
  for (k in 1:com) {
    S <- icacom$M[, k]  # Extract the k-th component
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
  return(list(C = t(C), S = icacom$S, M = icacom$M))
}

.ms2decr_ica <- function(x, ncomp) {
  ncomp <- min(ncomp, nrow(x), ncol(x))
  fit <- stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = ncomp)
  list(S = fit$x[, seq_len(ncomp), drop = FALSE],
       M = fit$rotation[, seq_len(ncomp), drop = FALSE])
}
