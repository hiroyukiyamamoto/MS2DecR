# Perform Independent Component Analysis (ICA) and preprocess the components
performICA <- function(Y0, com) {
  # Apply ICA to the transposed data
  icacom <- ica(t(Y0), com)
  
  # Initialize the component matrix
  C <- matrix(NA, com, nrow(Y0))
  
  # Iterate through each component to align polarity and apply non-negativity constraint
  for (k in 1:com) {
    S <- icacom$M[, k]  # Extract the k-th component
    R <- cor(Y0, S)     # Correlate the component with the original data
    
    # Adjust polarity if the maximum correlation is negative
    if (R[which.max(abs(R))] < 0) {
      S <- -S
    }
    # Apply non-negativity constraint
    S[S < 0] <- 0
    C[k, ] <- S  # Store the component
  }
  
  # Return the processed components and ICA results
  return(list(C = t(C), S = icacom$S, M = icacom$M))
}
