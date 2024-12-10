# Detect peaks in the component matrix based on matched filtering and intensity criteria
detectPeaks <- function(C, A, wid, gap) {
  C_fin <- NULL  # Finalized component matrix
  A_fin <- NULL  # Finalized mixing matrix
  
  # Iterate through each component
  for (i in 1:ncol(C)) {
    # Apply matched filter to smooth the component signal
    filtered_signal <- matchedfilter(C[, i], wid)[[1]]
    
    # Check if the peak positions in the filtered signal and original signal are close
    if (abs(which.max(filtered_signal) - which.max(C[, i])) <= gap && sum(A[i, ] != 0) >= 1) {
      C_fin <- cbind(C_fin, C[, i])  # Add component to the finalized matrix
      A_fin <- rbind(A_fin, A[i, ]) # Add mixing vector to the finalized matrix
    }
  }
  
  # Return the finalized matrices
  return(list(C_fin = C_fin, A_fin = A_fin))
}
