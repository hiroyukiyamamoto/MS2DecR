# Main function for ICA-based deconvolution with ALS optimization and peak detection
deconvICA <- function(Y0, com = 5, lambda = 0.1, maxiter = 50, wid = 0.5, gap = 3) {
  # Step 1: Perform ICA to initialize the component matrix
  ica_result <- performICA(Y0, com)
  C <- ica_result$C  # Initial component matrix
  X <- Y0            # Input signal matrix
  
  # Step 2: Optimize component and mixing matrices using ALS
  als_result <- optimizeALS(X, C, com, lambda, maxiter)
  A <- als_result$A  # Optimized mixing matrix
  C <- als_result$C  # Optimized component matrix
  
  # Step 3: Detect peaks in the optimized component matrix
  peak_result <- detectPeaks(C, A, wid, gap)
  
  # Step 4: Aggregate results
  deconv <- list()
  deconv$ALS$C <- peak_result$C_fin  # Finalized component matrix
  deconv$ALS$A <- peak_result$A_fin  # Finalized mixing matrix
  deconv$ICA$A <- ica_result$S       # Mixing matrix from ICA
  deconv$ICA$C <- ica_result$M       # Component matrix from ICA
  deconv$com <- nrow(peak_result$A_fin)  # Number of finalized components
  deconv$E <- als_result$E               # Error history
  
  # Return the deconvolution results
  return(deconv)
}
