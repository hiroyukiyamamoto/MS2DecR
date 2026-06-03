# Apply a Gaussian matched filter to the input signal
matchedfilter <- function(xt, fwhm) {
  
  # Calculate sigma from full width at half maximum (FWHM)
  sigma <- fwhm / 2.3548
  
  # Initialize variables
  M <- length(xt)                      # Length of the input signal
  filtered_signal <- rep(0, M)         # Initialize the filtered signal with zeros
  
  # Define the filter width based on sigma
  m <- ceiling(4 * sigma)              # Define the filter size (4-sigma range)
  
  # Create the Gaussian second derivative filter
  w <- -m:m                            # Window indices
  gaussian_filter <- (1 - w^2 / sigma^2) * exp(-w^2 / (2 * sigma^2))
  
  # Apply the filter to the input signal
  for (tx in (m + 1):(M - m)) {        # Iterate over valid range of indices
    filtered_signal[tx] <- sum(
      xt[(tx - m):(tx + m)] * gaussian_filter
    ) / sqrt(sum(gaussian_filter^2))   # Normalize by filter norm
  }
  
  # Return the filtered signal and the Gaussian filter
  result <- list(
    filtered_signal = filtered_signal, # Filtered signal (Ybis)
    gaussian_filter = gaussian_filter  # Gaussian filter (fbis)
  )
  return(result)
}
