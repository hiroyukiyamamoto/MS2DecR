# Perform Alternating Least Squares (ALS) optimization
optimizeALS <- function(X, C, com, lambda, maxiter) {
  # Initialize lists to store intermediate results and errors
  Alist <- list()
  Clist <- list()
  E <- numeric(maxiter)  # Error for each iteration
  
  # Iterative optimization loop
  for (k in 1:maxiter) {
    # Update A (mixing matrix) using least squares with L2 regularization
    A <- solve(t(C) %*% C + lambda * diag(1, com)) %*% t(C) %*% X
    A[A < 0] <- 0  # Apply non-negativity constraint
    
    # Normalize each row of A
    for (i in 1:com) {
      if (sum(A[i, ] != 0) > 0) {
        A[i, ] <- A[i, ] / sqrt(sum(A[i, ]^2))
      }
    }
    
    # Update C (component matrix) using least squares with L2 regularization
    C <- X %*% t(A) %*% solve(A %*% t(A) + lambda * diag(1, com))
    C[C < 0] <- 0  # Apply non-negativity constraint
    
    # Calculate Frobenius norm of the error
    E[k] <- norm(X - C %*% A, type = "F")
    
    # Store intermediate results
    Alist[[k]] <- A
    Clist[[k]] <- C
  }
  
  # Find the iteration with the minimum error
  min_index <- which.min(E)
  
  # Return the optimized matrices and error history
  return(list(A = Alist[[min_index]], C = Clist[[min_index]], E = E))
}
