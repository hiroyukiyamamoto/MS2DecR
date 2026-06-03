\name{optimizeALS}
\alias{optimizeALS}
\title{Alternating Least Squares (ALS) Optimization for Deconvolution}
\description{
Performs ALS optimization to refine loadings and scores in ICA-based deconvolution.
}
\usage{
optimizeALS(X, C, com, lambda, maxiter)
}
\arguments{
  \item{X}{
    A numeric matrix where rows represent samples and columns represent features.
  }
  \item{C}{
    Initial component matrix for the deconvolution process.
  }
  \item{com}{
    Number of components.
  }
  \item{lambda}{
    Regularization parameter for ALS.
  }
  \item{maxiter}{
    Maximum number of iterations for ALS optimization.
  }
}
\value{
A list containing:
\itemize{
  \item \code{A}: Optimized loadings matrix.
  \item \code{C}: Refined component matrix.
  \item \code{error}: Convergence error over iterations.
}
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example usage:
X <- matrix(runif(5000, 0, 1), nrow = 100, ncol = 50)
C <- matrix(runif(500, 0, 1), nrow = 100, ncol = 5)
result <- optimizeALS(X, C, com = 5, lambda = 0.1, maxiter = 50)
}
\keyword{ALS}
\keyword{optimization}
