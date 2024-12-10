\name{deconvICA}
\alias{deconvICA}
\title{ICA-Based Deconvolution of MS2 Data}
\description{
Performs ICA-based deconvolution on MS2 data to resolve overlapping components.
}
\usage{
deconvICA(Y, com = 5, lambda = 0.1, maxiter = 50, wid = 5, gap = 5)
}
\arguments{
  \item{Y}{
    A numeric matrix where rows represent time points and columns represent m/z bins.
  }
  \item{com}{
    Number of components to deconvolve. Default is \code{5}.
  }
  \item{lambda}{
    Regularization parameter for ALS optimization. Default is \code{0.1}.
  }
  \item{maxiter}{
    Maximum number of iterations for ALS optimization. Default is \code{50}.
  }
  \item{wid}{
    Peak width parameter for deconvolution. Default is \code{5}.
  }
  \item{gap}{
    Gap parameter for peak detection. Default is \code{5}.
  }
}
\value{
A list containing:
\itemize{
  \item \code{ALS}: ALS deconvolution results, including loadings and scores.
  \item \code{com}: Number of components resolved.
}
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example usage:
Y <- matrix(runif(5000, 0, 1), nrow = 100, ncol = 50)
result <- deconvICA(Y, com = 5)
}
\seealso{
\code{\link[Matrix]{Matrix}}, \code{\link[irlba]{irlba}}
}
\keyword{ICA}
\keyword{deconvolution}
