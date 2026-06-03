\name{deconvICA}
\alias{deconvICA}
\title{ICA-Based Deconvolution of MS2 Data}
\description{
Performs ICA-based deconvolution on MS2 data to resolve overlapping components.
}
\usage{
deconvICA(Y0, com = 5, lambda = 0.1, maxiter = 50, wid = 0.5, gap = 3)
}
\arguments{
  \item{Y0}{
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
    Peak width parameter used in matched-filter-based peak detection.
  }
  \item{gap}{
    Allowed gap between filtered and original peak tops in peak detection.
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
library(MS2DecR)

data(msdial)

# The msdial dataset is a precomputed MS-DIAL-like example derived from
# the HILIC positive SWATH example corresponding to the original MS-DIAL paper.
res <- deconvICA(
  Y0 = msdial$Y,
  com = 5,
  lambda = 0.1,
  maxiter = 50,
  wid = 5,
  gap = 1
)

res$com
dim(res$ALS$C)
dim(res$ALS$A)
}
\seealso{
\code{\link[Matrix]{Matrix}}, \code{\link[irlba]{irlba}}
}
\keyword{ICA}
\keyword{deconvolution}
