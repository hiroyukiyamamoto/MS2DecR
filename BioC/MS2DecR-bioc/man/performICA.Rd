\name{performICA}
\alias{performICA}
\title{Independent Component Analysis (ICA) for Spectral Data}
\description{
Applies \pkg{fastICA}-based ICA to decompose spectral data into independent
components.
}
\usage{
performICA(Y0, com)
}
\arguments{
  \item{Y0}{
    A numeric matrix where rows represent scans and columns represent m/z bins.
  }
  \item{com}{
    Number of components to extract.
  }
}
\value{
A list containing:
\itemize{
  \item \code{S}: Independent components matrix.
  \item \code{M}: Mixing matrix.
}
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example usage:
Y0 <- matrix(runif(5000, 0, 1), nrow = 100, ncol = 50)
result <- performICA(Y0, com = 5)
}
\seealso{
\code{\link[fastICA]{fastICA}}
}
\keyword{ICA}
\keyword{dimensionality reduction}
\keyword{spectral}
