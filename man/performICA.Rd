\name{performICA}
\alias{performICA}
\title{Independent Component Analysis (ICA) for Spectral Data}
\description{
Applies ICA to decompose spectral data into independent components.
}
\usage{
performICA(X, n_components = 5)
}
\arguments{
  \item{X}{
    A numeric matrix where rows represent samples and columns represent features.
  }
  \item{n_components}{
    Number of independent components to extract. Default is \code{5}.
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
X <- matrix(runif(5000, 0, 1), nrow = 100, ncol = 50)
result <- performICA(X, n_components = 5)
}
\seealso{
\code{\link[fastICA]{fastICA}}
}
\keyword{ICA}
\keyword{dimensionality reduction}
\keyword{spectral}
