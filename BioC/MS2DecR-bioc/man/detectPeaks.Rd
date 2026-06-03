\name{detectPeaks}
\alias{detectPeaks}
\title{Peak Detection in MS Spectral Data}
\description{
Filters deconvolved components using a matched filter and peak-position criteria.
}
\usage{
detectPeaks(C, A, wid, gap)
}
\arguments{
  \item{C}{
    Component chromatogram matrix with components in columns.
  }
  \item{A}{
    Component spectrum matrix with components in rows.
  }
  \item{wid}{
    Matched filter width.
  }
  \item{gap}{
    Maximum allowed scan gap between filtered and original peak tops.
  }
}
\value{
A list with filtered component chromatograms \code{C_fin} and spectra \code{A_fin}.
}
\author{
Hiroyuki Yamamoto
}
\examples{
C <- matrix(runif(100), nrow = 20, ncol = 5)
A <- matrix(runif(50), nrow = 5, ncol = 10)
result <- detectPeaks(C, A, wid = 3, gap = 5)
}
\keyword{peak detection}
\keyword{MS}
\keyword{SWATH}
