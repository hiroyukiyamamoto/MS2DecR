\name{detectPeaks}
\alias{detectPeaks}
\title{Peak Detection in MS Spectral Data}
\description{
Detects peaks in MS spectral data using centWave-based peak detection.
}
\usage{
detectPeaks(swath_data, params)
}
\arguments{
  \item{swath_data}{
    An \code{MSnExp} object containing the SWATH data.
  }
  \item{params}{
    A \code{CentWaveParam} object specifying the parameters for peak detection.
  }
}
\value{
An \code{MSnExp} object with detected peaks added to the feature data.
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example usage:
library(xcms)
swath_data <- readMSData("example.mzML", mode = "onDisk")
params <- CentWaveParam(snthresh = 10, noise = 100, ppm = 5, peakwidth = c(5, 20))
result <- detectPeaks(swath_data, params)
}
\seealso{
\code{\link[xcms]{CentWaveParam}}
}
\keyword{peak detection}
\keyword{MS}
\keyword{SWATH}
