\name{process_ms_data}
\alias{process_ms_data}
\title{Full Workflow for Processing MS Data}
\description{
Executes the complete workflow for processing MS data, including peak detection, binning, deconvolution, and spectral matching.
}
\usage{
process_ms_data(swath_data, sp, bin_size = 0.01, com = 5, lambda = 0.1,
                maxiter = 50, wid = 5, gap = 5, snthresh = 5,
                noise = 100, ppm = 10, peakwidth = c(3, 30))
}
\arguments{
  \item{swath_data}{
    An \code{MSnExp} object containing the SWATH data.
  }
  \item{sp}{
    A \code{Spectra} object containing the reference library.
  }
  \item{bin_size}{
    Size of bins for m/z values. Default is \code{0.01}.
  }
  \item{com}{
    Number of components for ICA deconvolution. Default is \code{5}.
  }
  \item{lambda}{
    Regularization parameter for ALS. Default is \code{0.1}.
  }
  \item{maxiter}{
    Maximum number of iterations for ALS. Default is \code{50}.
  }
  \item{wid}{
    Peak width parameter for deconvolution. Default is \code{5}.
  }
  \item{gap}{
    Gap parameter for peak detection. Default is \code{5}.
  }
  \item{snthresh}{
    Signal-to-noise threshold passed to \code{CentWaveParam}.
  }
  \item{noise}{
    Noise threshold passed to \code{CentWaveParam}.
  }
  \item{ppm}{
    Mass accuracy in parts per million.
  }
  \item{peakwidth}{
    Peak width range passed to \code{CentWaveParam}.
  }
}
\value{
A list containing:
\itemize{
  \item \code{peaks}: Detected peaks.
  \item \code{deconv_results}: Deconvolution results.
  \item \code{spectra_matches}: Spectral matching scores.
}
}
\author{
Hiroyuki Yamamoto
}
\examples{
\dontrun{
result <- process_ms_data(swath_data, sp)
}
}
\seealso{
\code{\link[MSnbase]{MSnExp}}, \code{\link[xcms]{CentWaveParam}}
}
\keyword{MS}
\keyword{SWATH}
\keyword{processing}
