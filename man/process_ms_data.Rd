\name{process_ms_data}
\alias{process_ms_data}
\title{Full Workflow for Processing MS Data}
\description{
Executes the complete workflow for processing MS data, including peak detection, binning, deconvolution, and spectral matching.
}
\usage{
process_ms_data(swath_data, peak_list, sp, params, bin_size = 0.01, com = 5,
                lambda = 0.1, maxiter = 50, wid = 5, gap = 5)
}
\arguments{
  \item{swath_data}{
    An \code{MSnExp} object containing the SWATH data.
  }
  \item{peak_list}{
    A data frame containing the list of peaks with columns \code{"mz"}, \code{"rtmin"}, and \code{"rtmax"}.
  }
  \item{sp}{
    A \code{Spectra} object containing the reference library.
  }
  \item{params}{
    A \code{CentWaveParam} object specifying the parameters for peak detection.
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
# Example usage:
swath_data <- readMSData("example.mzML", mode = "onDisk")
peak_list <- data.frame(mz = c(500.2, 501.3), rtmin = c(300, 400), rtmax = c(350, 450))
params <- CentWaveParam(snthresh = 10, noise = 100, ppm = 5, peakwidth = c(5, 20))
sp <- Spectra("library.msp", source = MsBackendMsp())
result <- process_ms_data(swath_data, peak_list, sp, params)
}
\seealso{
\code{\link[MSnbase]{MSnExp}}, \code{\link[xcms]{CentWaveParam}}
}
\keyword{MS}
\keyword{SWATH}
\keyword{processing}
