\name{process_ms2_data}
\alias{process_ms2_data}
\title{Process MS2 Data for Deconvolution}
\description{
This function processes MS2 data for a given peak, applies ICA-based deconvolution, and returns the deconvolution results and m/z values.
}
\usage{
process_ms2_data(swath_data, peak, bin_size = 0.01, com = 5, lambda = 0.1, 
                 maxiter = 50, wid = 5, gap = 5)
}
\arguments{
  \item{swath_data}{
    An \code{MSnExp} object containing the SWATH MS2 data.
  }
  \item{peak}{
    A named vector or list with keys \code{"mz"}, \code{"rtmin"}, and \code{"rtmax"} specifying the target m/z and retention time range.
  }
  \item{bin_size}{
    Size of bins for m/z values. Default is \code{0.01}.
  }
  \item{com}{
    Number of components for ICA deconvolution. Default is \code{5}.
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
  \item \code{deconv}: A list with deconvolution results, including loadings and components.
  \item \code{mz}: A vector of m/z values.
}
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example usage:
swath_data <- readMSData("example.mzML", mode = "onDisk")
peak <- list(mz = 500.2, rtmin = 300, rtmax = 400)
result <- process_ms2_data(swath_data, peak)
}
\seealso{
\code{\link[MSnbase]{filterRt}}, \code{\link[MSnbase]{bin}}
}
\keyword{deconvolution}
\keyword{MS2}
\keyword{SWATH}
