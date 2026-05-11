\name{msdial_model}
\alias{msdial_model}
\title{Build simplified MS-DIAL-like model chromatograms}
\description{
Build model chromatograms from an \code{msdial} object by detecting one
matched-filter-based peak candidate from each fragment chromatogram, grouping
nearby apex scans, and selecting representative model chromatograms. This is a
MS-DIAL-like procedure and not an exact reimplementation of MS-DIAL.
}
\usage{
msdial_model(msdial, int_min = 300, detect_fwhm = 5, ideal_min = 0.95,
  sharp_min = 0.02, scan_gap = 1L)
}
\arguments{
  \item{msdial}{
    A list-like object containing \code{Y}, \code{mz}, \code{rt}, and \code{focused_peak}.
  }
  \item{int_min}{
    Minimum apex intensity required for a fragment peak candidate.
  }
  \item{detect_fwhm}{
    Matched filter width used for peak detection within each fragment chromatogram.
  }
  \item{ideal_min}{
    Minimum ideal slope value required for a candidate peak.
  }
  \item{sharp_min}{
    Minimum sharpness value required for a candidate peak.
  }
  \item{scan_gap}{
    Maximum apex scan difference allowed when merging candidates into the same group.
  }
}
\value{
The input \code{msdial} object with \code{msdial$result$model} added.
}
\details{
The returned \code{msdial$result$model} list contains:
\itemize{
  \item \code{params}: model-building parameters.
  \item \code{candidates}: grouped candidate peaks.
  \item \code{triplet_summary}: summary of selected \code{M1}, \code{M2}, and \code{M3}.
  \item \code{model_triplet}: model chromatogram objects used downstream.
}
}
\examples{
\dontrun{
data(msdial)
msdial <- msdial_model(msdial)
msdial$result$model$triplet_summary
}
}
\seealso{
\code{\link{msdial_deconv}}
}
\keyword{models}
\keyword{mass spectrometry}
