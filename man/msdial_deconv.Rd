\name{msdial_deconv}
\alias{msdial_deconv}
\title{Perform simplified MS-DIAL-like MS2 deconvolution}
\description{
Perform a simplified MS-DIAL-like deconvolution step using model chromatograms
stored in \code{msdial$result$model}. Fragment chromatograms are fit by least
squares against model chromatograms plus simple baseline terms. This is
MS-DIAL-like and not an exact reimplementation of MS-DIAL.
}
\usage{
msdial_deconv(msdial)
}
\arguments{
  \item{msdial}{
    An object containing \code{Y}, \code{mz}, \code{focused_peak}, and
    \code{msdial$result$model}.
  }
}
\value{
The input \code{msdial} object with \code{msdial$result$deconv} added.
}
\details{
The returned \code{msdial$result$deconv} list contains:
\itemize{
  \item \code{spectra}: deconvoluted component spectra.
  \item \code{model_matrix}: matrix of model chromatograms used for fitting.
  \item \code{window}: scan indices used in the least-squares fit.
}
}
\examples{
\dontrun{
data(msdial)
msdial <- msdial_model(msdial)
msdial <- msdial_deconv(msdial)
msdial$result$deconv$spectra
}
}
\seealso{
\code{\link{msdial_model}}
}
\keyword{deconvolution}
\keyword{mass spectrometry}
