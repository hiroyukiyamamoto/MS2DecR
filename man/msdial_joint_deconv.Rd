\name{msdial_joint_deconv}
\alias{msdial_joint_deconv}
\title{Estimate joint MS/MS spectra from matched RT and ion mobility models}
\description{
Uses aligned model chromatograms and model mobilograms from
\code{\link{msdial_joint}} to estimate shared MS/MS spectra by solving
\code{X = C A} and \code{Y = M A} jointly.
}
\usage{
msdial_joint_deconv(joint, X, Y, lambda = 1e-08)
}
\arguments{
  \item{joint}{Result from \code{\link{msdial_joint}}.}
  \item{X}{RT x fragment m/z intensity matrix.}
  \item{Y}{Ion mobility x fragment m/z intensity matrix.}
  \item{lambda}{Ridge regularization strength.}
}
\value{
A list containing aligned \code{C}, \code{M}, estimated spectra \code{A}, and
fitted \code{X}/\code{Y} matrices.
}
\examples{
\dontrun{
fit <- msdial_joint_deconv(joint, diapasef_phospho_fig6$X, diapasef_phospho_fig6$Y)
}
}
\seealso{
\code{\link{msdial_joint_model}}, \code{\link{msdial_joint}}
}
\keyword{models}
\keyword{mass spectrometry}
