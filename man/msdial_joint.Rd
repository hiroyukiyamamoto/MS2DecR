\name{msdial_joint}
\alias{msdial_joint}
\title{Match RT and ion mobility models by preliminary MS/MS similarity}
\description{
Matches model chromatograms and model mobilograms produced by
\code{\link{msdial_joint_model}} using cosine similarity between preliminary
MS/MS spectra. The selected RT and ion-mobility model counts must match.
}
\usage{
msdial_joint(model, rt_models = c("center", "candidates"),
  im_models = c("candidates", "center"))
}
\arguments{
  \item{model}{Result from \code{\link{msdial_joint_model}}.}
  \item{rt_models}{Use RT center models or all RT candidates.}
  \item{im_models}{Use ion-mobility center models or all ion-mobility candidates.}
}
\value{
A list containing aligned \code{C} and \code{M} matrices, matched preliminary
spectra, a similarity matrix, and an alignment table.
}
\examples{
\dontrun{
joint <- msdial_joint(model, rt_models = "center", im_models = "candidates")
joint$alignment
}
}
\seealso{
\code{\link{msdial_joint_model}}, \code{\link{msdial_joint_deconv}}
}
\keyword{models}
\keyword{mass spectrometry}
