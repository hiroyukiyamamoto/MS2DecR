\name{refitJointALSSpectra}
\alias{refitJointALSSpectra}
\title{Refit shared MS/MS spectra from fixed joint ALS profiles}
\description{
Estimates shared MS/MS spectra while keeping retention-time and ion-mobility
component profiles fixed.
}
\usage{
refitJointALSSpectra(X, Y, C, M, lambda = 1e-08,
                     baseline = c("window", "none"))
}
\arguments{
  \item{X}{Numeric matrix with retention time points in rows and MS/MS fragment
  columns.}
  \item{Y}{Numeric matrix with ion-mobility points in rows and the same MS/MS
  fragment columns as \code{X}.}
  \item{C}{Fixed chromatographic component profile matrix.}
  \item{M}{Fixed ion-mobility component profile matrix.}
  \item{lambda}{Non-negative ridge regularization parameter.}
  \item{baseline}{If \code{"window"}, adds RT and IM slope/intercept nuisance
  terms over nonzero profile windows. If \code{"none"}, solves only for shared
  spectra.}
}
\value{
A list containing the refit spectra as \code{spectra} and \code{A}, fitted
matrices, and, for \code{baseline = "window"}, nuisance terms and the design
matrix.
}
\details{
This function is useful after component selection from overcomplete joint ALS.
It corresponds to the A-only refit step used in the dia-PASEF demonstration.
}
\seealso{
\code{\link{jointALS}}, \code{\link{selectJointALSComponents}},
\code{\link{msdial_joint_deconv}}
}
