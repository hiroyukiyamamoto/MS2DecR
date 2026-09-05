\name{optimizeJointALS}
\alias{optimizeJointALS}
\title{Optimize a joint RT/IM ALS model}
\description{
Optimizes chromatographic profiles, ion-mobility profiles, and shared MS/MS
spectra for two matrices with common fragment columns.
}
\usage{
optimizeJointALS(X, Y, C, M, com = ncol(C), lambda = 1e-06, maxiter = 500)
}
\arguments{
  \item{X}{Numeric matrix with retention time points in rows and MS/MS fragment m/z values in columns.}
  \item{Y}{Numeric matrix with ion-mobility points in rows and the same MS/MS fragment m/z columns as \code{X}.}
  \item{C}{Initial chromatographic component profile matrix.}
  \item{M}{Initial ion-mobility component profile matrix.}
  \item{com}{Number of components.}
  \item{lambda}{Non-negative ridge regularization parameter.}
  \item{maxiter}{Maximum number of ALS iterations.}
}
\value{
A list containing optimized \code{A}, \code{C}, \code{M}, error history
\code{E}, \code{best_iteration}, \code{fitted_X}, and \code{fitted_Y}.
}
\details{
The model is \code{X = C A} and \code{Y = M A}, where \code{A} is a shared
MS/MS spectrum matrix. Non-negativity is applied to \code{A}, \code{C}, and
\code{M} after each least-squares update.
}
\seealso{
\code{\link{jointALS}}, \code{\link{msdial_joint}}
}
