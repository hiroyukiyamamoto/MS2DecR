\name{jointALS}
\alias{jointALS}
\title{Joint ALS deconvolution for RT and ion-mobility matrices}
\description{
Runs joint ALS for retention-time and ion-mobility matrices that share MS/MS
fragment columns.
}
\usage{
jointALS(X, Y, com = 2L, C_init = NULL, M_init = NULL,
         lambda = 1e-06, maxiter = 500, align_m = TRUE,
         detect_peaks = TRUE, c_wid = 5, m_wid = 5,
         c_gap = 5, m_gap = 5,
         peak_mode = c("both", "either", "rt", "im"))
}
\arguments{
  \item{X}{Numeric matrix with retention time points in rows and MS/MS fragment m/z values in columns.}
  \item{Y}{Numeric matrix with ion-mobility points in rows and the same MS/MS fragment m/z columns as \code{X}.}
  \item{com}{Number of components.}
  \item{C_init}{Optional initial chromatographic component profile matrix.}
  \item{M_init}{Optional initial ion-mobility component profile matrix.}
  \item{lambda}{Non-negative ridge regularization parameter.}
  \item{maxiter}{Maximum number of ALS iterations.}
  \item{align_m}{Logical. If \code{TRUE}, IM components are reordered by
  preliminary MS/MS spectrum similarity before ALS optimization.}
  \item{detect_peaks}{Logical. If \code{TRUE}, run \code{\link{detectJointALSPeaks}} after ALS.}
  \item{c_wid, m_wid}{Matched-filter widths for RT and IM component profiles.}
  \item{c_gap, m_gap}{Allowed gaps between raw and filtered apex positions for RT and IM profiles.}
  \item{peak_mode}{Peak-selection mode passed to \code{\link{detectJointALSPeaks}}.}
}
\value{
A list containing optimized shared MS/MS spectra \code{A}, chromatographic
profiles \code{C}, ion-mobility profiles \code{M}, error history \code{E},
and fitted matrices. If \code{detect_peaks = TRUE}, the list also contains
\code{peaks} and \code{selected}.
}
\details{
If \code{C_init} and \code{M_init} are supplied, they are used as the initial
RT and IM profiles. If either initial matrix is omitted, component profiles are
initialized by ICA when \pkg{fastICA} is available, with an SVD fallback.
}
\examples{
\dontrun{
data(diapasef_phospho_fig6)
model <- msdial_joint_model(
  diapasef_phospho_fig6$X,
  diapasef_phospho_fig6$Y,
  diapasef_phospho_fig6$mz,
  diapasef_phospho_fig6$rt,
  diapasef_phospho_fig6$im,
  rt_range = c(17.6, 18.6)
)
joint <- msdial_joint(model, rt_models = "center", im_models = "candidates")
fit <- jointALS(
  diapasef_phospho_fig6$X,
  diapasef_phospho_fig6$Y,
  com = ncol(joint$C),
  C_init = joint$C,
  M_init = joint$M
)
}
}
\seealso{
\code{\link{optimizeJointALS}}, \code{\link{msdial_joint}}
}
