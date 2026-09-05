\name{msdial_joint_model}
\alias{msdial_joint_model}
\title{Build MS-DIAL-like RT and ion mobility model profiles}
\description{
Applies the existing \code{\link{msdial_model}} and
\code{\link{msdial_deconv}} workflow independently to prepared retention-time
and ion-mobility matrices. The function returns model chromatograms, model
mobilograms, and neighboring \code{M1}/\code{M2}/\code{M3} candidates for later
matching.
}
\usage{
msdial_joint_model(X, Y, mz, rt, im, rt_range = range(rt),
  im_range = range(im), rt_max_peaks = 5L, im_max_peaks = 5L,
  rt_min_rel_height = 0.02, im_min_rel_height = 0.2,
  rt_detect_fwhm = 5L, im_detect_fwhm = 3L,
  rt_ideal_ladder = c(0.95, 0.9, 0.8, 0.7, 0.6, 0.5),
  im_ideal_ladder = c(0.6, 0.5), int_min = 300, sharp_min = 0.02)
}
\arguments{
  \item{X}{RT x fragment m/z intensity matrix.}
  \item{Y}{Ion mobility x fragment m/z intensity matrix.}
  \item{mz}{Fragment m/z values matching columns of \code{X} and \code{Y}.}
  \item{rt}{Retention-time axis matching rows of \code{X}.}
  \item{im}{Ion-mobility axis matching rows of \code{Y}.}
  \item{rt_range, im_range}{Axis ranges used for peak-region detection.}
  \item{rt_max_peaks, im_max_peaks}{Maximum number of focused peak regions to keep on each axis.}
  \item{rt_min_rel_height, im_min_rel_height}{Minimum relative peak height for focused peak-region detection.}
  \item{rt_detect_fwhm, im_detect_fwhm}{Matched-filter widths for RT and ion mobility axes.}
  \item{rt_ideal_ladder, im_ideal_ladder}{Candidate \code{ideal_min} values tried when calling \code{\link{msdial_model}}.}
  \item{int_min, sharp_min}{Arguments passed to \code{\link{msdial_model}}.}
}
\value{
A list containing center model matrices \code{C} and \code{M}, candidate model
matrices \code{C_candidates} and \code{M_candidates}, preliminary spectra,
diagnostics, and summaries.
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
}
}
\seealso{
\code{\link{msdial_joint}}, \code{\link{msdial_joint_deconv}}
}
\keyword{models}
\keyword{mass spectrometry}
