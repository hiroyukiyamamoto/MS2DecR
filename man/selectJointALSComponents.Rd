\name{selectJointALSComponents}
\alias{selectJointALSComponents}
\title{Select components from an overcomplete joint ALS result}
\description{
Scores joint ALS components using retention-time and ion-mobility peak shapes,
then returns selected component matrices.
}
\usage{
selectJointALSComponents(fit, min_keep = 2L, max_keep = 5L,
                         c_wid = 5, c_gap = 5, c_max_local_peaks = 2L,
                         m_wid = 3, m_gap = 5, m_max_local_peaks = 2L,
                         m_smooth_width = 5, m_local_peak_frac = 0.25,
                         m_require_filtered_close = TRUE)
}
\arguments{
  \item{fit}{A result from \code{\link{jointALS}} or \code{\link{optimizeJointALS}}.}
  \item{min_keep}{Minimum number of components to keep. If fewer components pass
  the peak criteria, the highest-scoring components are retained.}
  \item{max_keep}{Maximum number of components to keep.}
  \item{c_wid, m_wid}{Matched-filter widths for RT and IM profiles.}
  \item{c_gap, m_gap}{Allowed gaps between raw and filtered apex positions.}
  \item{c_max_local_peaks, m_max_local_peaks}{Maximum allowed number of local peaks.}
  \item{m_smooth_width}{Smoothing width for IM profile scoring.}
  \item{m_local_peak_frac}{Relative height threshold for counting local IM peaks.}
  \item{m_require_filtered_close}{Logical. If \code{FALSE}, IM selection does not
  require close agreement between raw and filtered apex positions.}
}
\value{
A list containing selected \code{C}, \code{M}, and \code{A} matrices, the ALS
error history, selected indices, and a component scoring table.
}
\details{
This function is intended for overcomplete joint ALS results, for example
\code{jointALS(..., com = 5)}, where more components are estimated than are
expected to represent final analytes.
}
\seealso{
\code{\link{jointALS}}, \code{\link{refitJointALSSpectra}},
\code{\link{detectJointALSPeaks}}
}
