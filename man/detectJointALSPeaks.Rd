\name{detectJointALSPeaks}
\alias{detectJointALSPeaks}
\title{Detect valid components from joint ALS profiles}
\description{
Filters joint ALS components using matched-filter apex consistency in retention
time and ion-mobility profiles.
}
\usage{
detectJointALSPeaks(C, M, A, c_wid = 5, m_wid = 5,
                    c_gap = 5, m_gap = 5,
                    mode = c("both", "either", "rt", "im"))
}
\arguments{
  \item{C}{Chromatographic component profile matrix.}
  \item{M}{Ion-mobility component profile matrix.}
  \item{A}{Shared MS/MS spectrum matrix.}
  \item{c_wid, m_wid}{Matched-filter widths for RT and IM profiles.}
  \item{c_gap, m_gap}{Allowed gaps between raw and filtered apex positions.}
  \item{mode}{Component selection mode. \code{"both"} requires both RT and IM
  profiles to pass, \code{"either"} requires either profile to pass,
  \code{"rt"} checks only RT, and \code{"im"} checks only IM.}
}
\value{
A list containing filtered matrices \code{C_fin}, \code{M_fin}, \code{A_fin},
the retained component indices \code{keep_idx}, and a \code{component_table}.
}
\details{
This function is the joint RT/IM counterpart of \code{\link{detectPeaks}}.
For each component, it applies \code{\link{matchedfilter}} to the RT and IM
profiles and checks whether the filtered apex is close to the original apex.
Components with no nonzero MS/MS intensity are removed.
}
\seealso{
\code{\link{jointALS}}, \code{\link{detectPeaks}}
}
