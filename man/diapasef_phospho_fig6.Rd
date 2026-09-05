\name{diapasef_phospho_fig6}
\alias{diapasef_phospho_fig6}
\docType{data}
\title{Extracted dia-PASEF phosphoproteomics Fig. 6 matrices}
\description{
Extracted RT and ion-mobility MS/MS fragment matrices from the HeLa
phosphoproteomics dia-PASEF dataset used in PRIDE project PXD033904. The data
are prepared for examples of \code{\link{msdial_joint_model}},
\code{\link{msdial_joint}}, and \code{\link{msdial_joint_deconv}}.
}
\usage{
data(diapasef_phospho_fig6)
}
\format{
A list with elements:
\describe{
  \item{X}{RT x fragment m/z matrix.}
  \item{Y}{Ion mobility x fragment m/z matrix.}
  \item{mz}{Fragment m/z values.}
  \item{rt}{Retention-time axis in minutes.}
  \item{im}{Ion mobility axis in 1/K0.}
  \item{fragments}{Fragment annotation table, when available.}
  \item{source}{Dataset metadata.}
}
}
\examples{
data(diapasef_phospho_fig6)
dim(diapasef_phospho_fig6$X)
dim(diapasef_phospho_fig6$Y)
}
\keyword{datasets}
