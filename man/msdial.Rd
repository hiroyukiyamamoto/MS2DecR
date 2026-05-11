\name{msdial}
\docType{data}
\alias{msdial}
\title{MS-DIAL-like example data for simplified MS2 deconvolution}
\description{
\code{msdial} is an example dataset for the simplified MS-DIAL-like workflow in
\pkg{MS2DecR}. It contains a precomputed MS2 matrix and related metadata so that
model selection and deconvolution can be tested without re-reading the original
\code{mzML} file.
}
\usage{
data(msdial)
}
\format{
A list with the following components:
\itemize{
  \item \code{Y}: intensity matrix with RT scans in rows and fragment m/z bins in columns.
  \item \code{mz}: numeric vector of fragment m/z values.
  \item \code{rt}: numeric vector of retention times.
  \item \code{focused_peak}: a list with \code{left_idx}, \code{apex_idx}, and \code{right_idx}.
  \item \code{rt_range}: numeric vector of length 2.
  \item \code{result}: an initially empty list for model and deconvolution results.
}
}
\details{
This dataset is intended to be used with the simplified, MS-DIAL-like workflow
implemented by \code{\link{msdial_model}} and \code{\link{msdial_deconv}}.
The workflow is MS-DIAL-like and not an exact reimplementation of MS-DIAL.

The object was created from an \code{mzML} file by:
\itemize{
  \item reading DIA MS data with \code{MSnbase::readMSData()},
  \item extracting an MS2 matrix for a fixed precursor and retention-time range with \code{\link{generate_ms2_matrix}},
  \item storing the resulting \code{Y}, \code{mz}, and \code{rt},
  \item adding a simple \code{focused_peak} definition derived from the fixed retention-time range.
}

The source file name and acquisition settings are strongly consistent with the
human plasma HILIC positive-ion SWATH example described in the original
MS-DIAL paper, including the 25-Da SWATH window and the metoclopramide /
norcocaine deconvolution example shown in Figure 2. Therefore, this dataset is
best understood as a package-ready example corresponding to the MS-DIAL paper
example.

However, \code{msdial} is still a precomputed \pkg{MS2DecR} object, not raw
vendor data and not a direct MS-DIAL export. It is prepared for the simplified
MS-DIAL-like model-selection and deconvolution workflow in \pkg{MS2DecR}.
}
\examples{
data(msdial)
str(msdial, max.level = 1)
}
\keyword{datasets}
