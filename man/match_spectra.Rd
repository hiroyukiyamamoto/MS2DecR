\name{match_spectra}
\alias{match_spectra}
\title{Match Experimental Spectra Against a Reference Library}
\description{
This function calculates similarity scores between experimental spectra and a reference library using dot product or other metrics.
}
\usage{
match_spectra(mz, intensity, sp, target_mz, ppm = 10, fun = "dotproduct")
}
\arguments{
  \item{mz}{
    A numeric vector of m/z values for the experimental spectrum.
  }
  \item{intensity}{
    A numeric vector of intensity values for the experimental spectrum.
  }
  \item{sp}{
    A \code{Spectra} object containing the reference library.
  }
  \item{target_mz}{
    Target precursor m/z value for filtering the reference library.
  }
  \item{ppm}{
    Allowed m/z tolerance in parts per million. Default is \code{10}.
  }
  \item{fun}{
    Similarity metric to use. Default is \code{"dotproduct"}.
  }
}
\value{
A numeric value representing the highest similarity score.
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example usage:
library(MsBackendMsp)
sp <- Spectra("library.msp", source = MsBackendMsp())
mz <- c(100.1, 101.1, 102.1)
intensity <- c(10, 20, 30)
target_mz <- 500.2
score <- match_spectra(mz, intensity, sp, target_mz)
}
\seealso{
\code{\link[MSnbase]{Spectra}}, \code{\link[MSnbase]{compareSpectra}}
}
\keyword{spectra}
\keyword{matching}
\keyword{MS2}
