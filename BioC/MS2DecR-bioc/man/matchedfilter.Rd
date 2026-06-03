\name{matchedfilter}
\alias{matchedfilter}
\title{Apply a Gaussian Matched Filter}
\description{
Applies a second-derivative Gaussian matched filter to a numeric signal.
}
\usage{
matchedfilter(xt, fwhm)
}
\arguments{
  \item{xt}{
    Numeric signal.
  }
  \item{fwhm}{
    Full width at half maximum used to define the filter width.
  }
}
\value{
A list containing the filtered signal and Gaussian filter.
}
\examples{
x <- dnorm(seq(-3, 3, length.out = 50))
res <- matchedfilter(x, fwhm = 3)
plot(res$filtered_signal, type = "l")
}
