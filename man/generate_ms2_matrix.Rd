\name{generate_ms2_matrix}
\alias{generate_ms2_matrix}
\title{Generate MS2 Data Matrix}
\description{
This function processes MS2 data from a DIA dataset, filtering by m/z and retention time, binning the data, and converting it into a matrix format.
}
\usage{
generate_ms2_matrix(dia_data, premz, rt_range, bin_size = 0.01)
}
\arguments{
  \item{dia_data}{
    An \code{OnDiskMSnExp} object containing the DIA data to be processed.
  }
  \item{premz}{
    A numeric value specifying the target m/z range for filtering the MS2 data.
  }
  \item{rt_range}{
    A numeric vector of length 2 specifying the retention time range (\code{c(rt_min, rt_max)}).
  }
  \item{bin_size}{
    A numeric value specifying the bin size for binning the MS2 data. Default is \code{0.01}.
  }
}
\details{
This function performs the following steps:
\itemize{
  \item Filters the MS2 data by the specified m/z range and retention time range.
  \item Bins the intensity data using the specified bin size.
  \item Converts the binned intensity data into a matrix format where each row corresponds to a spectrum.
  \item Extracts the retention time (\code{rt}) and m/z axis (\code{mz}) values.
}
}
\value{
A list containing the following components:
\itemize{
  \item \code{mz}: A numeric vector of m/z values.
  \item \code{Y}: A matrix where each row corresponds to the intensity values of a spectrum.
  \item \code{rt}: A numeric vector of retention times for the spectra.
}
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example usage:
# Load DIA data (example file path, replace with actual path)
# dia_data <- MSnbase::readMSData("example.mzML", mode = "onDisk")

# Define parameters
premz <- 290.1387
rt_range <- c(170, 190)
bin_size <- 0.01

# Generate MS2 data matrix
ms2_data <- generate_ms2_matrix(dia_data, premz, rt_range, bin_size)

# Access the results
mz <- ms2_data$mz
Y <- ms2_data$Y
rt <- ms2_data$rt

# Inspect the results
print(dim(Y))  # Dimensions of the intensity matrix
print(mz[1:10])  # First 10 m/z values
print(rt[1:10])  # First 10 retention times
}
\seealso{
\code{\link[MSnbase]{filterMsLevel}},
\code{\link[MSnbase]{filterIsolationWindow}},
\code{\link[MSnbase]{filterRt}},
\code{\link[MSnbase]{bin}}
}
\keyword{MS2}
\keyword{DIA}
\keyword{binning}
\keyword{spectral}
