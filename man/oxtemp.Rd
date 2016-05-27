\name{oxtemp}
\alias{oxtemp}
\docType{data}
\title{ Oxford Temperature Data }
\description{
  Annual maximum temperatures collected at Oxford, UK.
}
\usage{data(oxtemp)}
\format{
  A data frame with 80 observations on the following 2 variables.
  \describe{
    \item{maxtemp}{Annual maximum temperatures (in degrees Fahrenheit). }
    \item{year}{The values 1901 to 1980. }
  }
}
\details{
  The data were collected from 1901 to 1980.


}
% zz:
\source{
  Unknown.


}
% \references{
% }
\examples{
\dontrun{ fit <- vglm(maxtemp ~ 1, gevff, data = oxtemp, trace = TRUE) }
}
\keyword{datasets}

