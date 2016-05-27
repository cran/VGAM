\name{marital.nz}
\alias{marital.nz}
\docType{data}
\title{
  New Zealand Marital Data
}
\description{
  Some marital data mainly from a large NZ company collected in the
  early 1990s.

}
\usage{data(marital.nz)}
\format{
  A data frame with 6053 observations on the following 3 variables.
  \describe{
    \item{\code{age}}{a numeric vector, age in years}
    \item{\code{ethnicity}}{a factor with levels \code{European}
    \code{Maori} \code{Other} \code{Polynesian}.
    Only Europeans are included in the data set.
  }

    \item{\code{mstatus}}{a factor with levels
   \code{Divorced/Separated}, \code{Married/Partnered},
   \code{Single}, \code{Widowed}.
  } }
}

\details{
This is a subset of a data set collected from a
self-administered questionnaire administered in a large
New Zealand workforce observational study conducted
during 1992--3. The data were augmented by a second study
consisting of retirees. The data can be considered a
reasonable representation of the white male New Zealand
population in the early 1990s.


}
\source{
  Clinical Trials Research Unit, University of Auckland, New Zealand.


}

\references{
  See \code{\link{bmi.nz}} and \code{\link{chest.nz}}.


}
\examples{
summary(marital.nz)
}
\keyword{datasets}



