\name{auuc}
\alias{auuc}
\docType{data}
\title{ Auckland University Undergraduate Counts Data}
\description{
Undergraduate student enrolments
at the University of Auckland in 1990.
}
\usage{data(auuc)}
\format{
  A data frame with 4 observations on the following 5 variables.
  \describe{
    \item{Commerce}{a numeric vector of counts.}
    \item{Arts}{a numeric vector of counts.}
    \item{SciEng}{a numeric vector of counts.}
    \item{Law}{a numeric vector of counts.}
    \item{Medicine}{a numeric vector of counts.}
  }
}
\details{
Each student is cross-classified by their colleges (Science
and Engineering have been combined) and the socio-economic
status (SES) of their fathers
(1 = highest, down to 4 = lowest).


}
\source{
    Dr Tony Morrison.


}
\references{

  Wild, C. J. and Seber, G. A. F. (2000)
  \emph{Chance Encounters: A First Course in Data Analysis and Inference},
  New York: Wiley.


}
\examples{
auuc
\dontrun{
round(fitted(grc(auuc)))
round(fitted(grc(auuc, Rank = 2)))
}
}
\keyword{datasets}
