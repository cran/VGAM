\name{gew}
\alias{gew}
\docType{data}
\title{ General Electric and Westinghouse Data }
\description{
  General Electric and Westinghouse capital data.
}
\usage{data(gew)}
\format{
  A data frame with 20 observations on the following 6 variables.
  \describe{
    \item{y1}{a numeric vector which may be regarded as investment
              figures for the two companies}
    \item{x1}{market values}
    \item{x2}{capital stocks}
    \item{y2}{a numeric vector which may be regarded as investment
              figures for the two companies}
    \item{x3}{market values}
    \item{x4}{capital stocks}
  }
}
\details{
  The period is 1934 to 1953.
}
\source{
  Unknown.
}
\references{
Zellner, A. (1962)
An efficient method of estimating seemingly unrelated regressions
and tests for aggregation bias.
\emph{Journal of the American Statistical Association},
\bold{57}, 348--368.

}
\examples{
data(gew)
str(gew)
}
\keyword{datasets}
