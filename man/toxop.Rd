\name{toxop}
\alias{toxop}
\docType{data}
\title{ Toxoplasmosis Data }
\description{
  Toxoplasmosis data in 34 cities in El Salvador.

}
\usage{data(toxop)}
\format{
  A data frame with 34 observations on the following 4 variables.
  \describe{
    \item{\code{rainfall}}{
  a numeric vector; the amount of rainfall
  in each city.

  }
    \item{\code{ssize}}{a numeric vector; sample size.}
    \item{\code{cityNo}}{a numeric vector; the city number.}
    \item{\code{positive}}{a numeric vector; the number of subjects
          testing positive for the disease. }
  }
}
\details{
  See the references for details.



}
\source{
  See the references for details.



}

\seealso{
  \code{\link[VGAM]{double.expbinomial}}.


}

\references{
  Efron, B. (1978)
  Regression and ANOVA With zero-one data: measures of
  residual variation.
  \emph{Journal of the American Statistical Association},
  \bold{73}, 113--121.


  Efron, B. (1986)
  Double exponential families and their use in generalized linear regression.
  \emph{Journal of the American Statistical Association},
  \bold{81}, 709--721.


}
\examples{
\dontrun{ with(toxop, plot(rainfall, positive / ssize, col = "blue"))
plot(toxop, col = "blue") }
}
\keyword{datasets}
