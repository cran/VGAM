\name{grain.us}
\alias{grain.us}
\docType{data}
\title{Grain Prices Data in USA }
\description{
A 4-column matrix.

}
\usage{data(grain.us)}
\format{
  The columns are:
  \describe{
    \item{wheat.flour}{numeric}
    \item{corn}{numeric}
    \item{wheat}{numeric}
    \item{rye}{numeric}
  }
}
\details{
  Monthly averages of grain prices
  in the United States for wheat flour, corn, wheat, and rye for the
  period January 1961 through October 1972.
  The units are US dollars per 100 pound sack for wheat flour, and
  per bushel for corn, wheat and rye.


}
\source{
Ahn and Reinsel (1988).


}
\references{

Ahn, S. K  and Reinsel, G. C. (1988)
Nested reduced-rank autoregressive models for multiple time series.
\emph{Journal of the American Statistical Association},
\bold{83}, 849--856.

}
\examples{
\dontrun{
cgrain <- scale(grain.us, scale = FALSE)  # Center the time series only
fit <- vglm(cgrain ~ 1, rrar(Rank = c(4, 1)),
            epsilon = 1e-3, stepsize = 0.5, trace = TRUE, maxit = 50)
summary(fit)
}
}
\keyword{datasets}

