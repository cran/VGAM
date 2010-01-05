\name{Zapois}
\alias{Zapois}
\alias{dzapois}
\alias{pzapois}
\alias{qzapois}
\alias{rzapois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Altered Poisson Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the zero-altered Poisson distribution with parameter \code{p0}.

}
\usage{
dzapois(x, lambda, p0 = 0, log = FALSE)
pzapois(q, lambda, p0 = 0)
qzapois(p, lambda, p0 = 0)
rzapois(n, lambda, p0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required. }
  \item{lambda}{ Vector of positive means. }
  \item{p0}{
  Probability of zero, called \eqn{p0}.
  The default value of \code{p0=0} corresponds
  to the response having a positive Poisson distribution.

  }
  \item{log}{ Logical. Return the logarithm of the answer? }
}
\details{
  The probability function of \eqn{Y} is 0 with probability
  \code{p0}, else a positive \eqn{Poisson(\lambda)}{Poisson(lambda)}.

}
\value{
  \code{dzapois} gives the density,
  \code{pzapois} gives the distribution function,
  \code{qzapois} gives the quantile function, and
  \code{rzapois} generates random deviates.
}
%\references{ }
\author{ Thomas W. Yee }
\note{ 
    The argument \code{p0} is recycled to the required length, and
    must have values which lie in the interval [0,1].

}

\seealso{ 
    \code{\link{zapoisson}}.

}
\examples{
lambda = 3; p0 = 0.2; x = (-1):7
(ii = dzapois(x, lambda, p0))
max(abs(cumsum(ii) - pzapois(x, lambda, p0)))  # Should be 0
table(rzapois(100, lambda, p0))
table(qzapois(runif(100), lambda, p0))
round(dzapois(0:10, lambda, p0) * 100) # Should be similar

\dontrun{
x = 0:10
barplot(rbind(dzapois(x, lambda, p0), dpois(x, lambda)),
        beside = TRUE, col = c("blue","green"),
        main=paste("ZAP(", lambda, ", p0=", p0, ") [blue] vs",
                   " Poisson(", lambda, ") [green] densities", sep=""),
        names.arg = as.character(x), ylab="Probability") }
}
\keyword{distribution}
