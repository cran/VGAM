\name{Zipois}
\alias{Zipois}
\alias{dzipois}
\alias{pzipois}
\alias{qzipois}
\alias{rzipois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Poisson Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the zero-inflated Poisson distribution with parameter \code{phi}.

}
\usage{
dzipois(x, lambda, phi = 0)
pzipois(q, lambda, phi = 0)
qzipois(p, lambda, phi = 0)
rzipois(n, lambda, phi = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Must be a single positive integer. }
  \item{lambda}{ Vector of positive means. }
  \item{phi}{ 
  Probability of zero (ignoring the Poisson distribution), called
  \eqn{\phi}{phi}.  The default value of \eqn{\phi=0}{phi=0} corresponds
  to the response having an ordinary Poisson distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 0 with probability
  \eqn{\phi}{phi}, and \eqn{Poisson(\lambda)}{Poisson(lambda)} with
  probability \eqn{1-\phi}{1-phi}. Thus 
  \deqn{P(Y=0) =\phi + (1-\phi) P(W=0)}{%
        P(Y=0) = phi + (1-phi) * P(W=0)}
  where \eqn{W} is distributed \eqn{Poisson(\lambda)}{Poisson(lambda)}.
}
\value{
  \code{dzipois} gives the density,
  \code{pzipois} gives the distribution function,
  \code{qzipois} gives the quantile function, and
  \code{rzipois} generates random deviates.
}
%\references{ }
\author{ Thomas W. Yee }
\note{ 
    The argument \code{phi} is recycled to the required length, and
    must have values which lie in the interval [0,1].

}

\seealso{ 
    \code{\link[stats:Poisson]{dpois}},
    \code{\link{zipoisson}}.
}
\examples{
lambda = 4
phi = 0.8
(i = dzipois(0:7, lambda, phi))
cumsum(i) - pzipois(0:7, lambda, phi)  # Should be 0s
table(rzipois(100, lambda, phi))

table(qzipois(runif(100), lambda, phi))
round(dzipois(0:10, lambda, phi) * 100) # Should be similar

\dontrun{
x = 0:10
plot(x, dzipois(x, lambda, phi), type="h", ylab="Probability", 
     main=paste("ZIP(", lambda, ", phi=", phi, ") (blue) vs",
     " Poisson(", lambda, ") (red and shifted slightly)", sep=""),
     lwd=2, col="blue", las=1)
lines(x+0.05, dpois(x, lambda), type="h", lwd=2, col="red")
}
}
\keyword{distribution}
