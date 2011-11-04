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
dzipois(x, lambda, phi = 0, log = FALSE)
pzipois(q, lambda, phi = 0)
qzipois(p, lambda, phi = 0)
rzipois(n, lambda, phi = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles. }
  \item{p}{vector of probabilities. }
  \item{n}{number of observations. Must be a single positive integer. }
  \item{lambda}{ Vector of positive means. }
  \item{phi}{ 
  Probability of zero (ignoring the Poisson distribution), called
  \eqn{\phi}{phi}.  The default value of \eqn{\phi = 0}{phi = 0} corresponds
  to the response having an ordinary Poisson distribution.

  }
  \item{log}{ Logical. Return the logarithm of the answer? }
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
    must have values which lie in the interval \eqn{[0,1]}.

}

\seealso{ 
    \code{\link{zipoisson}},
    \code{\link[stats:Poisson]{dpois}},
    \code{\link{rzinegbin}}.


}
\examples{
lambda = 3; phi = 0.2; x = (-1):7
(ii = dzipois(x, lambda, phi))
max(abs(cumsum(ii) - pzipois(x, lambda, phi))) # Should be 0
table(rzipois(100, lambda, phi))

table(qzipois(runif(100), lambda, phi))
round(dzipois(0:10, lambda, phi) * 100) # Should be similar

\dontrun{ x = 0:10
barplot(rbind(dzipois(x, lambda, phi), dpois(x, lambda)),
        beside = TRUE, col = c("blue","orange"),
        main = paste("ZIP(", lambda, ", phi = ", phi, ") (blue) vs",
                     " Poisson(", lambda, ") (orange)", sep = ""),
        names.arg = as.character(x)) }
}
\keyword{distribution}
