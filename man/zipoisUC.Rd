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
  for the zero-inflated and zero-deflated
  Poisson distribution with parameter \code{pstr0}.

}
\usage{
dzipois(x, lambda, pstr0 = 0, log = FALSE)
pzipois(q, lambda, pstr0 = 0)
qzipois(p, lambda, pstr0 = 0)
rzipois(n, lambda, pstr0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles. }
  \item{p}{vector of probabilities. }
  \item{n}{number of observations. Must be a single positive integer. }
  \item{lambda}{ Vector of positive means. }
  \item{pstr0}{
  Probability of a structural zero
  (i.e., ignoring the Poisson distribution),
  called \eqn{\phi}{phi}.
  The default value of \eqn{\phi = 0}{phi = 0} corresponds
  to the response having an ordinary Poisson distribution.
  This argument may be negative to allow for 0-deflation, hence
  its interpretation as a probability ceases.

  }
  \item{log}{ Logical. Return the logarithm of the answer? }
}
\details{
  The probability function of \eqn{Y} is 0 with probability
  \eqn{\phi}{phi}, and \eqn{Poisson(\lambda)}{Poisson(lambda)} with
  probability \eqn{1-\phi}{1-phi}. Thus
  \deqn{P(Y=0) =\phi + (1-\phi)  P(W=0)}{%
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
\author{ T. W. Yee }
\note{
  The argument \code{pstr0} is recycled to the required length, and
  must have values which lie in the interval \eqn{[0,1]}.


  These functions actually allow for the \emph{zero-deflated
  Poisson} distribution. Here, \code{pstr0} is also permitted
  to lie in the interval \code{[-1/expm1(lambda), 0]}. The
  resulting probability of a zero count is \emph{less than}
  the nominal Poisson value, and the use of \code{pstr0} to
  stand for the probability of a structural zero loses its
  meaning. When \code{pstr0} equals \code{-1/expm1(lambda)}
  this corresponds to the positive-Poisson distribution
  (e.g., see \code{\link{dpospois}}).


}

\seealso{
    \code{\link{zipoisson}},
    \code{\link[stats:Poisson]{dpois}},
    \code{\link{rzinegbin}}.


}
\examples{
lambda <- 3; pstr0 <- 0.2; x <- (-1):7
(ii <- dzipois(x, lambda, pstr0 = pstr0))
max(abs(cumsum(ii) - pzipois(x, lambda, pstr0 = pstr0)))  # Should be 0
table(rzipois(100, lambda, pstr0 = pstr0))

table(qzipois(runif(100), lambda, pstr0))
round(dzipois(0:10, lambda, pstr0 = pstr0) * 100)  # Should be similar

\dontrun{ x <- 0:10
par(mfrow = c(2, 1))  # Zero-inflated Poisson
barplot(rbind(dzipois(x, lambda, pstr0 = pstr0), dpois(x, lambda)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("ZIP(", lambda, ", pstr0 = ", pstr0, ") (blue) vs",
                     " Poisson(", lambda, ") (orange)", sep = ""),
        names.arg = as.character(x))

deflat.limit <- -1 / expm1(lambda)  # Zero-deflated Poisson
newpstr0 <- round(deflat.limit / 1.5, 3)
barplot(rbind(dzipois(x, lambda, pstr0 = newpstr0),
                dpois(x, lambda)),
        beside = TRUE, col = c("blue","orange"),
        main = paste("ZDP(", lambda, ", pstr0 = ", newpstr0, ") (blue) vs",
                     " Poisson(", lambda, ") (orange)", sep = ""),
        names.arg = as.character(x)) }
}
\keyword{distribution}
