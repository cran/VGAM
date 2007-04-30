\name{Zibinom}
\alias{Zibinom}
\alias{dzibinom}
\alias{pzibinom}
\alias{qzibinom}
\alias{rzibinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Binomial Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the zero-inflated binomial distribution with parameter \code{phi}.

}
\usage{
dzibinom(x, size, prob, log = FALSE, phi = 0)
pzibinom(q, size, prob, lower.tail = TRUE, log.p = FALSE, phi = 0)
qzibinom(p, size, prob, lower.tail = TRUE, log.p = FALSE, phi = 0)
rzibinom(n, size, prob, phi = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{size}{number of trials. It is the \eqn{N} symbol in the formula
   given in \code{\link{zibinomial}}. }
  \item{prob}{probability of success on each trial. }
  \item{n}{number of observations. Must be a single positive integer. }
  \item{log, log.p, lower.tail}{ Arguments that are passed on to 
      \code{\link[stats:Binomial]{pbinom}}.}
  \item{phi}{ 
  Probability of zero (ignoring the binomial distribution), called
  \eqn{\phi}{phi}.  The default value of \eqn{\phi=0}{phi=0} corresponds
  to the response having an ordinary binomial distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 0 with probability
  \eqn{\phi}{phi}, and \eqn{Binomial(size, prob)}{Binomial(size, prob)} with
  probability \eqn{1-\phi}{1-phi}. Thus 
  \deqn{P(Y=0) =\phi + (1-\phi) P(W=0)}{%
        P(Y=0) = phi + (1-phi) * P(W=0)}
  where \eqn{W} is distributed \eqn{Binomial(size, prob)}{Binomial(size, prob)}.
}
\value{
  \code{dzibinom} gives the density,
  \code{pzibinom} gives the distribution function,
  \code{qzibinom} gives the quantile function, and
  \code{rzibinom} generates random deviates.
}
%\references{ }
\author{ Thomas W. Yee }
\note{ 
    The argument \code{phi} is recycled to the required length, and
    must have values which lie in the interval [0,1].

}

\seealso{ 
    \code{\link{zibinomial}},
    \code{\link[stats:Binomial]{dbinom}}.
}
\examples{
prob = 0.2
size = 10
phi = 0.5
(i = dzibinom(0:size, size, prob, phi=phi))
cumsum(i) - pzibinom(0:size, size, prob, phi=phi)  # Should be 0s
table(rzibinom(100, size, prob, phi=phi))

table(qzibinom(runif(100), size, prob, phi=phi))
round(dzibinom(0:10, size, prob, phi=phi) * 100) # Should be similar

\dontrun{
x = 0:size
plot(x, dzibinom(x, size, prob, phi=phi), type="h", ylab="Probability", 
     main=paste("ZIB(", size, ", ", prob, ", phi=", phi, ") (blue) vs",
     " Binomial(", size, ", ", prob, ") (red and shifted slightly)", sep=""),
     lwd=2, col="blue", las=1)
lines(x+0.05, dbinom(x, size, prob), type="h", lwd=2, col="red")
}
}
\keyword{distribution}