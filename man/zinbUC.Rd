\name{Zinb}
\alias{Zinb}
\alias{dzinb}
\alias{pzinb}
\alias{qzinb}
\alias{rzinb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Negative Binomial Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the zero-inflated negative binomial distribution with parameter \code{phi}.

}
\usage{
dzinb(x, phi, size, prob=NULL, munb=NULL, log=FALSE)
pzinb(q, phi, size, prob=NULL, munb=NULL)
qzinb(p, phi, size, prob=NULL, munb=NULL)
rzinb(n, phi, size, prob=NULL, munb=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Must be a single positive integer. }
  \item{size, prob, munb, log}{
  Arguments matching \code{\link[stats:NegBinomial]{dnbinom}}.
  The argument \code{munb} corresponds to \code{mu} in
  \code{\link[stats:NegBinomial]{dnbinom}} and has been renamed
  to emphasize the fact that it is the mean of the negative binomial
  \emph{component}.

  }
  \item{phi}{ 
  Probability of zero (ignoring the negative binomial distribution), called
  \eqn{\phi}{phi}.

  }
}
\details{
  The probability function of \eqn{Y} is 0 with probability
  \eqn{\phi}{phi}, and a negative binomial distribution with
  probability \eqn{1-\phi}{1-phi}. Thus 
  \deqn{P(Y=0) =\phi + (1-\phi) P(W=0)}{%
        P(Y=0) = phi + (1-phi) * P(W=0)}
  where \eqn{W} is distributed as a negative binomial distribution
  (see \code{\link[stats:NegBinomial]{rnbinom}}.)
  See \code{\link{negbinomial}}, a \pkg{VGAM} family function,
  for the formula of the probability density function and other details
  of the negative binomial distribution.

}
\value{
  \code{dzinb} gives the density,
  \code{pzinb} gives the distribution function,
  \code{qzinb} gives the quantile function, and
  \code{rzinb} generates random deviates.
}
%\references{ }
\author{ Thomas W. Yee }
\note{ 
    The argument \code{phi} is recycled to the required length, and
    must have values which lie in the interval [0,1].

}

\seealso{ 
    \code{\link{zinegbinomial}},
    \code{\link[stats:NegBinomial]{rnbinom}},
    \code{\link{rzipois}}.

}
\examples{
munb = 3; phi = 0.2; k = 10
x = 0:10
(i = dzinb(x, phi=phi, mu=munb, size=k))
cumsum(i) - pzinb(x, phi=phi, mu=munb, size=k)  # Should be 0s
table(rzinb(100, phi=phi, mu=munb, size=k))

table(qzinb(runif(1000), phi=phi, mu=munb, size=k))
round(dzinb(x, phi=phi, mu=munb, size=k) * 1000) # Should be similar

\dontrun{
barplot(rbind(dzinb(x, phi=phi, mu=munb, size=k),
              dnbinom(x, mu=munb, size=k)),
        beside = TRUE, col = c("blue","green"),
        main=paste("ZINB(mu=", munb, ", k=", k, ", phi=", phi,
                   ") (blue) vs negative binomial(mu=", munb,
                   ", k=", k, ") (green)", sep=""),
        names.arg = as.character(x))
}
}
\keyword{distribution}
