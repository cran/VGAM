\name{Zanegbin}
\alias{Zanegbin}
\alias{dzanegbin}
\alias{pzanegbin}
\alias{qzanegbin}
\alias{rzanegbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Altered Negative Binomial Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the zero-altered negative binomial distribution with parameter \code{p0}.

}
\usage{
dzanegbin(x, p0, size, prob = NULL, munb = NULL, log = FALSE)
pzanegbin(q, p0, size, prob = NULL, munb = NULL)
qzanegbin(p, p0, size, prob = NULL, munb = NULL)
rzanegbin(n, p0, size, prob = NULL, munb = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. 
  If \code{length(n) > 1} then the length is taken to be the number required. }
  \item{size, prob, munb, log}{
  Parameters from the ordinary negative binomial distribution
  (see \code{\link[stats:NegBinomial]{dnbinom}}).
  Some arguments have been renamed slightly.

  }
  \item{p0}{
  Probability of zero, called \eqn{p0}.
  The default value of \code{p0 = 0} corresponds
  to the response having a positive negative binomial distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 0 with probability \code{p0},
  else a positive
  negative binomial(\eqn{\mu_{nb}}{munb}, size)
  distribution.

}
\value{
  \code{dzanegbin} gives the density and
  \code{pzanegbin} gives the distribution function,
  \code{qzanegbin} gives the quantile function, and
  \code{rzanegbin} generates random deviates.

}
%\references{ }
\author{ Thomas W. Yee }
\note{ 
    The argument \code{p0} is recycled to the required length, and
    must have values which lie in the interval \eqn{[0,1]}.

}

\seealso{ 
    \code{\link{zanegbinomial}},
    \code{\link{rposnegbin}}.

}
\examples{
munb <- 3; size <- 4; p0 <- 0.3; x <- (-1):7
dzanegbin(x,   p0 = p0, munb = munb, size = size)
table(rzanegbin(100, p0 = p0, munb = munb, size = size))

\dontrun{ x = 0:10
barplot(rbind(dzanegbin(x, p0 = p0, munb = munb, size = size),
                dnbinom(x, mu = munb, size = size)),
        beside = TRUE, col = c("blue","green"), cex.main = 0.7, las = 1,
        ylab = "Probability",names.arg = as.character(x),
        main = paste("ZANB(p0 =", p0, ", munb =", munb, ", size =", size,
                   ") [blue] vs",  " NB(mu =", munb, ", size =", size,
                   ") [green] densities", sep="")) }
}
\keyword{distribution}
