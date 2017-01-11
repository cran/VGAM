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
  for the zero-altered negative binomial distribution with parameter \code{pobs0}.

}
\usage{
dzanegbin(x, size, prob = NULL, munb = NULL, pobs0 = 0, log = FALSE)
pzanegbin(q, size, prob = NULL, munb = NULL, pobs0 = 0)
qzanegbin(p, size, prob = NULL, munb = NULL, pobs0 = 0)
rzanegbin(n, size, prob = NULL, munb = NULL, pobs0 = 0)
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
  \item{pobs0}{
  Probability of zero, called \eqn{pobs0}.
  The default value of \code{pobs0 = 0} corresponds
  to the response having a positive negative binomial distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 0 with probability \code{pobs0},
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
\author{ T. W. Yee }
\note{
    The argument \code{pobs0} is recycled to the required length, and
    must have values which lie in the interval \eqn{[0,1]}.


}

\seealso{
    \code{\link{zanegbinomial}},
    \code{\link{rposnegbin}}.



}
\examples{
munb <- 3; size <- 4; pobs0 <- 0.3; x <- (-1):7
dzanegbin(x, munb = munb, size = size, pobs0 = pobs0)
table(rzanegbin(100, munb = munb, size = size, pobs0 = pobs0))

\dontrun{ x <- 0:10
barplot(rbind(dzanegbin(x, munb = munb, size = size, pobs0 = pobs0),
                dnbinom(x, mu   = munb, size = size)),
        beside = TRUE, col = c("blue", "green"), cex.main = 0.7, las = 1,
        ylab = "Probability", names.arg = as.character(x),
        main = paste("ZANB(munb = ", munb, ", size = ", size,",
                     pobs0 = ", pobs0,
                   ") [blue] vs",  " NB(mu = ", munb, ", size = ", size,
                   ") [green] densities", sep = "")) }
}
\keyword{distribution}
