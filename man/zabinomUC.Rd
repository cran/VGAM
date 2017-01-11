\name{Zabinom}
\alias{Zabinom}
\alias{dzabinom}
\alias{pzabinom}
\alias{qzabinom}
\alias{rzabinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Altered Binomial Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the zero-altered binomial distribution with parameter \code{pobs0}.

}
\usage{
dzabinom(x, size, prob, pobs0 = 0, log = FALSE)
pzabinom(q, size, prob, pobs0 = 0)
qzabinom(p, size, prob, pobs0 = 0)
rzabinom(n, size, prob, pobs0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.
  }
  \item{size, prob, log}{
  Parameters from the ordinary binomial distribution
  (see \code{\link[stats:Binomial]{dbinom}}).

  }
  \item{pobs0}{
  Probability of (an observed) zero, called \eqn{pobs0}.
  The default value of \code{pobs0 = 0} corresponds
  to the response having a positive binomial distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 0 with probability \code{pobs0},
  else a positive
  binomial(size, prob)
  distribution.


}
\value{
  \code{dzabinom} gives the density and
  \code{pzabinom} gives the distribution function,
  \code{qzabinom} gives the quantile function, and
  \code{rzabinom} generates random deviates.


}
%\references{ }
\author{ T. W. Yee }
\note{
    The argument \code{pobs0} is recycled to the required length, and
    must have values which lie in the interval \eqn{[0,1]}.


}

\seealso{
    \code{\link{zibinomial}},
    \code{\link{rposbinom}}.


%   \code{\link{zabinomial}},


}
\examples{
size <- 10; prob <- 0.15; pobs0 <- 0.05; x <- (-1):7
dzabinom(x, size = size, prob = prob, pobs0 = pobs0)
table(rzabinom(100, size = size, prob = prob, pobs0 = pobs0))

\dontrun{ x <- 0:10
barplot(rbind(dzabinom(x, size = size, prob = prob, pobs0 = pobs0),
                dbinom(x, size = size, prob = prob)),
        beside = TRUE, col = c("blue", "orange"), cex.main = 0.7, las = 1,
        ylab = "Probability", names.arg = as.character(x),
        main = paste("ZAB(size = ", size, ", prob = ", prob, ", pobs0 = ", pobs0,
                   ") [blue] vs",  " Binom(size = ", size, ", prob = ", prob,
                   ") [orange] densities", sep = "")) }
}
\keyword{distribution}
