\name{Zageom}
\alias{Zageom}
\alias{dzageom}
\alias{pzageom}
\alias{qzageom}
\alias{rzageom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Altered Geometric Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the zero-altered geometric distribution with parameter \code{pobs0}.

}
\usage{
dzageom(x, prob, pobs0 = 0, log = FALSE)
pzageom(q, prob, pobs0 = 0)
qzageom(p, prob, pobs0 = 0)
rzageom(n, prob, pobs0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.
  }
  \item{prob, log}{
  Parameters from the ordinary geometric distribution
  (see \code{\link[stats:Geometric]{dgeom}}).

  }
  \item{pobs0}{
  Probability of (an observed) zero, called \eqn{pobs0}.
  The default value of \code{pobs0 = 0} corresponds
  to the response having a positive geometric distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 0 with probability \code{pobs0},
  else a positive
  geometric(prob)
  distribution.

}
\value{
  \code{dzageom} gives the density and
  \code{pzageom} gives the distribution function,
  \code{qzageom} gives the quantile function, and
  \code{rzageom} generates random deviates.

}
%\references{ }
\author{ T. W. Yee }
\note{
    The argument \code{pobs0} is recycled to the required length, and
    must have values which lie in the interval \eqn{[0,1]}.

}

\seealso{
    \code{\link{zageometric}},
    \code{\link{zigeometric}},
    \code{\link{rposgeom}}.


}
\examples{
prob <- 0.35; pobs0 <- 0.05; x <- (-1):7
dzageom(x, prob = prob, pobs0 = pobs0)
table(rzageom(100, prob = prob, pobs0 = pobs0))

\dontrun{ x <- 0:10
barplot(rbind(dzageom(x, prob = prob, pobs0 = pobs0),
                dgeom(x, prob = prob)),
        beside = TRUE, col = c("blue", "orange"), cex.main = 0.7, las = 1,
        ylab = "Probability", names.arg = as.character(x),
        main = paste("ZAG(prob = ", prob, ", pobs0 = ", pobs0,
                   ") [blue] vs",  " Geometric(prob = ", prob,
                   ") [orange] densities", sep = "")) }
}
\keyword{distribution}
