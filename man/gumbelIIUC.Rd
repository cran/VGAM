\name{Gumbel-II}
\alias{Gumbel-II}
\alias{dgumbelII}
\alias{pgumbelII}
\alias{qgumbelII}
\alias{rgumbelII}
\title{The Gumbel-II Distribution}
\description{
  Density, cumulative distribution function,
  quantile function
  and
  random generation for
  the Gumbel-II distribution.

}
\usage{
dgumbelII(x, scale = 1, shape, log = FALSE)
pgumbelII(q, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
qgumbelII(p, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
rgumbelII(n, scale = 1, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as in \code{\link[stats]{runif}}.


  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }
  \item{shape, scale}{positive shape and scale parameters. }

}
\value{
  \code{dgumbelII} gives the density,
  \code{pgumbelII} gives the cumulative distribution function,
  \code{qgumbelII} gives the quantile function, and
  \code{rgumbelII} generates random deviates.


}
\author{ T. W. Yee and Kai Huang
}

\details{
  See \code{\link{gumbelII}} for details.

}
%\note{
%
%}
\seealso{
  \code{\link{gumbelII}},
  \code{\link{dgumbel}}.


}
\examples{
probs <- seq(0.01, 0.99, by = 0.01)
Scale <- exp(1); Shape <- exp( 0.5);
max(abs(pgumbelII(qgumbelII(p = probs, shape = Shape, Scale),
                  shape = Shape, Scale) - probs))  # Should be 0

\dontrun{ x <- seq(-0.1, 10, by = 0.01);
plot(x, dgumbelII(x, shape = Shape, Scale), type = "l", col = "blue", las = 1,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles",
     ylab = "", ylim = 0:1)
abline(h = 0, col = "blue", lty = 2)
lines(x, pgumbelII(x, shape = Shape, Scale), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qgumbelII(probs, shape = Shape, Scale)
lines(Q, dgumbelII(Q, Scale, Shape), col = "purple", lty = 3, type = "h")
pgumbelII(Q, shape = Shape, Scale) - probs # Should be all zero
abline(h = probs, col = "purple", lty = 3) }
}
\keyword{distribution}


