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
dgumbelII(x, shape, scale = 1, log = FALSE)
pgumbelII(q, shape, scale = 1)
qgumbelII(p, shape, scale = 1)
rgumbelII(n, shape, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }
  \item{shape, scale}{positive shape and scale parameters. }

}
\value{
  \code{dgumbelII} gives the density,
  \code{pgumbelII} gives the cumulative distribution function,
  \code{qgumbelII} gives the quantile function, and
  \code{rgumbelII} generates random deviates.


}
\author{ T. W. Yee }
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
Shape <- exp( 0.5); Scale <- exp(1);
max(abs(pgumbelII(qgumbelII(p = probs, Shape, Scale),
                  Shape, Scale) - probs)) # Should be 0

\dontrun{ x <- seq(-0.1, 10, by = 0.01);
plot(x, dgumbelII(x, Shape, Scale), type = "l", col = "blue", las = 1,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles",
     ylab = "", ylim = 0:1)
abline(h = 0, col = "blue", lty = 2)
lines(x, pgumbelII(x, Shape, Scale), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qgumbelII(probs, Shape, Scale)
lines(Q, dgumbelII(Q, Shape, Scale), col = "purple", lty = 3, type = "h")
pgumbelII(Q, Shape, Scale) - probs # Should be all zero
abline(h = probs, col = "purple", lty = 3) }
}
\keyword{distribution}


