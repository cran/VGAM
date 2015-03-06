\name{Gompertz}
\alias{Gompertz}
\alias{dgompertz}
\alias{pgompertz}
\alias{qgompertz}
\alias{rgompertz}
\title{The Gompertz Distribution}
\description{
  Density, cumulative distribution function,
  quantile function
  and
  random generation for
  the Gompertz distribution.


}
\usage{
dgompertz(x, scale = 1, shape, log = FALSE)
pgompertz(q, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
qgompertz(p, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
rgompertz(n, scale = 1, shape)
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
  \item{scale, shape}{positive scale and shape parameters. }

}
\value{
  \code{dgompertz} gives the density,
  \code{pgompertz} gives the cumulative distribution function,
  \code{qgompertz} gives the quantile function, and
  \code{rgompertz} generates random deviates.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{gompertz}} for details.


}
%\note{
%
%}
\seealso{
  \code{\link{gompertz}},
  \code{\link{dgumbel}},
  \code{\link{dmakeham}}.


}
\examples{
probs <- seq(0.01, 0.99, by = 0.01)
Shape <- exp(1); Scale <- exp(1)
max(abs(pgompertz(qgompertz(p = probs, Scale, shape = Shape),
                  Scale, shape = Shape) - probs))  # Should be 0

\dontrun{ x <- seq(-0.1, 1.0, by = 0.001)
plot(x, dgompertz(x, Scale,shape = Shape), type = "l", col = "blue", las = 1,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles",
     ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(x, pgompertz(x, Scale, shape = Shape), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qgompertz(probs, Scale, shape = Shape)
lines(Q, dgompertz(Q, Scale, shape = Shape), col = "purple",
      lty = 3, type = "h")
pgompertz(Q, Scale, shape = Shape) - probs  # Should be all zero
abline(h = probs, col = "purple", lty = 3) }
}
\keyword{distribution}


