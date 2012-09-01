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
dgompertz(x, shape, scale = 1, log = FALSE)
pgompertz(q, shape, scale = 1)
qgompertz(p, shape, scale = 1)
rgompertz(n, shape, scale = 1)
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
  \code{dgompertz} gives the density,
  \code{pgompertz} gives the cumulative distribution function,
  \code{qgompertz} gives the quantile function, and
  \code{rgompertz} generates random deviates.


}
\author{ T. W. Yee }
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
Shape <- exp(1); Scale <- exp(1);
max(abs(pgompertz(qgompertz(p = probs, Shape, Scale),
                  Shape, Scale) - probs)) # Should be 0

\dontrun{ x <- seq(-0.1, 1.0, by = 0.01);
plot(x, dgompertz(x, Shape, Scale), type = "l", col = "blue", las = 1,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles",
     ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(x, pgompertz(x, Shape, Scale), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qgompertz(probs, Shape, Scale)
lines(Q, dgompertz(Q, Shape, Scale), col = "purple", lty = 3, type = "h")
pgompertz(Q, Shape, Scale) - probs # Should be all zero
abline(h = probs, col = "purple", lty = 3) }
}
\keyword{distribution}


