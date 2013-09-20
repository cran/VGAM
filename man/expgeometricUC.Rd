\name{expgeom}
\alias{expgeom}
\alias{dexpgeom}
\alias{pexpgeom}
\alias{qexpgeom}
\alias{rexpgeom}
\title{The Exponential Geometric Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the exponential geometric distribution.


}
\usage{
dexpgeom(x, scale = 1, shape, log = FALSE)
pexpgeom(q, scale = 1, shape)
qexpgeom(p, scale = 1, shape) 
rexpgeom(n, scale = 1, shape) 
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    If \code{length(n) > 1} then the length is taken to be the number required. }
  \item{scale, shape}{
   positive scale and shape parameters. }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dexpgeom} gives the density,
  \code{pexpgeom} gives the distribution function,
  \code{qexpgeom} gives the quantile function, and
  \code{rexpgeom} generates random deviates.

}
\author{ J. G. Lauder and T. W. Yee }
\details{
  See \code{\link{expgeometric}}, the \pkg{VGAM} family function
  for estimating the parameters, 
  for the formula of the probability density function and other details.


}
\note{
  We define \code{scale} as the reciprocal of the scale parameter
  used by Adamidis and Loukas (1998).


}
\seealso{
  \code{\link{expgeometric}},
  \code{\link{exponential}}, 
  \code{\link{geometric}}.


}
\examples{
\dontrun{
shape <- 0.5; scale <- 1; nn <- 501
x <- seq(-0.10, 3.0, len = nn)
plot(x, dexpgeom(x, scale, shape), type = "l", las = 1, ylim = c(0, 2),
     ylab = paste("[dp]expgeom(shape = ", shape, ", scale = ", scale, ")"),
     col = "blue", cex.main = 0.8,
     main = "Blue is density, red is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles")
lines(x, pexpgeom(x, scale, shape), col = "red")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qexpgeom(probs, scale, shape)
lines(Q, dexpgeom(Q, scale, shape), col = "purple", lty = 3, type = "h")
lines(Q, pexpgeom(Q, scale, shape), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(pexpgeom(Q, scale, shape) - probs))  # Should be 0
}
}
\keyword{distribution}


