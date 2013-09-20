\name{explog}
\alias{explog}
\alias{dexplog}
\alias{pexplog}
\alias{qexplog}
\alias{rexplog}
\title{The Exponential Logarithmic Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the exponential logarithmic distribution.

}
\usage{
dexplog(x, scale = 1, shape, log = FALSE)
pexplog(q, scale = 1, shape)
qexplog(p, scale = 1, shape) 
rexplog(n, scale = 1, shape) 
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    If \code{length(n) > 1} then the length is taken to be the number
    required. }
  \item{scale, shape}{
   positive scale and shape parameters. }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }

}
\value{
  \code{dexplog} gives the density,
  \code{pexplog} gives the distribution function,
  \code{qexplog} gives the quantile function, and
  \code{rexplog} generates random deviates.


}
\author{ J. G. Lauder and T. W. Yee }
\details{
  See \code{\link{explogff}}, the \pkg{VGAM} family function
  for estimating the parameters, 
  for the formula of the probability density function and other details.


}
\note{
  We define \code{scale} as the reciprocal of the scale parameter
  used by Tahmasabi and Rezaei (2008).


}
\seealso{
  \code{\link{explogff}},
  \code{\link{exponential}}.  


}
\examples{
\dontrun{
shape <- 0.5; scale <- 2; nn <- 501
x <- seq(-0.50, 6.0, len = nn)
plot(x, dexplog(x, scale, shape), type = "l", las = 1, ylim = c(0, 1.1),
     ylab = paste("[dp]explog(shape = ", shape, ", scale = ", scale, ")"),
     col = "blue", cex.main = 0.8,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles")
lines(x, pexplog(x, scale, shape), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qexplog(probs, scale, shape = shape)
lines(Q, dexplog(Q, scale, shape = shape), col = "purple", lty = 3, type = "h")
lines(Q, pexplog(Q, scale, shape = shape), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(pexplog(Q, scale, shape = shape) - probs)) # Should be 0
}
}
\keyword{distribution}


