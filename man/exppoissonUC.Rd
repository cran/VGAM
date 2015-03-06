\name{exppois}
\alias{exppois}
\alias{dexppois}
\alias{pexppois}
\alias{qexppois}
\alias{rexppois}
\title{The Exponential Poisson Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the exponential poisson distribution.

}
\usage{
dexppois(x, rate = 1, shape, log = FALSE)
pexppois(q, rate = 1, shape, lower.tail = TRUE, log.p = FALSE)
qexppois(p, rate = 1, shape, lower.tail = TRUE, log.p = FALSE)
rexppois(n, rate = 1, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.


}
  \item{shape, rate}{ positive parameters. }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }
}
\value{
  \code{dexppois} gives the density,
  \code{pexppois} gives the distribution function,
  \code{qexppois} gives the quantile function, and
  \code{rexppois} generates random deviates.


}
\author{ Kai Huang and J. G. Lauder }
\details{
  See \code{\link{exppoisson}}, the \pkg{VGAM} family function
  for estimating the parameters, 
  for the formula of the probability density function and other details.


}
%\note{
%}
\seealso{
  \code{\link{exppoisson}}.


}
\examples{
\dontrun{ rate <- 2; shape <- 0.5; nn <- 201
x <- seq(-0.05, 1.05, len = nn)
plot(x, dexppois(x, rate = rate, shape), type = "l", las = 1, ylim = c(0, 3),
     ylab = paste("fexppoisson(rate = ", rate, ", shape = ", shape, ")"),
     col = "blue", cex.main = 0.8,
     main = "Blue is the density, orange the cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles")
lines(x, pexppois(x, rate = rate, shape), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qexppois(probs, rate = rate, shape)
lines(Q, dexppois(Q, rate = rate, shape), col = "purple", lty = 3, type = "h")
lines(Q, pexppois(Q, rate = rate, shape), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3); abline(h = 0, col = "gray50")
max(abs(pexppois(Q, rate = rate, shape) - probs))  # Should be 0
}
}
\keyword{distribution}


