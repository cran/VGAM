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
pexppois(q, rate = 1, shape)
qexppois(p, rate = 1, shape)  
rexppois(n, rate = 1, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    If \code{length(n) > 1} then the length is taken to be the number required. }
  \item{shape, rate}{ both positive parameters. }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }

}
\value{
  \code{dexppois} gives the density,
  \code{pexppois} gives the distribution function,
  \code{qexppois} gives the quantile function, and
  \code{rexppois} generates random deviates.


}
\author{ J. G. Lauder, jamesglauder@gmail.com


}
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
\dontrun{
shape <- 2; rate <- 2; nn <- 201
x <- seq(-0.05, 1.05, len = nn)
plot(x, dexppois(x, shape, rate = rate), type = "l", las = 1, ylim = c(0, 5),
     ylab = paste("[dp]exppoisson(shape = ", shape, ", rate = ", rate, ")"),
     col = "blue", cex.main = 0.8,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles")
lines(x, pexppois(x, shape, rate = rate), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qexppois(probs, shape, rate = rate)
lines(Q, dexppois(Q, shape, rate = rate), col = "purple", lty = 3, type = "h")
lines(Q, pexppois(Q, shape, rate = rate), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(pexppois(Q, shape, rate = rate) - probs))  # Should be 0
}
}
\keyword{distribution}


