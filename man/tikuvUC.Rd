\name{Tikuv}
\alias{Tikuv}
\alias{dtikuv}
\alias{ptikuv}
\alias{qtikuv}
\alias{rtikuv}
\title{A Short-tailed Symmetric Distribution }
\description{
  Density, cumulative distribution function, quantile function and
  random generation for
  the short-tailed symmetric distribution of Tiku and Vaughan (1999).


}
\usage{
dtikuv(x, d, mean = 0, sigma = 1, log = FALSE)
ptikuv(q, d, mean = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)
qtikuv(p, d, mean = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE, ...)
rtikuv(n, d, mean = 0, sigma = 1, Smallno = 1.0e-6)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
           Same as in \code{\link[stats]{runif}}.


  }
  \item{d, mean, sigma }{
  arguments for the parameters of the distribution.
  See \code{\link{tikuv}} for more details.
  For \code{rtikuv}, arguments \code{mean} and \code{sigma} must be of
  length 1.


  }
  \item{Smallno}{
  Numeric, a small value used by the rejection method for determining
  the lower and upper limits of the distribution.
  That is, \code{ptikuv(L) < Smallno} and \code{ptikuv(U) > 1-Smallno}
  where \code{L} and \code{U} are the lower and upper limits respectively.


  }
  \item{\ldots}{
  Arguments that can be passed into \code{\link[stats]{uniroot}}.


  }
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
  \code{dtikuv} gives the density,
  \code{ptikuv} gives the cumulative distribution function,
  \code{qtikuv} gives the quantile function, and
  \code{rtikuv} generates random deviates.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{tikuv}} for more details.


}
%\note{
%
%}
\seealso{
  \code{\link{tikuv}}.


}
\examples{
\dontrun{ par(mfrow = c(2, 1))
x <- seq(-5, 5, len = 401)
plot(x, dnorm(x), type = "l", col = "black", ylab = "", las = 1,
     main = "Black is standard normal, others are dtikuv(x, d)")
lines(x, dtikuv(x, d = -10), col = "orange")
lines(x, dtikuv(x, d = -1 ), col = "blue")
lines(x, dtikuv(x, d =  1 ), col = "green")
legend("topleft", col = c("orange","blue","green"), lty = rep(1, len = 3),
       legend = paste("d =", c(-10, -1, 1)))

plot(x, pnorm(x), type = "l", col = "black", ylab = "", las = 1,
     main = "Black is standard normal, others are ptikuv(x, d)")
lines(x, ptikuv(x, d = -10), col = "orange")
lines(x, ptikuv(x, d = -1 ), col = "blue")
lines(x, ptikuv(x, d =  1 ), col = "green")
legend("topleft", col = c("orange","blue","green"), lty = rep(1, len = 3),
       legend = paste("d =", c(-10, -1, 1))) }

probs <- seq(0.1, 0.9, by = 0.1)
ptikuv(qtikuv(p = probs, d =  1), d = 1) - probs  # Should be all 0
}
\keyword{distribution}


