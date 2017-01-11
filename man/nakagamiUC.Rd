\name{Nakagami}
\alias{Nakagami}
\alias{dnaka}
\alias{pnaka}
\alias{qnaka}
\alias{rnaka}
\title{Nakagami Distribution }
\description{
  Density, cumulative distribution function, quantile function and
  random generation for
  the Nakagami distribution.

}
\usage{
dnaka(x, scale = 1, shape, log = FALSE)
pnaka(q, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
qnaka(p, scale = 1, shape, ...)
rnaka(n, scale = 1, shape, Smallno = 1.0e-6)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Same as in \code{\link[stats]{runif}}.


%   Must be a positive integer of length 1.


  }
  \item{scale, shape}{
  arguments for the parameters of the distribution.
  See \code{\link{nakagami}} for more details.
  For \code{rnaka}, arguments \code{shape} and \code{scale} must be of
  length 1.


  }
  \item{Smallno}{
  Numeric, a small value used by the rejection method for determining
  the upper limit of the distribution.
  That is, \code{pnaka(U) > 1-Smallno} where \code{U} is the upper limit.


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
  \code{dnaka} gives the density,
  \code{pnaka} gives the cumulative distribution function,
  \code{qnaka} gives the quantile function, and
  \code{rnaka} generates random deviates.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{nakagami}} for more details.


}
%\note{
%
%}
\seealso{
  \code{\link{nakagami}}.


}
\examples{
\dontrun{ x <- seq(0, 3.2, len = 200)
plot(x, dgamma(x, shape = 1), type = "n", col = "black", ylab = "",
     ylim = c(0,1.5), main = "dnaka(x, shape = shape)")
lines(x, dnaka(x, shape = 1), col = "orange")
lines(x, dnaka(x, shape = 2), col = "blue")
lines(x, dnaka(x, shape = 3), col = "green")
legend(2, 1.0, col = c("orange","blue","green"), lty = rep(1, len = 3),
       legend = paste("shape =", c(1, 2, 3)))

plot(x, pnorm(x), type = "n", col = "black", ylab = "",
     ylim = 0:1, main = "pnaka(x, shape = shape)")
lines(x, pnaka(x, shape = 1), col = "orange")
lines(x, pnaka(x, shape = 2), col = "blue")
lines(x, pnaka(x, shape = 3), col = "green")
legend(2, 0.6, col = c("orange","blue","green"), lty = rep(1, len = 3),
       legend = paste("shape =", c(1, 2, 3))) }

probs <- seq(0.1, 0.9, by = 0.1)
pnaka(qnaka(p = probs, shape = 2), shape = 2) - probs  # Should be all 0
}
\keyword{distribution}


