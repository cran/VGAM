\name{Card}
\alias{Card}
\alias{dcard}
\alias{pcard}
\alias{qcard}
\alias{rcard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cardioid Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the cardioid distribution.

}
\usage{
dcard(x, mu, rho, log = FALSE)
pcard(q, mu, rho, lower.tail = TRUE, log.p = FALSE)
qcard(p, mu, rho, tolerance = 1e-07, maxits = 500,
      lower.tail = TRUE, log.p = FALSE)
rcard(n, mu, rho, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as in \code{\link[stats]{runif}}.


  }
  \item{mu, rho}{
  See \code{\link{cardioid}} for more information.


  }
  \item{tolerance, maxits, ...}{
  The first two are control parameters for the algorithm
  used to solve for the roots of a nonlinear system of equations;
  \code{tolerance} controls for the accuracy and
  \code{maxits} is the maximum number of iterations.
  \code{rcard} calls \code{qcard} so the \code{...} can be used
  to vary the two arguments.


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
\details{
  See \code{\link{cardioid}}, the \pkg{VGAM} family function
  for estimating the two parameters by maximum likelihood estimation,
  for the formula of the probability density function and other details.


}
\value{
  \code{dcard} gives the density,
  \code{pcard} gives the distribution function,
  \code{qcard} gives the quantile function, and
  \code{rcard} generates random deviates.


}
%\references{ }
\author{ Thomas W. Yee and Kai Huang }
\note{
  Convergence problems might occur with \code{rcard}.


}

\seealso{
  \code{\link{cardioid}}.


}
\examples{
\dontrun{
mu <- 4; rho <- 0.4; x <- seq(0, 2*pi, len = 501)
plot(x, dcard(x, mu, rho), type = "l", las = 1, ylim = c(0, 1), col = "blue",
     ylab = paste("[dp]card(mu=", mu, ", rho=", rho, ")"),
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles")
lines(x, pcard(x, mu, rho), col = "orange")

probs <- seq(0.1, 0.9, by = 0.1)
Q <- qcard(probs, mu, rho)
lines(Q, dcard(Q, mu, rho), col = "purple", lty = 3, type = "h")
lines(Q, pcard(Q, mu, rho), col = "purple", lty = 3, type = "h")
abline(h = c(0,probs, 1), v = c(0, 2*pi), col = "purple", lty = 3)
max(abs(pcard(Q, mu, rho) - probs))  # Should be 0
}
}
\keyword{distribution}
