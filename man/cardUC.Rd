\name{Card}
\alias{Card}
\alias{dcard}
\alias{pcard}
\alias{qcard}
\alias{rcard}
\alias{dcard2}
\alias{pcard2}
\alias{qcard2}
\alias{rcard2}
\title{ Cardioid Distribution }
\description{
  Density, distribution function, quantile function
  and random
  generation for the cardioid distribution.

}
\usage{
dcard(x, mu, rho, log = FALSE)
pcard(q, mu, rho, lower.tail = TRUE, log.p = FALSE)
qcard(p, mu, rho, lower.tail = TRUE, log.p = FALSE, tol0 = 1e-6)
rcard(n, mu, rho)
dcard2(x, mu, rho2 = 0, log = FALSE)
pcard2(q, mu, rho2 = 0, lower.tail = TRUE, log.p = FALSE)
qcard2(p, mu, rho2 = 0, lower.tail = TRUE, log.p = FALSE,
       tol0 = 1e-6)
rcard2(n, mu, rho2, tol0 = 1e-4)
}
%- maybe also 'usage' for other objs documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as in \code{\link[stats]{runif}}.


  }
  \item{mu, rho, rho2}{
    See \code{\link{cardioid}} and \code{\link{cardioid2}}
    for more information.
    Note: \code{rho2 = 2 * rho}, or equivalently,
    \eqn{\rho_2 = 2 \rho}{rho2 = 2 * rho}
    so is in \eqn{[0, 1]}.


  }
  \item{tol0}{
    Numeric, positive and small.
    Used under different contexts, any value
    within a certain tolerance is treated as being
    equivalent to 0.
    For example, if \code{0 <= rho2 < tol0}
    then the distribution is treated as uniform.
    Another context is for the 
    algorithm to solve for the roots of a
    nonlinear system of equations.


  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of
  the density is returned.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }

}
\details{
  See \code{\link{cardioid}} and \code{\link{cardioid2}},
  the \pkg{VGAM}
  family functions
  for estimating the two parameters by
  maximum likelihood
  estimation, for the formula of the
  probability density
  function and other details.


}
\value{
  \code{dcard}/\code{dcard2} give the density,
  \code{pcard}/\code{pcard2} give the distribution function,
  \code{qcard}/\code{qcard2} give the quantile function, and
  \code{rcard}/\code{rcard2} generates random deviates.


}
\references{

Pewsey, A. (2025).
On Jeffreys's cardioid distribution.
\emph{Computational Statistics and Data Analysis},
\bold{82}, in press.


%\bold{82}, in press.




}
\author{ Thomas W. Yee}
\note{
  Convergence problems might occur with \code{rcard}.


}

\seealso{
  \code{\link{cardioid}},
  \code{\link{cardioid2}}.


}
\examples{
\dontrun{
mu <- 4; rho2 <- 0.4 * 2; x <- seq(-0.5, 2*pi+0.5, len = 1001)
plot(x, dcard2(x, mu, rho2), type = "l", las = 1, ylim = 0:1,
     ylab = paste0("[dp]card2(mu=", mu, ", rho2=", rho2, ")"),
     main = "Blue is density, orange the CDF", col = "blue",
     sub = "Purple lines are the 10,20,...,90 percentiles")
lines(x, pcard2(x, mu, rho2), col = "orange")

probs <- seq(0, 1, by = 0.1)
Q <- qcard2(probs, mu, rho2)
lines(Q, dcard2(Q, mu, rho2), col = 3, lty = 3, type = "h")
lines(Q, pcard2(Q, mu, rho2), col = 3, lty = 3, type = "h")
abline(h = probs, v = c(0, 2*pi), col = 3, lty = 3)
max(abs(pcard2(Q, mu, rho2) - probs))  # Should be 0
}
}
\keyword{distribution}
