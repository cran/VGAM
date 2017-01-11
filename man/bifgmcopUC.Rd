\name{Bifgmcop}
\alias{Bifgmcop}
\alias{dbifgmcop}
\alias{pbifgmcop}
\alias{rbifgmcop}
\title{Farlie-Gumbel-Morgenstern's Bivariate Distribution}
\description{
  Density, distribution function, and random
  generation for the (one parameter) bivariate
  Farlie-Gumbel-Morgenstern's distribution.


}
\usage{
dbifgmcop(x1, x2, apar, log = FALSE)
pbifgmcop(q1, q2, apar)
rbifgmcop(n, apar)
}
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.}
  \item{n}{number of observations.
    Same as in \code{\link[stats]{runif}}.



  }
  \item{apar}{the association parameter.}
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm is returned.

  }
}
\value{
  \code{dbifgmcop} gives the density,
  \code{pbifgmcop} gives the distribution function, and
  \code{rbifgmcop} generates random deviates (a two-column matrix).


}
%\references{
%
%}
\author{ T. W. Yee }
\details{
  See \code{\link{bifgmcop}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation, for the formula of the
  cumulative distribution function and other details.


}
%\note{
%}
\seealso{
  \code{\link{bifgmcop}}.


}
\examples{
\dontrun{ N <- 101; x <- seq(0.0, 1.0, len = N); apar <- 0.7
ox <- expand.grid(x, x)
zedd <- dbifgmcop(ox[, 1], ox[, 2], apar = apar)
contour(x, x, matrix(zedd, N, N), col = "blue")
zedd <- pbifgmcop(ox[, 1], ox[, 2], apar = apar)
contour(x, x, matrix(zedd, N, N), col = "blue")

plot(r <- rbifgmcop(n = 3000, apar = apar), col = "blue")
par(mfrow = c(1, 2))
hist(r[, 1])  # Should be uniform
hist(r[, 2])  # Should be uniform
}
}
\keyword{distribution}


