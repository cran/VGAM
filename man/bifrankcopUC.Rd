\name{Frank}
\alias{Frank}
\alias{dbifrankcop}
\alias{pbifrankcop}
\alias{rbifrankcop}
\title{Frank's Bivariate Distribution}
\description{
  Density, distribution function, and random
  generation for the (one parameter) bivariate Frank distribution.

}
\usage{
dbifrankcop(x1, x2, apar, log = FALSE)
pbifrankcop(q1, q2, apar)
rbifrankcop(n, apar)
}
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.}
  \item{n}{number of observations.
    Same as in \code{\link[stats]{runif}}.


  }
  \item{apar}{the positive association parameter. }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }

}
\value{
  \code{dbifrankcop} gives the density,
  \code{pbifrankcop} gives the distribution function, and
  \code{rbifrankcop} generates random deviates (a two-column matrix).


}
\references{

Genest, C. (1987)
Frank's family of bivariate distributions.
\emph{Biometrika},
\bold{74}, 549--555.

}
\author{ T. W. Yee }
\details{
  See \code{\link{bifrankcop}}, the \pkg{VGAM}
  family functions for estimating the association
  parameter by maximum likelihood estimation, for the formula of the
  cumulative distribution function and other details.


}
%\note{
%}
\seealso{
  \code{\link{bifrankcop}}.


}
\examples{
\dontrun{N <- 100; apar <- exp(2)
xx <- seq(-0.30, 1.30, len = N)
ox <- expand.grid(xx, xx)
zedd <- dbifrankcop(ox[, 1], ox[, 2], apar = apar)
contour(xx, xx, matrix(zedd, N, N))
zedd <- pbifrankcop(ox[, 1], ox[, 2], apar = apar)
contour(xx, xx, matrix(zedd, N, N))

plot(rr <- rbifrankcop(n = 3000, apar = exp(4)))
par(mfrow = c(1, 2))
hist(rr[, 1]); hist(rr[, 2])  # Should be uniform
}
}
\keyword{distribution}


