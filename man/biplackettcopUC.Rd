\name{Biplackett}
\alias{Biplackett}
\alias{dbiplackcop}
\alias{pbiplackcop}
\alias{rbiplackcop}
\title{Plackett's Bivariate Copula }
\description{
  Density, distribution function, and random
  generation for the (one parameter) bivariate Plackett copula.  %distribution.

}
\usage{
dbiplackcop(x1, x2, oratio, log = FALSE)
pbiplackcop(q1, q2, oratio)
rbiplackcop(n, oratio)
}
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.}
  \item{n}{number of observations.
  Same as in \code{\link[stats]{runif}}.


  }
  \item{oratio}{the positive odds ratio \eqn{\psi}{psi}.}
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm is returned.


  }
}
\value{
  \code{dbiplackcop} gives the density,
  \code{pbiplackcop} gives the distribution function, and
  \code{rbiplackcop} generates random deviates (a two-column matrix).


}
\references{


Mardia, K. V. (1967)
Some contributions to contingency-type distributions.
\emph{Biometrika},
\bold{54}, 235--249.


}
\author{ T. W. Yee }
\details{
  See \code{\link{biplackettcop}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation, for the formula of the
  cumulative distribution function and other details.


}
%\note{
%}
\seealso{
  \code{\link{biplackettcop}},
  \code{\link{bifrankcop}}.


}
\examples{
\dontrun{ N <- 101; oratio <- exp(1)
x <- seq(0.0, 1.0, len = N)
ox <- expand.grid(x, x)
zedd <- dbiplackcop(ox[, 1], ox[, 2], oratio = oratio)
contour(x, x, matrix(zedd, N, N), col = "blue")
zedd <- pbiplackcop(ox[, 1], ox[, 2], oratio = oratio)
contour(x, x, matrix(zedd, N, N), col = "blue")

plot(rr <- rbiplackcop(n = 3000, oratio = oratio))
par(mfrow = c(1, 2))
hist(rr[, 1])  # Should be uniform
hist(rr[, 2])  # Should be uniform
}
}
\keyword{distribution}


