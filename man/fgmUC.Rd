\name{Fgm}
\alias{Fgm}
\alias{dfgm}
\alias{pfgm}
\alias{rfgm}
\title{Farlie-Gumbel-Morgenstern's Bivariate Distribution}
\description{
  Density, distribution function, and random
  generation for the (one parameter) bivariate 
  Farlie-Gumbel-Morgenstern's distribution.


}
\usage{
dfgm(x1, x2, alpha, log = FALSE)
pfgm(q1, q2, alpha)
rfgm(n, alpha)
}
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{alpha}{the association parameter.}
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm is returned.

  }
}
\value{
  \code{dfgm} gives the density,
  \code{pfgm} gives the distribution function, and
  \code{rfgm} generates random deviates (a two-column matrix).

}
%\references{
%
%}
\author{ T. W. Yee }
\details{
  See \code{\link{fgm}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation, for the formula of the
  cumulative distribution function and other details.


}
%\note{
%}
\seealso{
  \code{\link{fgm}}.


}
\examples{
\dontrun{ N <- 101; x <- seq(0.0, 1.0, len = N); alpha <- 0.7
ox <- expand.grid(x, x)
zedd <- dfgm(ox[, 1], ox[, 2], alpha = alpha)
contour(x, x, matrix(zedd, N, N), col = "blue")
zedd <- pfgm(ox[, 1], ox[, 2], alpha = alpha)
contour(x, x, matrix(zedd, N, N), col = "blue")

plot(r <- rfgm(n = 3000, alpha = alpha), col = "blue")
par(mfrow = c(1, 2))
hist(r[, 1])  # Should be uniform
hist(r[, 2])  # Should be uniform
}
}
\keyword{distribution}


