\name{Biamhcop}
\alias{Biamhcop}
\alias{dbiamhcop}
\alias{pbiamhcop}
\alias{rbiamhcop}
\title{Ali-Mikhail-Haq Bivariate Distribution}
\description{
  Density, distribution function, and random
  generation for the (one parameter) bivariate
  Ali-Mikhail-Haq distribution.


}
\usage{
dbiamhcop(x1, x2, apar, log = FALSE)
pbiamhcop(q1, q2, apar)
rbiamhcop(n, apar)
}
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.}
  \item{n}{number of observations.
  Same as \code{\link[stats]{runif}}


  }
  \item{apar}{the association parameter.}
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm is returned.


  }
}
\value{
  \code{dbiamhcop} gives the density,
  \code{pbiamhcop} gives the distribution function, and
  \code{rbiamhcop} generates random deviates (a two-column matrix).


}
%\references{
%
%}
\author{ T. W. Yee and C. S. Chee}
\details{
  See \code{\link{biamhcop}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation, for the formula of the
  cumulative distribution function and other details.


}
%\note{
%}
\seealso{
  \code{\link{biamhcop}}.


}
\examples{ x <- seq(0, 1, len = (N <- 101)); apar <- 0.7
ox <- expand.grid(x, x)
zedd <- dbiamhcop(ox[, 1], ox[, 2], apar = apar)
\dontrun{
contour(x, x, matrix(zedd, N, N), col = "blue")
zedd <- pbiamhcop(ox[, 1], ox[, 2], apar = apar)
contour(x, x, matrix(zedd, N, N), col = "blue")

plot(r <- rbiamhcop(n = 1000, apar = apar), col = "blue")
par(mfrow = c(1, 2))
hist(r[, 1])  # Should be uniform
hist(r[, 2])  # Should be uniform
}
}
\keyword{distribution}



