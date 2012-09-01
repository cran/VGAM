\name{Amh}
\alias{Amh}
\alias{damh}
\alias{pamh}
\alias{ramh}
\title{Ali-Mikhail-Haq Distribution's Bivariate Distribution}
\description{
  Density, distribution function, and random
  generation for the (one parameter) bivariate
  Ali-Mikhail-Haq distribution.

}
\usage{
damh(x1, x2, alpha, log = FALSE)
pamh(q1, q2, alpha)
ramh(n, alpha)
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
  \code{damh} gives the density,
  \code{pamh} gives the distribution function, and
  \code{ramh} generates random deviates (a two-column matrix).
}
%\references{
%
%}
\author{ T. W. Yee and C. S. Chee}
\details{
  See \code{\link{amh}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation, for the formula of the
  cumulative distribution function and other details.

}
%\note{
%}
\seealso{
  \code{\link{amh}}.

}
\examples{
x = seq(0, 1, len=(N <- 101))
alpha = 0.7
ox = expand.grid(x, x)
z = damh(ox[,1], ox[,2], alpha=alpha)
\dontrun{
contour(x, x, matrix(z, N, N), col="blue")
z = pamh(ox[,1], ox[,2], alpha=alpha)
contour(x, x, matrix(z, N, N), col="blue")

plot(r <- ramh(n=1000, alpha=alpha), col="blue")
par(mfrow=c(1,2))
hist(r[,1]) # Should be uniform
hist(r[,2]) # Should be uniform
}
}
\keyword{distribution}



