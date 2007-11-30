\name{Frank}
\alias{Frank}
\alias{dfrank}
\alias{pfrank}
\alias{rfrank}
\title{Frank's Bivariate Distribution}
\description{
  Density, distribution function, and random
  generation for the (one parameter) bivariate Frank distribution.

}
\usage{
dfrank(x1, x2, alpha)
pfrank(q1, q2, alpha)
rfrank(n, alpha)
}
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{alpha}{the positive association parameter \eqn{\alpha}{alpha}.}
}
\value{
  \code{dfrank} gives the density,
  \code{pfrank} gives the distribution function, and
  \code{rfrank} generates random deviates (a two-column matrix).
}
\references{

Genest, C. (1987)
Frank's family of bivariate distributions.
\emph{Biometrika},
\bold{74}, 549--555.

}
\author{ T. W. Yee }
\details{
  See \code{\link{frank}}, the \pkg{VGAM}
  family functions for estimating the association
  parameter by maximum likelihood estimation, for the formula of the
  cumulative distribution function and other details.

}
%\note{
%}
\seealso{
  \code{\link{frank}}.
}
\examples{
\dontrun{
N = 100
x = seq(-0.30, 1.30, len=N)
alpha = 8
ox = expand.grid(x, x)
z = dfrank(ox[,1], ox[,2], alp=alpha)
contour(x, x, matrix(z, N, N))
z = pfrank(ox[,1], ox[,2], alp=alpha)
contour(x, x, matrix(z, N, N))

alpha = exp(4)
plot(r <- rfrank(n=3000, alpha=alpha))
par(mfrow=c(1,2))
hist(r[,1]) # Should be uniform
hist(r[,2]) # Should be uniform
}
}
\keyword{distribution}


