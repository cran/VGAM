\name{Plackett}
\alias{Plackett}
\alias{dplack}
\alias{pplack}
\alias{rplack}
\title{Plackett's Bivariate Distribution}
\description{
  Density, distribution function, and random
  generation for the (one parameter) bivariate Plackett distribution.

}
\usage{
dplack(x1, x2, oratio, log=FALSE)
pplack(q1, q2, oratio)
rplack(n, oratio)
}
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{oratio}{the positive odds ratio \eqn{\psi}{psi}.}
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm is returned.

  }
}
\value{
  \code{dplack} gives the density,
  \code{pplack} gives the distribution function, and
  \code{rplack} generates random deviates (a two-column matrix).
}
\references{

Mardia, K. V. (1967)
Some contributions to contingency-type distributions.
\emph{Biometrika},
\bold{54}, 235--249.

}
\author{ T. W. Yee }
\details{
  See \code{\link{plackett}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation, for the formula of the
  cumulative distribution function and other details.

}
%\note{
%}
\seealso{
  \code{\link{plackett}},
  \code{\link{frank}}.

}
\examples{
\dontrun{
N = 101
x = seq(0.0, 1.0, len=N)
oratio = exp(1)
ox = expand.grid(x, x)
z = dplack(ox[,1], ox[,2], oratio=oratio)
contour(x, x, matrix(z, N, N), col="blue")
z = pplack(ox[,1], ox[,2], oratio=oratio)
contour(x, x, matrix(z, N, N), col="blue")

plot(r <- rplack(n=3000, oratio=oratio))
par(mfrow=c(1,2))
hist(r[,1]) # Should be uniform
hist(r[,2]) # Should be uniform
}
}
\keyword{distribution}


