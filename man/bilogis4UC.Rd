\name{bilogis4}
\alias{bilogis4}
\alias{dbilogis4}
\alias{pbilogis4}
\alias{rbilogis4}
\title{Bivariate Logistic Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the 4-parameter bivariate logistic distribution.

}
\usage{
dbilogis4(x1, x2, loc1=0, scale1=1, loc2=0, scale2=1, log=FALSE)
pbilogis4(q1, q2, loc1=0, scale1=1, loc2=0, scale2=1)
rbilogis4(n, loc1=0, scale1=1, loc2=0, scale2=1)
}
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{loc1, loc2}{the location parameters \eqn{l_1}{l1} and \eqn{l_2}{l2}.}
  \item{scale1, scale2}{the scale parameters \eqn{s_1}{s1} and \eqn{s_2}{s2}.}
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dbilogis4} gives the density,
  \code{pbilogis4} gives the distribution function, and
  \code{rbilogis4} generates random deviates (a two-column matrix).
}
\references{

Gumbel, E. J. (1961)
Bivariate logistic distributions.
\emph{Journal of the American Statistical Association},
\bold{56}, 335--349.

}
\author{ T. W. Yee }
\details{
  See \code{\link{bilogis4}}, the \pkg{VGAM} family function for
  estimating the four parameters by maximum likelihood estimation, for
  the formula of the cumulative distribution function and other details.

}
%\note{
%}
\seealso{
  \code{\link{bilogistic4}}.
}
\examples{
\dontrun{
par(mfrow=c(1,3))
n = 2000
ymat = rbilogis4(n, loc1=5, loc2=7, scale2=exp(1))
myxlim = c(-2,15)
myylim = c(-10,30)
plot(ymat, xlim=myxlim, ylim=myylim)

N = 100
x1 = seq(myxlim[1], myxlim[2], len=N)
x2 = seq(myylim[1], myylim[2], len=N)
ox = expand.grid(x1, x2)
z = dbilogis4(ox[,1], ox[,2], loc1=5, loc2=7, scale2=exp(1))
contour(x1, x2, matrix(z, N, N), main="density")
z = pbilogis4(ox[,1], ox[,2], loc1=5, loc2=7, scale2=exp(1))
contour(x1, x2, matrix(z, N, N), main="cdf")
}
}
\keyword{distribution}


