\name{Inv.gaussian}
\alias{Inv.gaussian}
\alias{dinv.gaussian}
\alias{pinv.gaussian}
\title{The Inverse Gaussian Distribution}
\description{
  Density and distribution function
  for the inverse Gaussian distribution with parameters
  \code{mu} and \code{lambda}.
}
\usage{
dinv.gaussian(x, mu, lambda)
pinv.gaussian(q, mu, lambda)
}
\arguments{
  \item{x, q}{vector of quantiles.}
%%\item{p}{vector of probabilities.}
%%\item{n}{number of observations. If \code{length(n) > 1}, the length
%%  is taken to be the number required.}
  \item{mu}{the mean parameter.}
  \item{lambda}{the \eqn{\lambda}{lambda} parameter.}
}
\value{
  \code{dinv.gaussian} gives the density,
  \code{pinv.gaussian} gives the distribution function.
% \code{qinv.gaussian} gives the quantile function, and
% \code{rinv.gaussian} generates random deviates.
}
\references{
Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
\emph{Continuous Univariate Distributions},
2nd edition,
Volume 1,
New York: Wiley.
}
\author{ T. W. Yee }
\details{
  See \code{\link{inv.gaussianff}}, the \pkg{VGAM} family function
  for estimating both parameters by maximum likelihood estimation,
  for the formula of the probability density function.
}
\note{
  Currently \code{qinv.gaussian} and \code{rinv.gaussian} are unavailable.
}
\seealso{
  \code{\link{inv.gaussianff}}.
}
\examples{
\dontrun{
x = seq(-0.05, 4, len=300)
plot(x, dinv.gaussian(x, mu=1, lambda=1), type="l", col="blue", las=1,
     main="blue is density, red is cumulative distribution function")
abline(h=0, col="blue", lty=2)
lines(x, pinv.gaussian(x, mu=1, lambda=1), type="l", col="red")
}
}
\keyword{distribution}


