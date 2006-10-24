\name{gpdUC}
\alias{gpdUC}
\alias{dgpd}
\alias{pgpd}
\alias{qgpd}
\alias{rgpd}
\title{The Generalized Pareto Distribution }
\description{
  Density, distribution function, quantile function and random
  generation for the generalized Pareto distribution (GPD) with
  location parameter \code{location},
  scale parameter \code{scale} and
  shape parameter \code{shape}.
}
\usage{
dgpd(x, location=0, scale=1, shape=0)
pgpd(q, location=0, scale=1, shape=0)
qgpd(p, location=0, scale=1, shape=0)
rgpd(n, location=0, scale=1, shape=0)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Positive integer of length 1.}
  \item{location}{the location parameter \eqn{\mu}{mu}.}
  \item{scale}{the scale parameter \eqn{\sigma}{sigma}.}
  \item{shape}{the shape parameter \eqn{\xi}{xi}.}
}
\value{
  \code{dgpd} gives the density,
  \code{pgpd} gives the distribution function,
  \code{qgpd} gives the quantile function, and
  \code{rgpd} generates random deviates.
}
\references{
Coles, S. (2001)
\emph{An Introduction to Statistical Modeling of Extreme Values}.
London: Springer-Verlag.
}
\author{ T. W. Yee }
\details{
  See \code{\link{gpd}}, the \pkg{VGAM} family function
  for estimating the two parameters by maximum likelihood estimation,
  for formulae and other details.
  Apart from \code{n}, all the above arguments may be vectors and
  are recyled to the appropriate length if necessary.
}
\note{
  The default values of all three parameters, especially
  \eqn{\xi=0}{xi=0}, means the default distribution is the exponential.

  Currently, these functions have different argument names compared
  with those in the \pkg{evd} package.
}
\seealso{
  \code{\link{gpd}}.
}
\examples{
\dontrun{
x = seq(-0.2, 3, by=0.01)
loc = 0; sigma = 1; xi = -0.4
plot(x, dgpd(x, loc, sigma, xi), type="l", col="blue", ylim=c(0,1),
     main="Blue is density, red is cumulative distribution function",
     sub="Purple are 5,10,...,95 percentiles", ylab="", las=1)
abline(h=0, col="blue", lty=2)
lines(qgpd(seq(0.05,0.95,by=0.05), loc, sigma, xi), 
      dgpd(qgpd(seq(0.05,0.95,by=0.05), loc, sigma, xi), loc, sigma, xi),
      col="purple", lty=3, type="h")
lines(x, pgpd(x, loc, sigma, xi), type="l", col="red")
abline(h=0, lty=2)
}
}
\keyword{distribution}


