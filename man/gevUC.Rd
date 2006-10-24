\name{gevUC}
\alias{gevUC}
\alias{dgev}
\alias{pgev}
\alias{qgev}
\alias{rgev}
\title{The Generalized Extreme Value Distribution }
\description{
  Density, distribution function, quantile function and random
  generation for the generalized extreme value distribution (GEV) with
  location parameter \code{location},
  scale parameter \code{scale} and
  shape parameter \code{shape}.
}
\usage{
dgev(x, location=0, scale=1, shape=0)
pgev(q, location=0, scale=1, shape=0)
qgev(p, location=0, scale=1, shape=0)
rgev(n, location=0, scale=1, shape=0)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Positive integer of length 1.}
  \item{location}{the location parameter \eqn{\mu}{mu}.}
  \item{scale}{the scale parameter \eqn{\sigma}{sigma}.
      Must consist of positive values. }
  \item{shape}{the shape parameter \eqn{\xi}{xi}.}
}
\value{
  \code{dgev} gives the density,
  \code{pgev} gives the distribution function,
  \code{qgev} gives the quantile function, and
  \code{rgev} generates random deviates.
}
\references{
Coles, S. (2001)
\emph{An Introduction to Statistical Modeling of Extreme Values}.
London: Springer-Verlag.
}
\author{ T. W. Yee }
\details{
  See \code{\link{gev}}, the \pkg{VGAM} family function
  for estimating the two parameters by maximum likelihood estimation,
  for formulae and other details.
  Apart from \code{n}, all the above arguments may be vectors and
  are recyled to the appropriate length if necessary.
}
\note{
  The default value of \eqn{\xi=0}{xi=0} means the default distribution
  is the Gumbel.

  Currently, these functions have different argument names compared with
  those in the \pkg{evd} package.

}
\seealso{
  \code{\link{gev}}.
}
\examples{
\dontrun{
x = seq(-3, 3, by=0.01)
loc = 0; sigma = 1; xi = -0.4
plot(x, dgev(x, loc, sigma, xi), type="l", col="blue", ylim=c(0,1),
     main="Blue is density, red is cumulative distribution function",
     sub="Purple are 5,10,...,95 percentiles", ylab="", las=1)
abline(h=0, col="blue", lty=2)
lines(qgev(seq(0.05,0.95,by=0.05), loc, sigma, xi), 
      dgev(qgev(seq(0.05,0.95,by=0.05), loc, sigma, xi), loc, sigma, xi),
      col="purple", lty=3, type="h")
lines(x, pgev(x, loc, sigma, xi), type="l", col="red")
abline(h=0, lty=2)

pgev(qgev(seq(0.05,0.95,by=0.05), loc, sigma, xi), loc, sigma, xi)
}
}
\keyword{distribution}


