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
dgpd(x, location = 0, scale = 1, shape = 0, log = FALSE,
     tolshape0 = sqrt(.Machine$double.eps))
pgpd(q, location = 0, scale = 1, shape = 0,
     lower.tail = TRUE, log.p = FALSE)
qgpd(p, location = 0, scale = 1, shape = 0,
     lower.tail = TRUE, log.p = FALSE)
rgpd(n, location = 0, scale = 1, shape = 0)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.}
  \item{location}{the location parameter \eqn{\mu}{mu}.}
  \item{scale}{the (positive) scale parameter \eqn{\sigma}{sigma}.}
  \item{shape}{the shape parameter \eqn{\xi}{xi}.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Uniform]{punif}}
  or \code{\link[stats:Uniform]{qunif}}.


  }
  \item{tolshape0}{
  Positive numeric.
  Threshold/tolerance value for resting whether \eqn{\xi}{xi} is zero.
  If the absolute value of the estimate of \eqn{\xi}{xi} is less than
  this value then it will be assumed zero and an exponential distribution will
  be used.


  }


% \item{oobounds.log, giveWarning}{
% Numeric and logical.
% The GPD distribution has support in the region satisfying
% \code{(x-location)/scale > 0}
% and
% \code{1+shape*(x-location)/scale > 0}.
% Outside that region, the
% logarithm of the density is assigned \code{oobounds.log}, which
% equates to a zero density.
% It should not be assigned a positive number, and ideally is very negative.
% Since \code{\link{gpd}} uses this function it is necessary
% to return a finite value outside this region  so as to allow
% for half-stepping. Both arguments are in support of this.
% This argument and others match those of \code{\link{gpd}}.
% }

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
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{gpd}}, the \pkg{VGAM} family function
  for estimating the two parameters by maximum likelihood estimation,
  for formulae and other details.
  Apart from \code{n}, all the above arguments may be vectors and
  are recyled to the appropriate length if necessary.


}
\note{
  The default values of all three parameters, especially
  \eqn{\xi = 0}{xi = 0}, means the default distribution is the exponential.


  Currently, these functions have different argument names compared
  with those in the \pkg{evd} package.


}
\seealso{
  \code{\link{gpd}},
  \code{\link[stats]{Exponential}}.


}
\examples{
\dontrun{ loc <- 2; sigma <- 1; xi <- -0.4
x <- seq(loc - 0.2, loc + 3, by = 0.01)
plot(x, dgpd(x, loc, sigma, xi), type = "l", col = "blue", ylim = c(0, 1),
     main = "Blue is density, red is cumulative distribution function",
     sub = "Purple are 5,10,...,95 percentiles", ylab = "", las = 1)
abline(h = 0, col = "blue", lty = 2)
lines(qgpd(seq(0.05, 0.95, by = 0.05), loc, sigma, xi),
      dgpd(qgpd(seq(0.05, 0.95, by = 0.05), loc, sigma, xi), loc, sigma, xi),
      col = "purple", lty = 3, type = "h")
lines(x, pgpd(x, loc, sigma, xi), type = "l", col = "red")
abline(h = 0, lty = 2)

pgpd(qgpd(seq(0.05, 0.95, by = 0.05), loc, sigma, xi), loc, sigma, xi)
}
}
\keyword{distribution}


%    oobounds.log = -Inf, giveWarning = FALSE
