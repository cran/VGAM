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
dgev(x, location = 0, scale = 1, shape = 0, log = FALSE,
     tolshape0 = sqrt(.Machine$double.eps))
pgev(q, location = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
qgev(p, location = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
rgev(n, location = 0, scale = 1, shape = 0)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required. }
  \item{location}{the location parameter \eqn{\mu}{mu}.}
  \item{scale}{the (positive) scale parameter \eqn{\sigma}{sigma}.
      Must consist of positive values. }
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
  this value then it will be assumed zero and a Gumbel distribution will
  be used.


  }



% 20160412; Depreciated:
% \item{oobounds.log, giveWarning}{
% Numeric and logical.
% The GEV distribution has support in the region satisfying
% \code{1+shape*(x-location)/scale > 0}. Outside that region, the
% logarithm of the density is assigned \code{oobounds.log}, which
% equates to a zero density.
% It should not be assigned a positive number, and ideally is very negative.
% Since \code{\link{egev}} uses this function it is necessary
% to return a finite value outside this region  so as to allow
% for half-stepping. Both arguments are in support of this.
% This argument and others match those of \code{\link{egev}}.
% }




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
  for estimating the 3 parameters by maximum likelihood estimation,
  for formulae and other details.
  Apart from \code{n}, all the above arguments may be vectors and
  are recyled to the appropriate length if necessary.


}
\note{
  The default value of \eqn{\xi = 0}{xi = 0} means the default distribution
  is the Gumbel.


  Currently, these functions have different argument names compared with
  those in the \pkg{evd} package.


}
\seealso{
  \code{\link{gev}},
  \code{\link{gevff}},
  \code{\link{vglm.control}}.


}
\examples{ loc <- 2; sigma <- 1; xi <- -0.4
pgev(qgev(seq(0.05, 0.95, by = 0.05), loc, sigma, xi), loc, sigma, xi)
\dontrun{ x <- seq(loc - 3, loc + 3, by = 0.01)
plot(x, dgev(x, loc, sigma, xi), type = "l", col = "blue", ylim = c(0,1),
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple are 5,10,...,95 percentiles", ylab = "", las = 1)
abline(h = 0, col = "blue", lty = 2)
lines(qgev(seq(0.05, 0.95, by = 0.05), loc, sigma, xi),
      dgev(qgev(seq(0.05, 0.95, by = 0.05), loc, sigma, xi), loc, sigma, xi),
      col = "purple", lty = 3, type = "h")
lines(x, pgev(x, loc, sigma, xi), type = "l", col = "orange")
abline(h = 0, lty = 2, col = "gray50")
}
}
\keyword{distribution}


%dgev(x, location = 0, scale = 1, shape = 0, log = FALSE, tolshape0 =
%     sqrt(.Machine$double.eps), oobounds.log = -Inf, giveWarning = FALSE)






