\name{Inv.gaussian}
\alias{Inv.gaussian}
\alias{dinv.gaussian}
\alias{pinv.gaussian}
\alias{rinv.gaussian}
\title{The Inverse Gaussian Distribution}
\description{
  Density, distribution function and random generation
  for the inverse Gaussian distribution.

}
\usage{
dinv.gaussian(x, mu, lambda, log = FALSE)
pinv.gaussian(q, mu, lambda)
rinv.gaussian(n, mu, lambda)
}
\arguments{
  \item{x, q}{vector of quantiles.}
%%\item{p}{vector of probabilities.}
  \item{n}{number of observations.
    If \code{length(n) > 1} then the length
    is taken to be the number required.

  }
  \item{mu}{the mean parameter.}
  \item{lambda}{the \eqn{\lambda}{lambda} parameter.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dinv.gaussian} gives the density,
  \code{pinv.gaussian} gives the distribution function, and
  \code{rinv.gaussian} generates random deviates.


% \code{qinv.gaussian} gives the quantile function, and


}
\references{
Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
\emph{Continuous Univariate Distributions},
2nd edition,
Volume 1,
New York: Wiley.


  Taraldsen, G. and Lindqvist, B. H. (2005)
  The multiple roots simulation algorithm,
  the inverse Gaussian distribution, and the
  sufficient conditional Monte Carlo method.
  \emph{Preprint Statistics No. 4/2005},
  Norwegian University of Science and Technology,
  Trondheim, Norway.

}
\author{ T. W. Yee }
\details{
  See \code{\link{inv.gaussianff}}, the \pkg{VGAM} family function
  for estimating both parameters by maximum likelihood estimation,
  for the formula of the probability density function.

}
\note{
  Currently \code{qinv.gaussian} is unavailable.

}
\seealso{
  \code{\link{inv.gaussianff}},
  \code{\link{waldff}}.


}
\examples{
\dontrun{ x <- seq(-0.05, 4, len = 300)
plot(x, dinv.gaussian(x, mu = 1, lambda = 1), type = "l",
     col = "blue",las = 1, main =
     "blue is density, orange is cumulative distribution function")
abline(h = 0, col = "gray", lty = 2)
lines(x, pinv.gaussian(x, mu = 1, lambda = 1), type = "l", col = "orange") }
}
\keyword{distribution}


