\name{skewnorm}
\alias{skewnorm}
\alias{dskewnorm}
%\alias{pskewnorm}
%\alias{qskewnorm}
\alias{rskewnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Skew-Normal Distribution }
\description{
  Density and
  random generation
  for the univariate skew-normal distribution.


% , distribution function, quantile function and


}
\usage{
dskewnorm(x, location = 0, scale = 1, shape = 0, log = FALSE)
rskewnorm(n, location = 0, scale = 1, shape = 0)
}
%pskewnorm(q, lambda)
%qskewnorm(p, lambda)
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of quantiles.}
% \item{x, q}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as \code{\link[stats]{runif}}.

  }

  \item{location}{
  The location parameter \eqn{\xi}{xi}. A vector.


  }
  \item{scale}{
  The scale parameter \eqn{\omega}{w}. A positive vector.


  }
  \item{shape}{
  The shape parameter. It is called \eqn{\alpha}{alpha} in
  \code{\link{skewnormal}}.


  }
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.


  }

}
\details{
  See \code{\link{skewnormal}}, which currently only estimates the shape
  parameter.
  More generally here, \eqn{Z = \xi + \omega Y}{Z = xi + w * Y} where
  \eqn{Y} has a standard skew-normal distribution
  (see \code{\link{skewnormal}}),
  \eqn{\xi}{xi} is the location parameter and
  \eqn{\omega}{w} is the scale parameter.


}
\value{
  \code{dskewnorm} gives the density,
  \code{rskewnorm} generates random deviates.


% \code{pskewnorm} gives the distribution function,
% \code{qskewnorm} gives the quantile function, and


}
\references{
    \url{http://tango.stat.unipd.it/SN}.


}

\author{ T. W. Yee }
\note{
  The default values of all three parameters corresponds to the
  skew-normal being the standard normal distribution.


}

\seealso{
  \code{\link{skewnormal}}.


}
\examples{
\dontrun{ N <- 200  # Grid resolution
shape <- 7; x <- seq(-4, 4, len = N)
plot(x, dskewnorm(x, shape = shape), type = "l", col = "blue", las = 1,
     ylab = "", lty = 1, lwd = 2)
abline(v = 0, h = 0, col = "grey")
lines(x, dnorm(x), col = "orange", lty = 2, lwd = 2)
legend("topleft", leg = c(paste("Blue = dskewnorm(x, ", shape,")", sep = ""),
       "Orange = standard normal density"), lty = 1:2, lwd = 2,
       col = c("blue", "orange")) }
}
\keyword{distribution}

