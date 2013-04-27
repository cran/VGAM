\name{snorm}
\alias{snorm}
\alias{dsnorm}
%\alias{psnorm}
%\alias{qsnorm}
\alias{rsnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Skew-Normal Distribution }
\description{
  Density and 
% , distribution function, quantile function and
  random generation
  for the skew-normal distribution.

}
\usage{
dsnorm(x, location = 0, scale = 1, shape = 0, log = FALSE)
%psnorm(q, lambda)
%qsnorm(p, lambda)
rsnorm(n, location = 0, scale = 1, shape = 0)
}
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
  \code{\link{skewnormal1}}.

  }
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.

  }

}
\details{
  See \code{\link{skewnormal1}}, which currently only estimates the shape
  parameter.
  More generally here, \eqn{Z = \xi + \omega Y}{Z = xi + w * Y} where
  \eqn{Y} has a standard skew-normal distribution
  (see \code{\link{skewnormal1}}),
  \eqn{\xi}{xi} is the location parameter and
  \eqn{\omega}{w} is the scale parameter.

}
\value{
  \code{dsnorm} gives the density,
% \code{psnorm} gives the distribution function,
% \code{qsnorm} gives the quantile function, and
  \code{rsnorm} generates random deviates.

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
  \code{\link{skewnormal1}}.


}
\examples{
\dontrun{ N <- 200 # grid resolution
shape <- 7
x <- seq(-4, 4, len = N)
plot(x, dsnorm(x, shape = shape), type = "l", col = "blue", las = 1,
     ylab = "", lty = 1, lwd = 2)
abline(v = 0, h = 0, col = "grey")
lines(x, dnorm(x), col = "orange", lty = 2, lwd = 2)
legend("topleft", leg = c(paste("Blue = dsnorm(x, ", shape,")", sep = ""),
       "Orange = standard normal density"), lty = 1:2, lwd = 2,
       col = c("blue", "orange")) }
}
\keyword{distribution}

