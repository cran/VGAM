\name{Dagum}
\alias{Dagum}
\alias{ddagum}
\alias{pdagum}
\alias{qdagum}
\alias{rdagum}
\title{The Dagum Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Dagum distribution with shape parameters \code{a}
  and \code{p}, and scale parameter \code{scale}.
}
\usage{
ddagum(x, a, scale, p.arg, log=FALSE)
pdagum(q, a, scale, p.arg)
qdagum(p, a, scale, p.arg)
rdagum(n, a, scale, p.arg)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{a, p.arg}{shape parameters.}
  \item{scale}{scale parameter.}
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{ddagum} gives the density,
  \code{pdagum} gives the distribution function,
  \code{qdagum} gives the quantile function, and
  \code{rdagum} generates random deviates.
}
\references{
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ: Wiley-Interscience.
}
\author{ T. W. Yee }
\details{
  See \code{\link{dagum}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.
}
\note{
  The Dagum distribution is a special case of the 4-parameter
  generalized beta II distribution.
}
\seealso{
  \code{\link{dagum}},
  \code{\link{genbetaII}}.
}
\examples{
y = rdagum(n=3000, 4, 6, 2)
fit = vglm(y ~ 1, dagum(init.a=2.1), trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
}
\keyword{distribution}
