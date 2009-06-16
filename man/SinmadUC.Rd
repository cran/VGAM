\name{Sinmad}
\alias{Sinmad}
\alias{dsinmad}
\alias{psinmad}
\alias{qsinmad}
\alias{rsinmad}
\title{The Singh-Maddala Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Singh-Maddala distribution with shape parameters \code{a}
  and \code{q}, and scale parameter \code{scale}.
}
\usage{
dsinmad(x, a, scale=1, q.arg, log=FALSE)
psinmad(q, a, scale=1, q.arg)
qsinmad(p, a, scale=1, q.arg)
rsinmad(n, a, scale=1, q.arg)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{a, q.arg}{shape parameters.}
  \item{scale}{scale parameter.}
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dsinmad} gives the density,
  \code{psinmad} gives the distribution function,
  \code{qsinmad} gives the quantile function, and
  \code{rsinmad} generates random deviates.
}
\references{
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ: Wiley-Interscience.
}
\author{ T. W. Yee }
\details{
  See \code{\link{sinmad}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.
}
\note{
  The Singh-Maddala distribution is a special case of the 4-parameter
  generalized beta II distribution.
}
\seealso{
  \code{\link{sinmad}},
  \code{\link{genbetaII}}.
}
\examples{
y = rsinmad(n=3000, 4, 6, 2)
fit = vglm(y ~ 1, sinmad(init.a=2.1), trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
}
\keyword{distribution}
