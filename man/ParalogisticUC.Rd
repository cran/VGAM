\name{Paralogistic}
\alias{Paralogistic}
\alias{dparalogistic}
\alias{pparalogistic}
\alias{qparalogistic}
\alias{rparalogistic}
\title{The Paralogistic Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the paralogistic distribution with shape parameter \code{a}
  and scale parameter \code{scale}.
}
\usage{
dparalogistic(x, a, scale)
pparalogistic(q, a, scale)
qparalogistic(p, a, scale)
rparalogistic(n, a, scale)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{a}{shape parameter.}
  \item{scale}{scale parameter.}
}
\value{
  \code{dparalogistic} gives the density,
  \code{pparalogistic} gives the distribution function,
  \code{qparalogistic} gives the quantile function, and
  \code{rparalogistic} generates random deviates.
}
\references{
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ: Wiley-Interscience.
}
\author{ T. W. Yee }
\details{
  See \code{\link{paralogistic}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.
}
\note{
  The paralogistic distribution is a special case of the 4-parameter
  generalized beta II distribution.
}
\seealso{
  \code{\link{paralogistic}},
  \code{\link{genbetaII}}.
}
\examples{
y = rparalogistic(n=3000, 4, 6)
fit = vglm(y ~ 1, paralogistic(init.a=2.1), trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
}
\keyword{distribution}
