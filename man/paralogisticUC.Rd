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
dparalogistic(x, shape1.a, scale = 1, log = FALSE)
pparalogistic(q, shape1.a, scale = 1)
qparalogistic(p, shape1.a, scale = 1)
rparalogistic(n, shape1.a, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape1.a}{shape parameter.}
  \item{scale}{scale parameter.}
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.

  }

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
Hoboken, NJ, USA: Wiley-Interscience.


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
pdata <- data.frame(y = rparalogistic(n = 3000, exp(1), exp(2)))
fit <- vglm(y ~ 1, paralogistic(ishape1.a = 2.1), data = pdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{distribution}
