\name{Invparalogistic}
\alias{Invparalogistic}
\alias{dinvparalogistic}
\alias{pinvparalogistic}
\alias{qinvparalogistic}
\alias{rinvparalogistic}
\title{The Inverse Paralogistic Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the inverse paralogistic distribution with
  shape parameters \code{a} and \code{p}, and scale parameter \code{scale}.


}
\usage{
dinvparalogistic(x, shape1.a, scale = 1, log = FALSE)
pinvparalogistic(q, shape1.a, scale = 1)
qinvparalogistic(p, shape1.a, scale = 1)
rinvparalogistic(n, shape1.a, scale = 1)
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
  If \code{log = TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dinvparalogistic} gives the density,
  \code{pinvparalogistic} gives the distribution function,
  \code{qinvparalogistic} gives the quantile function, and
  \code{rinvparalogistic} generates random deviates.


}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.

}
\author{ T. W. Yee }
\details{
  See \code{\link{invparalogistic}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.


}
\note{
  The inverse paralogistic distribution is a special case of the 4-parameter
  generalized beta II distribution.


}
\seealso{
  \code{\link{invparalogistic}},
  \code{\link{genbetaII}}.


}
\examples{
idata <- data.frame(y = rinvparalogistic(n = 3000, exp(1), exp(2)))
fit <- vglm(y ~ 1, invparalogistic(ishape1.a = 2.1),
            data = idata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{distribution}
