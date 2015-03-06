\name{Inv.paralogistic}
\alias{Inv.paralogistic}
\alias{dinv.paralogistic}
\alias{pinv.paralogistic}
\alias{qinv.paralogistic}
\alias{rinv.paralogistic}
\title{The Inverse Paralogistic Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the inverse paralogistic distribution with
  shape parameters \code{a} and \code{p}, and scale parameter \code{scale}.


}
\usage{
dinv.paralogistic(x, scale = 1, shape1.a, log = FALSE)
pinv.paralogistic(q, scale = 1, shape1.a, lower.tail = TRUE, log.p = FALSE)
qinv.paralogistic(p, scale = 1, shape1.a, lower.tail = TRUE, log.p = FALSE)
rinv.paralogistic(n, scale = 1, shape1.a)
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
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }

}
\value{
  \code{dinv.paralogistic} gives the density,
  \code{pinv.paralogistic} gives the distribution function,
  \code{qinv.paralogistic} gives the quantile function, and
  \code{rinv.paralogistic} generates random deviates.


}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.

}
\author{ T. W. Yee }
\details{
  See \code{\link{inv.paralogistic}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.


}
\note{
  The inverse paralogistic distribution is a special case of the 4-parameter
  generalized beta II distribution.


}
\seealso{
  \code{\link{inv.paralogistic}},
  \code{\link{genbetaII}}.


}
\examples{
idata <- data.frame(y = rinv.paralogistic(n = 3000, exp(1), scale = exp(2)))
fit <- vglm(y ~ 1, inv.paralogistic(lss = FALSE, ishape1.a = 2.1),
            data = idata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{distribution}
