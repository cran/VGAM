\name{Inv.lomax}
\alias{Inv.lomax}
\alias{dinv.lomax}
\alias{pinv.lomax}
\alias{qinv.lomax}
\alias{rinv.lomax}
\title{The Inverse Lomax Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the inverse Lomax distribution with shape parameter
  \code{p} and scale parameter \code{scale}.

}
\usage{
dinv.lomax(x, scale = 1, shape2.p, log = FALSE)
pinv.lomax(q, scale = 1, shape2.p, lower.tail = TRUE, log.p = FALSE)
qinv.lomax(p, scale = 1, shape2.p, lower.tail = TRUE, log.p = FALSE)
rinv.lomax(n, scale = 1, shape2.p)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape2.p}{shape parameter.}
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
  \code{dinv.lomax} gives the density,
  \code{pinv.lomax} gives the distribution function,
  \code{qinv.lomax} gives the quantile function, and
  \code{rinv.lomax} generates random deviates.



}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.

}
\author{ T. W. Yee }
\details{
  See \code{\link{inv.lomax}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.


}
\note{
  The inverse Lomax distribution is a special case of the 4-parameter
  generalized beta II distribution.


}
\seealso{
  \code{\link{inv.lomax}},
  \code{\link{genbetaII}}.


}
\examples{
idata <- data.frame(y = rinv.lomax(n = 1000, exp(2), exp(1)))
fit <- vglm(y ~ 1, inv.lomax, data = idata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{distribution}
