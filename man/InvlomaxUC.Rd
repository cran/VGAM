\name{Invlomax}
\alias{Invlomax}
\alias{dinvlomax}
\alias{pinvlomax}
\alias{qinvlomax}
\alias{rinvlomax}
\title{The Inverse Lomax Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the inverse Lomax distribution with shape parameter
  \code{p} and scale parameter \code{scale}.

}
\usage{
dinvlomax(x, scale = 1, shape2.p, log = FALSE)
pinvlomax(q, scale = 1, shape2.p)
qinvlomax(p, scale = 1, shape2.p)
rinvlomax(n, scale = 1, shape2.p)
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

}
\value{
  \code{dinvlomax} gives the density,
  \code{pinvlomax} gives the distribution function,
  \code{qinvlomax} gives the quantile function, and
  \code{rinvlomax} generates random deviates.

}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.

}
\author{ T. W. Yee }
\details{
  See \code{\link{invlomax}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.


}
\note{
  The inverse Lomax distribution is a special case of the 4-parameter
  generalized beta II distribution.


}
\seealso{
  \code{\link{invlomax}},
  \code{\link{genbetaII}}.


}
\examples{
idata <- data.frame(y = rinvlomax(n = 1000, exp(2), exp(1)))
fit <- vglm(y ~ 1, invlomax, idata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{distribution}
