\name{Fisk}
\alias{Fisk}
\alias{dfisk}
\alias{pfisk}
\alias{qfisk}
\alias{rfisk}
\title{The Fisk Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Fisk distribution with shape parameter \code{a}
  and scale parameter \code{scale}.
}
\usage{
dfisk(x, shape1.a, scale = 1, log = FALSE)
pfisk(q, shape1.a, scale = 1)
qfisk(p, shape1.a, scale = 1)
rfisk(n, shape1.a, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.}
  \item{shape1.a}{shape parameter.}
  \item{scale}{scale parameter.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dfisk} gives the density,
  \code{pfisk} gives the distribution function,
  \code{qfisk} gives the quantile function, and
  \code{rfisk} generates random deviates.

}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


}
\author{ T. W. Yee }
\details{
  See \code{\link{fisk}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.


}
\note{
  The Fisk distribution is a special case of the 4-parameter
  generalized beta II distribution.


}
\seealso{
  \code{\link{fisk}},
  \code{\link{genbetaII}}.

}
\examples{
fdata = data.frame(y = rfisk(n = 1000, 4, 6))
fit = vglm(y ~ 1, fisk, data = fdata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{distribution}
