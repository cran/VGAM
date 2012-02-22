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
dsinmad(x, shape1.a, scale = 1, shape3.q, log = FALSE)
psinmad(q, shape1.a, scale = 1, shape3.q)
qsinmad(p, shape1.a, scale = 1, shape3.q)
rsinmad(n, shape1.a, scale = 1, shape3.q)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape1.a, shape3.q}{shape parameters.}
  \item{scale}{scale parameter.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

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
sdata = data.frame(y = rsinmad(n = 3000, 4, 6, 2))
fit = vglm(y ~ 1, sinmad(ishape1.a = 2.1), sdata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{distribution}
