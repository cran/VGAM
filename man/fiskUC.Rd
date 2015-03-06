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
dfisk(x, scale = 1, shape1.a, log = FALSE)
pfisk(q, scale = 1, shape1.a, lower.tail = TRUE, log.p = FALSE)
qfisk(p, scale = 1, shape1.a, lower.tail = TRUE, log.p = FALSE)
rfisk(n, scale = 1, shape1.a)
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
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


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
\author{ T. W. Yee and Kai Huang }
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
fdata <- data.frame(y = rfisk(n = 1000, shape = exp(1), scale = exp(2)))
fit <- vglm(y ~ 1, fisk(lss = FALSE), data = fdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{distribution}
