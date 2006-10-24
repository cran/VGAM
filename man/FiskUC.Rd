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
dfisk(x, a, scale)
pfisk(q, a, scale)
qfisk(p, a, scale)
rfisk(n, a, scale)
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
  \code{dfisk} gives the density,
  \code{pfisk} gives the distribution function,
  \code{qfisk} gives the quantile function, and
  \code{rfisk} generates random deviates.
}
\references{
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ: Wiley-Interscience.
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
y = rfisk(n=1000, 4, 6)
fit = vglm(y ~ 1, fisk, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
}
\keyword{distribution}
