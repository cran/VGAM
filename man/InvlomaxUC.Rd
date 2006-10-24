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
dinvlomax(x, scale, p.arg)
pinvlomax(q, scale, p.arg)
qinvlomax(p, scale, p.arg)
rinvlomax(n, scale, p.arg)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{p.arg}{shape parameter.}
  \item{scale}{scale parameter.}
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
Hoboken, NJ: Wiley-Interscience.
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
y = rinvlomax(n=1000, 6, 2)
fit = vglm(y ~ 1, invlomax, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
}
\keyword{distribution}
