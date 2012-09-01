\name{Lomax}
\alias{Lomax}
\alias{dlomax}
\alias{plomax}
\alias{qlomax}
\alias{rlomax}
\title{The Lomax Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Lomax distribution with scale parameter \code{scale}
  and shape parameter \code{q}.
}
\usage{
dlomax(x, scale = 1, shape3.q, log = FALSE)
plomax(q, scale = 1, shape3.q)
qlomax(p, scale = 1, shape3.q)
rlomax(n, scale = 1, shape3.q)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{scale}{scale parameter.}
  \item{shape3.q}{shape parameter.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dlomax} gives the density,
  \code{plomax} gives the distribution function,
  \code{qlomax} gives the quantile function, and
  \code{rlomax} generates random deviates.
}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.

}
\author{ T. W. Yee }
\details{
  See \code{\link{lomax}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.

}
\note{
  The Lomax distribution is a special case of the 4-parameter
  generalized beta II distribution.

}
\seealso{
  \code{\link{lomax}},
  \code{\link{genbetaII}}.

}
\examples{
ldata <- data.frame(y = rlomax(n = 2000, 6, 2))
fit <- vglm(y ~ 1, lomax(ishape3.q = 2.1), ldata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{distribution}
