\name{Rayleigh}
\alias{Rayleigh}
\alias{drayleigh}
\alias{prayleigh}
\alias{qrayleigh}
\alias{rrayleigh}
\title{The Rayleigh Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Rayleigh distribution with parameter
  \code{a}.


}
\usage{
drayleigh(x, scale = 1, log = FALSE)
prayleigh(q, scale = 1, lower.tail = TRUE, log.p = FALSE)
qrayleigh(p, scale = 1, lower.tail = TRUE, log.p = FALSE)
rrayleigh(n, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Fed into \code{\link[stats]{runif}}.


  }
  \item{scale}{the scale parameter \eqn{b}.}
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
  \code{drayleigh} gives the density,
  \code{prayleigh} gives the distribution function,
  \code{qrayleigh} gives the quantile function, and
  \code{rrayleigh} generates random deviates.


}
\references{

Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{rayleigh}}, the \pkg{VGAM} family
  function for estimating the scale parameter \eqn{b} by
  maximum likelihood estimation, for the formula of the
  probability density function and range restrictions on
  the parameter \eqn{b}.


}
\note{
  The Rayleigh distribution is related to the Maxwell distribution.


}
\seealso{
  \code{\link{rayleigh}},
  \code{\link{maxwell}}.


}
\examples{
\dontrun{ Scale <- 2; x <- seq(-1, 8, by = 0.1)
plot(x, drayleigh(x, scale = Scale), type = "l", ylim = c(0,1),
     las = 1, ylab = "",
     main = "Rayleigh density divided into 10 equal areas; orange = cdf")
abline(h = 0, col = "blue", lty = 2)
qq <- qrayleigh(seq(0.1, 0.9, by = 0.1), scale = Scale)
lines(qq, drayleigh(qq, scale = Scale), col = "purple", lty = 3, type = "h")
lines(x, prayleigh(x, scale = Scale), col = "orange") }
}
\keyword{distribution}


