\name{Benini}
\alias{Benini}
\alias{dbenini}
\alias{pbenini}
\alias{qbenini}
\alias{rbenini}
\title{The Benini Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Benini distribution with parameter
  \code{shape}.


}
\usage{
dbenini(x, y0, shape, log = FALSE)
pbenini(q, y0, shape, lower.tail = TRUE, log.p = FALSE)
qbenini(p, y0, shape, lower.tail = TRUE, log.p = FALSE)
rbenini(n, y0, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Same as \code{\link[stats]{runif}}.


  }
  \item{y0}{the scale parameter \eqn{y_0}{y0}.


  }
  \item{shape}{the shape parameter \eqn{b}.


  }
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
  \code{dbenini} gives the density,
  \code{pbenini} gives the distribution function,
  \code{qbenini} gives the quantile function, and
  \code{rbenini} generates random deviates.


}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{benini1}}, the \pkg{VGAM} family function
  for estimating the parameter \eqn{s} by maximum likelihood estimation,
  for the formula of the probability density function and other details.


}
%\note{
%
%}
\seealso{
  \code{\link{benini1}}.


}
\examples{
\dontrun{
y0 <- 1; shape <- exp(1)
xx <- seq(0.0, 4, len = 101)
plot(xx, dbenini(xx, y0 = y0, shape = shape), type = "l", col = "blue",
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles", ylim = 0:1,
     las = 1, ylab = "", xlab = "x")
abline(h = 0, col = "blue", lty = 2)
lines(xx, pbenini(xx, y0 = y0, shape = shape), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qbenini(probs, y0 = y0, shape = shape)
lines(Q, dbenini(Q, y0 = y0, shape = shape),
      col = "purple", lty = 3, type = "h")
pbenini(Q, y0 = y0, shape = shape) - probs  # Should be all zero
}
}
\keyword{distribution}


