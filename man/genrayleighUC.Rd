\name{genray}
\alias{genray}
\alias{dgenray}
\alias{pgenray}
\alias{qgenray}
\alias{rgenray}
\title{The Generalized Rayleigh Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the generalized Rayleigh distribution.


}
\usage{
dgenray(x, scale = 1, shape, log = FALSE)
pgenray(q, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
qgenray(p, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
rgenray(n, scale = 1, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    If \code{length(n) > 1} then the length is taken to be the number required. }
  \item{scale, shape}{
   positive scale and shape parameters. }
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
  \code{dgenray} gives the density,
  \code{pgenray} gives the distribution function,
  \code{qgenray} gives the quantile function, and
  \code{rgenray} generates random deviates.


}

\author{ Kai Huang and J. G. Lauder and T. W. Yee }

\details{
  See \code{\link{genrayleigh}}, the \pkg{VGAM} family function
  for estimating the parameters, 
  for the formula of the probability density function and other details.


}
\note{
  We define \code{scale} as the reciprocal of the scale parameter
  used by Kundu and Raqab (2005).


}
\seealso{
  \code{\link{genrayleigh}},
  \code{\link{rayleigh}}.


}
\examples{
\dontrun{
shape <- 0.5; Scale <- 1; nn <- 501
x <- seq(-0.10, 3.0, len = nn)
plot(x, dgenray(x, shape, scale = Scale), type = "l", las = 1, ylim = c(0, 1.2),
     ylab = paste("[dp]genray(shape = ", shape, ", scale = ", Scale, ")"),
     col = "blue", cex.main = 0.8,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles")
lines(x, pgenray(x, shape, scale = Scale), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qgenray(probs, shape, scale = Scale)
lines(Q, dgenray(Q, shape, scale = Scale), col = "purple", lty = 3, type = "h")
lines(Q, pgenray(Q, shape, scale = Scale), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(pgenray(Q, shape, scale = Scale) - probs))  # Should be 0
}
}
\keyword{distribution}


