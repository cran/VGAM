\name{Dagum}
\alias{Dagum}
\alias{ddagum}
\alias{pdagum}
\alias{qdagum}
\alias{rdagum}
\title{The Dagum Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Dagum distribution with shape parameters \code{a}
  and \code{p}, and scale parameter \code{scale}.


}
\usage{
ddagum(x, scale = 1, shape1.a, shape2.p, log = FALSE)
pdagum(q, scale = 1, shape1.a, shape2.p, lower.tail = TRUE, log.p = FALSE)
qdagum(p, scale = 1, shape1.a, shape2.p, lower.tail = TRUE, log.p = FALSE)
rdagum(n, scale = 1, shape1.a, shape2.p)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape1.a, shape2.p}{shape parameters.}
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
  \code{ddagum} gives the density,
  \code{pdagum} gives the distribution function,
  \code{qdagum} gives the quantile function, and
  \code{rdagum} generates random deviates.

}
\references{


Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{dagum}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.


}
\note{
  The Dagum distribution is a special case of the 4-parameter
  generalized beta II distribution.


}
\seealso{
  \code{\link{dagum}},
  \code{\link{genbetaII}}.


}
\examples{
probs <- seq(0.1, 0.9, by = 0.1)
shape1.a <- 1; shape2.p <- 2
# Should be 0:
max(abs(pdagum(qdagum(p = probs, shape1.a = shape1.a, shape2.p =  shape2.p),
                                 shape1.a = shape1.a, shape2.p = shape2.p) - probs))

\dontrun{ par(mfrow = c(1, 2))
x <- seq(-0.01, 5, len = 401)
plot(x, dexp(x), type = "l", col = "black", ylab = "", las = 1, ylim = c(0, 1),
     main = "Black is standard exponential, others are ddagum(x, ...)")
lines(x, ddagum(x, shape1.a = shape1.a, shape2.p = 1), col = "orange")
lines(x, ddagum(x, shape1.a = shape1.a, shape2.p = 2), col = "blue")
lines(x, ddagum(x, shape1.a = shape1.a, shape2.p = 5), col = "green")
legend("topright", col = c("orange","blue","green"), lty = rep(1, len = 3),
       legend = paste("shape1.a =", shape1.a, ", shape2.p =", c(1, 2, 5)))

plot(x, pexp(x), type = "l", col = "black", ylab = "", las = 1,
     main = "Black is standard exponential, others are pdagum(x, ...)")
lines(x, pdagum(x, shape1.a = shape1.a, shape2.p = 1), col = "orange")
lines(x, pdagum(x, shape1.a = shape1.a, shape2.p = 2), col = "blue")
lines(x, pdagum(x, shape1.a = shape1.a, shape2.p = 5), col = "green")
legend("bottomright", col = c("orange","blue","green"), lty = rep(1, len = 3),
       legend = paste("shape1.a =", shape1.a, ", shape2.p =", c(1, 2, 5)))
}
}
\keyword{distribution}
