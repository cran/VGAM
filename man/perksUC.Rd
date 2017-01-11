\name{Perks}
\alias{Perks}
\alias{dperks}
\alias{pperks}
\alias{qperks}
\alias{rperks}
\title{The Perks Distribution}
\description{
  Density, cumulative distribution function,
  quantile function
  and
  random generation for
  the Perks distribution.

}
\usage{
dperks(x, scale = 1, shape, log = FALSE)
pperks(q, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
qperks(p, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
rperks(n, scale = 1, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
           Same as in \code{\link[stats]{runif}}.

  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }

  \item{shape, scale}{positive shape and scale parameters. }

}
\value{
  \code{dperks} gives the density,
  \code{pperks} gives the cumulative distribution function,
  \code{qperks} gives the quantile function, and
  \code{rperks} generates random deviates.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{perks}} for details.


}
%\note{
%
%}
\seealso{
  \code{\link{perks}}.


}
\examples{
probs <- seq(0.01, 0.99, by = 0.01)
Shape <- exp(-1.0); Scale <- exp(1);
max(abs(pperks(qperks(p = probs, Shape, Scale),
                  Shape, Scale) - probs))  # Should be 0

\dontrun{ x <- seq(-0.1, 07, by = 0.01);
plot(x, dperks(x, Shape, Scale), type = "l", col = "blue", las = 1,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles",
     ylab = "", ylim = 0:1)
abline(h = 0, col = "blue", lty = 2)
lines(x, pperks(x, Shape, Scale), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qperks(probs, Shape, Scale)
lines(Q, dperks(Q, Shape, Scale), col = "purple", lty = 3, type = "h")
pperks(Q, Shape, Scale) - probs  # Should be all zero
abline(h = probs, col = "purple", lty = 3) }
}
\keyword{distribution}


