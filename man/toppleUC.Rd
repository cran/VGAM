\name{Topple}
\alias{Topple}
\alias{dtopple}
\alias{ptopple}
\alias{qtopple}
\alias{rtopple}
\title{The Topp-Leone Distribution}
\description{
  Density,
  distribution function,
  quantile function and random generation
  for the Topp-Leone distribution.

}
\usage{
dtopple(x, shape, log = FALSE)
ptopple(q, shape, lower.tail = TRUE, log.p = FALSE)
qtopple(p, shape)
rtopple(n, shape)
}
\arguments{
  \item{x, q, p, n}{
  Same as \code{\link[stats:Uniform]{Uniform}}.


  }
  \item{shape}{the (shape) parameter, which lies in \eqn{(0, 1)}.}
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
  \code{dtopple} gives the density,
  \code{ptopple} gives the distribution function,
  \code{qtopple} gives the quantile function, and
  \code{rtopple} generates random deviates.


}
\references{

  Topp, C. W. and F. C. Leone (1955)
  A family of J-shaped frequency functions.
  \emph{Journal of the American Statistical Association},
  \bold{50}, 209--219.



}
\author{ T. W. Yee }
\details{
  See \code{\link{topple}}, the \pkg{VGAM} family function for
  estimating the (shape) parameter \eqn{s} by maximum likelihood
  estimation, for the formula of the probability density function.


}
\note{
  The Topp-Leone distribution is related to the triangle distribution.


}
\seealso{
  \code{\link{topple}},
  \code{\link{Triangle}}.



}
\examples{
\dontrun{ shape <- 0.7; x <- seq(0.02, 0.999, length = 300)
plot(x, dtopple(x, shape = shape), type = "l", col = "blue", las = 1,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles", ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(x, ptopple(x, shape = shape), type = "l", col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qtopple(probs, shape = shape)
lines(Q, dtopple(Q, shape), col = "purple", lty = 3, type = "h")
lines(Q, ptopple(Q, shape), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(ptopple(Q, shape) - probs))  # Should be zero
}
}
\keyword{distribution}


