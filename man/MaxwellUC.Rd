\name{Maxwell}
\alias{Maxwell}
\alias{dmaxwell}
\alias{pmaxwell}
\alias{qmaxwell}
\alias{rmaxwell}
\title{The Maxwell Distribution}
\description{
  Density,
  distribution function,
  quantile function and random generation
  for the Maxwell distribution.

}
\usage{
dmaxwell(x, a, log = FALSE)
pmaxwell(q, a)
qmaxwell(p, a)
rmaxwell(n, a)
}
\arguments{
  \item{x, q, p, n}{
  Same as \code{\link[stats:Uniform]{Uniform}}.


  }
  \item{a}{the parameter.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dmaxwell} gives the density,
  \code{pmaxwell} gives the distribution function,
  \code{qmaxwell} gives the quantile function, and
  \code{rmaxwell} generates random deviates.


}
\references{

  Balakrishnan, N. and Nevzorov, V. B. (2003) 
  \emph{A Primer on Statistical Distributions}.
  Hoboken, New Jersey: Wiley.


}
\author{ T. W. Yee }
\details{
  See \code{\link{maxwell}}, the \pkg{VGAM} family function
  for estimating the parameter \eqn{a} by maximum likelihood estimation,
  for the formula of the probability density function.


}
\note{
  The Maxwell distribution is related to the Rayleigh distribution.


}
\seealso{
  \code{\link{maxwell}},
  \code{\link{Rayleigh}},
  \code{\link{rayleigh}}.


}
\examples{
\dontrun{ a <- 3; x <- seq(-0.5, 3, length = 100)
plot(x, dmaxwell(x, a = a), type = "l", col = "blue", las = 1, ylab = "",
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles")
abline(h = 0, col = "blue", lty = 2)
lines(x, pmaxwell(x, a = a), type = "l", col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qmaxwell(probs, a = a)
lines(Q, dmaxwell(Q, a), col = "purple", lty = 3, type = "h")
lines(Q, pmaxwell(Q, a), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(pmaxwell(Q, a) - probs))  # Should be zero
}
}
\keyword{distribution}


