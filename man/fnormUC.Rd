\name{Fnorm}
\alias{Fnorm}
\alias{dfnorm}
\alias{pfnorm}
\alias{qfnorm}
\alias{rfnorm}
\title{The Folded-Normal Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the (generalized) folded-normal distribution.

}
\usage{
dfnorm(x, mean = 0, sd = 1, a1 = 1, a2 = 1)
pfnorm(q, mean = 0, sd = 1, a1 = 1, a2 = 1)
qfnorm(p, mean = 0, sd = 1, a1 = 1, a2 = 1, ...)
rfnorm(n, mean = 0, sd = 1, a1 = 1, a2 = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{mean, sd}{ see \code{\link[stats:Normal]{rnorm}}. }
  \item{a1, a2}{ see \code{\link{fnormal1}}. }
  \item{\ldots}{
  Arguments that can be passed into \code{\link[stats]{uniroot}}.

  }
}
\value{
  \code{dfnorm} gives the density,
  \code{pfnorm} gives the distribution function,
  \code{qfnorm} gives the quantile function, and
  \code{rfnorm} generates random deviates.
}
\author{ T. W. Yee }
\details{
  See \code{\link{fnormal1}}, the \pkg{VGAM} family function
  for estimating the parameters, 
  for the formula of the probability density function and other details.

}
\note{
  \code{qfnorm} runs very slowly because it calls
  \code{\link[stats]{uniroot}} for each value of the argument \code{p}.
  The solution is consequently not exact; the \code{...} can be used
  to obtain a more accurate solution if necessary.

}
\seealso{
  \code{\link{fnormal1}},
  \code{\link[stats]{uniroot}}.

}
\examples{
\dontrun{
m <- 1.5; SD<-exp(0)
x <- seq(-1, 4, len = 501)
plot(x, dfnorm(x, m = m, sd = SD), type = "l", ylim = 0:1, las = 1,
     ylab = paste("fnorm(m = ", m, ", sd = ", round(SD, dig = 3), ")"),
     main = "Blue is density, red is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles", col = "blue")
lines(x, pfnorm(x, m = m, sd = SD), col = "red")
abline(h = 0)
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qfnorm(probs, m = m, sd = SD)
lines(Q, dfnorm(Q, m = m, sd = SD), col = "purple", lty = 3, type = "h")
lines(Q, pfnorm(Q, m = m, sd = SD), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(pfnorm(Q, m = m, sd = SD) - probs)) # Should be 0
}
}
\keyword{distribution}


