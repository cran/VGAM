\name{Foldnorm}
\alias{Foldnorm}
\alias{dfoldnorm}
\alias{pfoldnorm}
\alias{qfoldnorm}
\alias{rfoldnorm}
\title{The Folded-Normal Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the (generalized) folded-normal distribution.


}
\usage{
dfoldnorm(x, mean = 0, sd = 1, a1 = 1, a2 = 1, log = FALSE)
pfoldnorm(q, mean = 0, sd = 1, a1 = 1, a2 = 1,
          lower.tail = TRUE, log.p = FALSE)
qfoldnorm(p, mean = 0, sd = 1, a1 = 1, a2 = 1,
          lower.tail = TRUE, log.p = FALSE, ...)
rfoldnorm(n, mean = 0, sd = 1, a1 = 1, a2 = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Same as \code{\link[stats:Normal]{rnorm}}.


  }
  \item{mean, sd}{ see \code{\link[stats:Normal]{rnorm}}. }
  \item{a1, a2}{ see \code{\link{foldnormal}}. }
  \item{log}{
  Logical.
  If \code{TRUE} then the log density is returned.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }
  \item{\ldots}{
  Arguments that can be passed into \code{\link[stats]{uniroot}}.


  }
}
\value{
  \code{dfoldnorm} gives the density,
  \code{pfoldnorm} gives the distribution function,
  \code{qfoldnorm} gives the quantile function, and
  \code{rfoldnorm} generates random deviates.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{foldnormal}}, the \pkg{VGAM} family function
  for estimating the parameters,
  for the formula of the probability density function and other details.


}
\note{
  \code{qfoldnorm()} runs very slowly because it calls
  \code{\link[stats]{uniroot}} for each value of the argument \code{p}.
  The solution is consequently not exact; the \code{...} can be used
  to obtain a more accurate solution if necessary.


}
\seealso{
  \code{\link{foldnormal}},
  \code{\link[stats]{uniroot}}.


}
\examples{
\dontrun{
m <- 1.5; SD <- exp(0)
x <- seq(-1, 4, len = 501)
plot(x, dfoldnorm(x, m = m, sd = SD), type = "l", ylim = 0:1, las = 1,
     ylab = paste("foldnorm(m = ", m, ", sd = ", round(SD, digits = 3), ")"),
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles", col = "blue")
abline(h = 0, col = "gray50")
lines(x, pfoldnorm(x, m = m, sd = SD), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qfoldnorm(probs, m = m, sd = SD)
lines(Q, dfoldnorm(Q, m = m, sd = SD), col = "purple", lty = 3, type = "h")
lines(Q, pfoldnorm(Q, m = m, sd = SD), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(pfoldnorm(Q, m = m, sd = SD) - probs))  # Should be 0
}
}
\keyword{distribution}


