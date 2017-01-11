\name{Posnorm}
\alias{Posnorm}
\alias{dposnorm}
\alias{pposnorm}
\alias{qposnorm}
\alias{rposnorm}
\title{The Positive-Normal Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the univariate positive-normal distribution.

}
\usage{
dposnorm(x, mean = 0, sd = 1, log = FALSE)
pposnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
qposnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
rposnorm(n, mean = 0, sd = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.

  }
  \item{mean, sd, log, lower.tail, log.p}{
  see \code{\link[stats:Normal]{rnorm}}.


  }

}
\value{
  \code{dposnorm} gives the density,
  \code{pposnorm} gives the distribution function,
  \code{qposnorm} gives the quantile function, and
  \code{rposnorm} generates random deviates.


}
\author{ T. W. Yee }
\details{
  See \code{\link{posnormal}}, the \pkg{VGAM} family function
  for estimating the parameters,
  for the formula of the probability density function and other details.


}
%\note{

%}
\seealso{
  \code{\link{posnormal}}.

}
\examples{
\dontrun{ m <-  0.8; x <- seq(-1, 4, len = 501)
plot(x, dposnorm(x, m = m), type = "l", las = 1, ylim = 0:1,
     ylab = paste("posnorm(m = ", m, ", sd = 1)"), col = "blue",
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles")
abline(h = 0, col = "grey")
lines(x, pposnorm(x, m = m), col = "orange", type = "l")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qposnorm(probs, m = m)
lines(Q, dposnorm(Q, m = m), col = "purple", lty = 3, type = "h")
lines(Q, pposnorm(Q, m = m), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(pposnorm(Q, m = m) - probs))  # Should be 0
}
}
\keyword{distribution}


% 20150207; bug involving ifelse() picked up for qposnorm().




