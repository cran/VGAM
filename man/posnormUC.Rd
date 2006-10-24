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
dposnorm(x, mean=0, sd=1)
pposnorm(q, mean=0, sd=1)
qposnorm(p, mean=0, sd=1)
rposnorm(n, mean=0, sd=1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{mean, sd}{ see \code{\link[stats]{rnorm}}. }
}
\value{
  \code{dposnorm} gives the density,
  \code{pposnorm} gives the distribution function,
  \code{qposnorm} gives the quantile function, and
  \code{rposnorm} generates random deviates.
}
\author{ T. W. Yee }
\details{
  See \code{\link{posnormal1}}, the \pkg{VGAM} family function
  for estimating the parameters, 
  for the formula of the probability density function and other details.
}
\note{
  \code{rposnorm()} may run very slowly if the mean is very negative.
}
\seealso{
  \code{\link{posnormal1}}.
}
\examples{
\dontrun{
m =  0.8
x = seq(-1, 4, len=501)
plot(x, dposnorm(x, m=m), type="l", ylim=0:1, las=1,
     ylab=paste("posnorm(m=",m,", sd=1)"), col="blue",
     main="Blue is density, red is cumulative distribution function",
     sub="Purple lines are the 10,20,...,90 percentiles")
lines(x, pposnorm(x, m=m), col="red")
abline(h=0)
probs = seq(0.1, 0.9, by=0.1)
Q = qposnorm(probs, m=m)
lines(Q, dposnorm(Q, m=m), col="purple", lty=3, type="h")
lines(Q, pposnorm(Q, m=m), col="purple", lty=3, type="h")
abline(h=probs, col="purple", lty=3)
pposnorm(Q, m=m) - probs # Should be all 0
}
}
\keyword{distribution}


