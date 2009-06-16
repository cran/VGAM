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
dbenini(x, shape, y0, log=FALSE)
pbenini(q, shape, y0)
qbenini(p, shape, y0)
rbenini(n, shape, y0)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{shape}{the shape parameter \eqn{b}.}
  \item{y0}{the scale parameter \eqn{y_0}{y0}.}
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.

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
Hoboken, NJ: Wiley-Interscience.

}
\author{ T. W. Yee }
\details{
  See \code{\link{benini}}, the \pkg{VGAM} family function
  for estimating the parameter \eqn{b} by maximum likelihood estimation,
  for the formula of the probability density function and other details.
}
%\note{
%  
%}
\seealso{
  \code{\link{benini}}.
}
\examples{
\dontrun{
y0 = 1
shape = exp(1)
x = seq(-0.0, 4, len=101)
plot(x, dbenini(x, y0=y0,shape=shape), type="l", col="blue", las=1,
     main="Blue is density, red is cumulative distribution function",
     sub="Purple lines are the 10,20,...,90 percentiles", ylim=0:1,
     ylab="")
abline(h=0, col="blue", lty=2)
lines(x, pbenini(x, y0=y0,shape=shape), col="red")
probs = seq(0.1, 0.9, by=0.1)
Q = qbenini(probs, y0=y0,shape=shape)
lines(Q, dbenini(Q, y0=y0,shape=shape), col="purple", lty=3, type="h")
pbenini(Q, y0=y0,shape=shape) - probs    # Should be all zero
}
}
\keyword{distribution}


