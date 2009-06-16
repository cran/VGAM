\name{Triangle}
\alias{Triangle}
\alias{dtriangle}
\alias{ptriangle}
\alias{qtriangle}
\alias{rtriangle}
\title{The Triangle Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Triangle distribution with parameter
  \code{theta}.
}
\usage{
dtriangle(x, theta, lower=0, upper=1, log=FALSE)
ptriangle(q, theta, lower=0, upper=1)
qtriangle(p, theta, lower=0, upper=1)
rtriangle(n, theta, lower=0, upper=1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{theta}{the theta parameter which lies between \code{lower}
     and \code{upper}. }
  \item{lower, upper}{lower and upper limits of the distribution.
     Must be finite.
   }
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dtriangle} gives the density,
  \code{ptriangle} gives the distribution function,
  \code{qtriangle} gives the quantile function, and
  \code{rtriangle} generates random deviates.
}
%\references{
%
%}
\author{ T. W. Yee }
\details{
  See \code{\link{triangle}}, the \pkg{VGAM} family function
  for estimating the parameter \eqn{\theta}{theta} by
  maximum likelihood estimation.

}
%\note{
%  
%}
\seealso{
  \code{\link{triangle}}.
}
\examples{
\dontrun{
x = seq(-0.1, 1.1, by=0.01)
theta = 0.75
plot(x, dtriangle(x, theta=theta), type="l", col="blue", las=1,
     main="Blue is density, red is cumulative distribution function",
     sub="Purple lines are the 10,20,...,90 percentiles",
     ylim=c(0,2), ylab="")
abline(h=0, col="blue", lty=2)
lines(x, ptriangle(x, theta=theta), col="red")
probs = seq(0.1, 0.9, by=0.1)
Q = qtriangle(probs, theta=theta)
lines(Q, dtriangle(Q, theta=theta), col="purple", lty=3, type="h")
ptriangle(Q, theta=theta) - probs    # Should be all zero
abline(h=probs, col="purple", lty=3)
}
}
\keyword{distribution}


