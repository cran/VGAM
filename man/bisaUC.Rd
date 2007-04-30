\name{Bisa}
\alias{Bisa}
\alias{dbisa}
\alias{pbisa}
\alias{qbisa}
\alias{rbisa}
\title{The Birnbaum-Saunders Distribution}
\description{
  Density, distribution function, and random
  generation for the Birnbaum-Saunders distribution.

}
\usage{
dbisa(x, shape, scale=1, log=FALSE)
pbisa(q, shape, scale=1)
qbisa(p, shape, scale=1)
rbisa(n, shape, scale=1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Must be a positive integer of length 1.}
  \item{shape, scale}{
  the (positive) shape and scale parameters.

  }
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm of the density is returned.

  }
}
\value{
  \code{dbisa} gives the density,
  \code{pbisa} gives the distribution function, and
  \code{qbisa} gives the quantile function, and
  \code{rbisa} generates random deviates.
}
\author{ T. W. Yee }
\details{
  The Birnbaum-Saunders distribution
  is a distribution which is used in survival analysis.
  See \code{\link{bisa}}, the \pkg{VGAM} family function
  for estimating the parameters, 
  for more details.

}
%\note{
%}
\seealso{
  \code{\link{bisa}}.
}
\examples{
\dontrun{
x = seq(0, 6, len=400)
plot(x, dbisa(x, shape=1), type="l", col="blue", ylab="Density", lwd=2,
     main="X ~ Birnbaum-Saunders(shape, scale=1)", ylim=c(0,1.3), lty=3)
lines(x, dbisa(x, shape=2), col="red", lty=2, lwd=2)
lines(x, dbisa(x, shape=0.5), col="green", lty=1, lwd=2)
legend(x=3, y=0.9, legend=paste("shape =",c(0.5,1,2)),
       col=c("green","blue","red"), lty=1:3, lwd=2)

shape = 1
x = seq(0.0, 4, len=401)
plot(x, dbisa(x, shape=shape), type="l", col="blue", las=1, ylab="",
     main="Blue is density, red is cumulative distribution function",
     sub="Purple lines are the 10,20,...,90 percentiles", ylim=0:1)
abline(h=0, col="blue", lty=2)
lines(x, pbisa(x, shape=shape), col="red")
probs = seq(0.1, 0.9, by=0.1)
Q = qbisa(probs, shape=shape)
lines(Q, dbisa(Q, shape=shape), col="purple", lty=3, type="h")
pbisa(Q, shape=shape) - probs    # Should be all zero
abline(h=probs, col="purple", lty=3)
lines(Q, pbisa(Q, shape), col="purple", lty=3, type="h")
}
}
\keyword{distribution}

