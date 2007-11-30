\name{Kumar}
\alias{Kumar}
\alias{dkumar}
\alias{pkumar}
\alias{qkumar}
\alias{rkumar}
\title{The Kumaraswamy Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Kumaraswamy distribution.

}
\usage{
dkumar(x, shape1, shape2)
pkumar(q, shape1, shape2)
qkumar(p, shape1, shape2)
rkumar(n, shape1, shape2)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{shape1, shape2}{ positive shape parameters. }
}
\value{
  \code{dkumar} gives the density,
  \code{pkumar} gives the distribution function,
  \code{qkumar} gives the quantile function, and
  \code{rkumar} generates random deviates.
}
\author{ T. W. Yee }
\details{
  See \code{\link{kumar}}, the \pkg{VGAM} family function
  for estimating the parameters, 
  for the formula of the probability density function and other details.

}
%\note{
%}
\seealso{
  \code{\link{kumar}}.

}
\examples{
\dontrun{
# shape1 = shape2 = 0.5;
shape1 = 2; shape2 = 2;
nn = 201
x = seq(0.01, 0.99, len=nn)
plot(x, dkumar(x, shape1, shape2), type="l", las=1, ylim=c(0,1.5),
     ylab=paste("fkumar(shape1=", shape1, ", shape2=", shape2, ")"),
     col="blue",
     main="Blue is density, red is cumulative distribution function",
     sub="Purple lines are the 10,20,...,90 percentiles")
lines(x, pkumar(x, shape1, shape2), col="red")
abline(h=0)
probs = seq(0.1, 0.9, by=0.1)
Q = qkumar(probs, shape1, shape2)
lines(Q, dkumar(Q, shape1, shape2), col="purple", lty=3, type="h")
lines(Q, pkumar(Q, shape1, shape2), col="purple", lty=3, type="h")
abline(h=probs, col="purple", lty=3)
max(abs(pkumar(Q, shape1, shape2) - probs)) # Should be 0
}
}
\keyword{distribution}


