\name{Slash}
\alias{Slash}
\alias{dslash}
\alias{pslash}
\alias{rslash}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Slash Distribution }
\description{
  Density function, distribution function, and
  random generation for the slash distribution.

}
\usage{
dslash(x, mu=0, sigma=1, log=FALSE, smallno=.Machine$double.eps*1000)
pslash(q, mu=0, sigma=1)
rslash(n, mu=0, sigma=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{n}{number of observations. Must be a single positive integer. }
  \item{mu, sigma}{the mean and standard deviation of 
  the univariate normal distribution.
  
  }
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm of the density is returned.

  }
  \item{smallno}{
  See \code{\link{slash}}.

  }
}
\details{
  See \code{\link{slash}}, the \pkg{VGAM} family function
  for estimating the two parameters by maximum likelihood estimation,
  for the formula of the probability density function and other details.

}
\value{
  \code{dslash} gives the density, and
  \code{pslash} gives the distribution function,
  \code{rslash} generates random deviates.

}
%\references{ }
\author{ Thomas W. Yee and C. S. Chee}
\note{
  \code{pslash} is very slow.

}

\seealso{ 
    \code{\link{slash}}.

}
\examples{
\dontrun{
x = seq(-5,5,length=201)
plot(x, dslash(x), type="l", col="blue", ylab="f(x)", ylim=c(0,0.4),
     main="Standard Slash, Normal and Cauchy Densities", lwd=2)
lines(x, dnorm(x), col="black", lty=2, lwd=2)
lines(x, dcauchy(x), col="red", lty=3, lwd=2)
legend(x=2, y=0.3, c("slash","normal","cauchy"), lty=1:3,
       col=c("blue","black","red"), lwd=2)

plot(x, pslash(q=x), col="blue", type="l", ylim=0:1)
}
}
\keyword{distribution}
