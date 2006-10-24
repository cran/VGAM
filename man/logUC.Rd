\name{Log}
\alias{Log}
\alias{dlog}
\alias{plog}
\alias{qlog}
\alias{rlog}
\title{ Logarithmic Distribution }
\description{
  Density, distribution function,
% quantile function
  and random generation
  for the logarithmic distribution.

}
\usage{
dlog(x, prob)
plog(q, prob, log.p=FALSE)
% qlog(p, prob)
rlog(n, prob, Smallno=1.0e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{
   Vector of quantiles. For the density, it should be a vector with
   positive integer values in order for the probabilities to be positive.

  }
% \item{p}{vector of probabilities.}
  \item{n}{number of observations. A single positive integer.}
  \item{prob}{ 
   The parameter value \eqn{c} described in in \code{\link{logff}}.
   Here it is called \code{prob} because \eqn{0<c<1} is the range.
   For \code{rlog()} this parameter must be of length 1.

  }
  \item{log.p}{
  Logical.
  If \code{TRUE} then all probabilities \code{p} are given as \code{log(p)}.

  }
  \item{Smallno}{
  Numeric, a small value used by the rejection method for determining
  the upper limit of the distribution.
  That is, \code{plog(U, prob) > 1-Smallno} where \code{U} is the upper limit.

  }
}
\details{
   The details are given in \code{\link{logff}}.

}
\value{
  \code{dlog} gives the density,
  \code{plog} gives the distribution function, and
% \code{qlog} gives the quantile function, and
  \code{rlog} generates random deviates.

}
\references{ 

Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.

}
\author{ T. W. Yee }
\note{ 
   Given some response data, the \pkg{VGAM} family function
   \code{\link{logff}} estimates the parameter \code{prob}.
}

\seealso{
    \code{\link{logff}}.
}
\examples{
dlog(1:20, 0.5)
rlog(20, 0.5)

\dontrun{
prob = 0.8; x = 1:10
plot(x, dlog(x, prob=prob), type="h", ylim=0:1,
     sub="prob=0.8", las=1, col="blue", ylab="Probability",
     main="Logarithmic distribution: blue=density; red=distribution function")
lines(x+0.1, plog(x, prob=prob), col="red", lty=3, type="h")
}
}
\keyword{distribution}
