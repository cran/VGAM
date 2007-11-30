\name{Yules}
\alias{Yules}
\alias{dyules}
\alias{pyules}
%\alias{qyules}
\alias{ryules}
\title{ Yule-Simon Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the Yule-Simon distribution.

}
\usage{
dyules(x, rho, log.arg=FALSE)
pyules(q, rho)
%qyules(p, rho)
ryules(n, rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{
   Vector of quantiles. For the density, it should be a vector with
   positive integer values in order for the probabilities to be positive.

  }
% \item{p}{vector of probabilities.}
  \item{n}{number of observations. A single positive integer.}
  \item{rho}{ 
   See \code{\link{yulesimon}}.

  }
  \item{log.arg}{logical; if TRUE, the logarithm is returned. }
}
\details{
   See \code{\link{yulesimon}}, the \pkg{VGAM} family function
  for estimating the parameter,
  for the formula of the probability density function and other details.

}
\value{
  \code{dyules} gives the density,
  \code{pyules} gives the distribution function, and
% \code{qyules} gives the quantile function, and
  \code{ryules} generates random deviates.

}
%\references{ 
%
%}
\author{ T. W. Yee }
%\note{ 
%}

\seealso{
  \code{\link{yulesimon}}.

}
\examples{
dyules(1:20, 2.1)
ryules(20, 2.1)

round(1000 * dyules(1:8, 2))
table(ryules(1000, 2))

\dontrun{
x = 0:6
plot(x, dyules(x, rho=2.2), type="h", las=1, col="blue")
}
}
\keyword{distribution}
