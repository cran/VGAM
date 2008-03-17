\name{alaplaceUC}
\alias{dalaplace}
\alias{palaplace}
\alias{qalaplace}
\alias{ralaplace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Laplace Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the 3-parameter asymmetric Laplace distribution with location
  parameter \code{location}, scale parameter \code{scale},
  and asymmetry parameter \code{kappa}.

}
\usage{
dalaplace(x, location=0, scale=1, tau=0.5, kappa=sqrt(tau/(1-tau)))
palaplace(q, location=0, scale=1, tau=0.5, kappa=sqrt(tau/(1-tau)))
qalaplace(p, location=0, scale=1, tau=0.5, kappa=sqrt(tau/(1-tau)))
ralaplace(n, location=0, scale=1, tau=0.5, kappa=sqrt(tau/(1-tau)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Positive integer of length 1.}
  \item{location}{
    the location parameter \eqn{\xi}{xi}.
  }
  \item{scale}{
  the scale parameter \eqn{b}.
  Must consist of positive values.
  }
  \item{tau}{
  the quantile parameter \eqn{\tau}{tau}.
  Must consist of values in \eqn{(0,1)}.
  This argument is used to specify \code{kappa} and is ignored
  if \code{kappa} is assigned.

  }
  \item{kappa}{
  the asymmetry parameter \eqn{\kappa}{kappa}.
  Must consist of positive values.
  }
}
\details{
  There are several variants of asymmetric Laplace distributions (ALDs) and
  this one is known as \emph{the} ALD by Kotz et al. (2001).
  See \code{\link{alaplace3}}, the \pkg{VGAM} family function
  for estimating the three parameters by maximum likelihood estimation,
  for formulae and details.
  Apart from \code{n}, all the above arguments may be vectors and
  are recyled to the appropriate length if necessary.

}
\value{
  \code{dalaplace} gives the density,
  \code{palaplace} gives the distribution function,
  \code{qalaplace} gives the quantile function, and
  \code{ralaplace} generates random deviates.
}
\references{
Kotz, S., Kozubowski, T. J. and Podgorski, K. (2001)
\emph{The Laplace distribution and generalizations:
a revisit with applications to communications,
economics, engineering, and finance},
Boston: Birkhauser.

}
\author{ T. W. Yee }
%\note{
%  The \pkg{VGAM} family function \code{\link{alaplace3}}
%  estimates the three parameters by maximum likelihood estimation.
%}

\seealso{
  \code{\link{alaplace3}}.
}
\examples{
x = seq(-5, 5, by=0.01)
loc = 0; b = 1.5; kappa = 2
\dontrun{
plot(x, dalaplace(x, loc, b, kappa=kappa), type="l", col="blue",
     main="Blue is density, red is cumulative distribution function",
     ylim=c(0,1), sub="Purple are 5,10,...,95 percentiles", las=1, ylab="")
abline(h=0, col="blue", lty=2)
lines(qalaplace(seq(0.05,0.95,by=0.05), loc, b, kappa=kappa),
      dalaplace(qalaplace(seq(0.05,0.95,by=0.05), loc, b, kappa=kappa),
                loc, b, kappa=kappa), col="purple", lty=3, type="h")
lines(x, palaplace(x, loc, b, kappa=kappa), type="l", col="red")
abline(h=0, lty=2)
}
palaplace(qalaplace(seq(0.05,0.95,by=0.05), loc, b, kappa=kappa),
          loc, b, kappa=kappa)
}
\keyword{distribution}

