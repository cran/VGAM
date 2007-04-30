\name{lgammaUC}
\alias{lgammaUC}
\alias{dlgamma}
\alias{plgamma}
\alias{qlgamma}
\alias{rlgamma}
\title{The Log-Gamma Distribution }
\description{
  Density, distribution function, quantile function and random
  generation for the log-gamma distribution with
  location parameter \code{location},
  scale parameter \code{scale} and
  shape parameter \code{k}.

}
\usage{
dlgamma(x, location=0, scale=1, k=1)
plgamma(q, location=0, scale=1, k=1)
qlgamma(p, location=0, scale=1, k=1)
rlgamma(n, location=0, scale=1, k=1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Positive integer of length 1.}
  \item{location}{the location parameter \eqn{a}.}
  \item{scale}{the (positive) scale parameter \eqn{b}.}
  \item{k}{the (positive) shape parameter \eqn{k}.}
}
\value{
  \code{dlgamma} gives the density,
  \code{plgamma} gives the distribution function,
  \code{qlgamma} gives the quantile function, and
  \code{rlgamma} generates random deviates.
}
\references{
Kotz, S. and Nadarajah, S. (2000)
\emph{Extreme Value Distributions: Theory and Applications},
pages 48--49,
London: Imperial College Press.

}
\author{ T. W. Yee }
\details{
  See \code{\link{lgammaff}}, the \pkg{VGAM} family function for
  estimating the one parameter standard log-gamma distribution by maximum
  likelihood estimation, for formulae and other details.  Apart from
  \code{n}, all the above arguments may be vectors and are recyled to
  the appropriate length if necessary.

}
\note{
  The \pkg{VGAM} family function \code{\link{lgamma3ff}} is
  for the three parameter (nonstandard) log-gamma distribution.
}
\seealso{
  \code{\link{lgammaff}},
  \code{\link{prentice74}}.
}
\examples{
\dontrun{
x = seq(-3.2, 5, by=0.01)
loc = 1
Scale = 1.5
k = 1.4
plot(x, dlgamma(x, loc, Scale, k), type="l", col="blue", ylim=c(0,1),
     main="Blue is density, red is cumulative distribution function",
     sub="Purple are 5,10,...,95 percentiles", las=1, ylab="")
abline(h=0, col="blue", lty=2)
lines(qlgamma(seq(0.05,0.95,by=0.05), loc, Scale, k), 
      dlgamma(qlgamma(seq(0.05,0.95,by=0.05), loc, Scale, k), loc, Scale, k),
      col="purple", lty=3, type="h")
lines(x, plgamma(x, loc, Scale, k), type="l", col="red")
abline(h=0, lty=2)
}
}
\keyword{distribution}

