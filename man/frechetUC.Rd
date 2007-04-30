\name{Frechet}
\alias{Frechet}
\alias{dfrechet}
\alias{pfrechet}
\alias{qfrechet}
\alias{rfrechet}
\title{The Frechet Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the three parameter Frechet distribution.
}
\usage{
dfrechet(x, location=0, scale=1, shape)
pfrechet(q, location=0, scale=1, shape)
qfrechet(p, location=0, scale=1, shape)
rfrechet(n, location=0, scale=1, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{location, scale, shape}{the location parameter \eqn{a},
  scale parameter \eqn{b}, and shape parameter \eqn{s}.}
}
\value{
  \code{dfrechet} gives the density,
  \code{pfrechet} gives the distribution function,
  \code{qfrechet} gives the quantile function, and
  \code{rfrechet} generates random deviates.
}
\references{
Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
\emph{Extreme Value and Related Models with Applications in Engineering and Science},
Hoboken, N.J.: Wiley-Interscience.

}
\author{ T. W. Yee }
\details{
  See \code{\link{frechet2}} and \code{\link{frechet3}}, the \pkg{VGAM}
  family functions for estimating the 2 (without location parameter) and
  3 parameters by maximum likelihood estimation, for the formula of the
  probability density function and range restrictions on the parameters.

}
%\note{
%}
\seealso{
  \code{\link{frechet2}},
  \code{\link{frechet3}}.
}
\examples{
\dontrun{
shape = 5
x = seq(-0.1, 3.5, len=100)
plot(x, dfrechet(x, shape=shape), type="l", ylab="", las=1,
     main="Frechet density divided into 10 equal areas; red=cdf")
abline(h=0, col="blue", lty=2)
qq = qfrechet(seq(0.1,0.9,by=0.1), shape=shape)
lines(qq, dfrechet(qq, shape=shape), col="purple", lty=3, type="h")
lines(x, pfrechet(q=x, shape=shape), col="red")
}
}
\keyword{distribution}

