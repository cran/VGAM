\name{Pareto}
\alias{Pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
\title{The Pareto Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Pareto(I) distribution with parameters
  \code{location} and \code{shape}.

}
\usage{
dpareto(x, location, shape)
ppareto(q, location, shape)
qpareto(p, location, shape)
rpareto(n, location, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Must be a single positive integer. }
  \item{location, shape}{the \eqn{\alpha}{alpha} and \eqn{k} parameters.}
}
\value{
  \code{dpareto} gives the density,
  \code{ppareto} gives the distribution function,
  \code{qpareto} gives the quantile function, and
  \code{rpareto} generates random deviates.
}
\references{
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.
}
\author{ T. W. Yee }
\details{
  See \code{\link{pareto1}}, the \pkg{VGAM} family function
  for estimating the parameter \eqn{k} by maximum likelihood estimation,
  for the formula of the probability density function and the
  range restrictions imposed on the parameters.
}
%%\note{
%%  The Pareto distribution is 
%%}
\seealso{
  \code{\link{pareto1}},
  \code{\link{ParetoIV}}.

}
\examples{
alpha = 3; k = exp(1); x = seq(2.8, 8, len=300)
\dontrun{
plot(x, dpareto(x, location=alpha, shape=k), type="l",
     main="Pareto density split into 10 equal areas")
abline(h=0, col="blue", lty=2)
qq = qpareto(seq(0.1,0.9,by=0.1),location=alpha,shape=k)
lines(qq, dpareto(qq, loc=alpha, shape=k), col="purple", lty=3, type="h")
}
pp = seq(0.1,0.9,by=0.1)
qq = qpareto(pp, location=alpha, shape=k)
ppareto(qq, location=alpha, shape=k)
qpareto(ppareto(qq,loc=alpha,shape=k),loc=alpha,shape=k) - qq # Should be 0
}
\keyword{distribution}


