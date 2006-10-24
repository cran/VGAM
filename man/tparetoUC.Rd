\name{Tpareto}
\alias{Tpareto}
\alias{dtpareto}
\alias{ptpareto}
\alias{qtpareto}
\alias{rtpareto}
\title{The Truncated Pareto Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the upper truncated Pareto(I) distribution with parameters
  \code{lower}, \code{upper} and \code{shape}.

}
\usage{
dtpareto(x, lower, upper, shape)
ptpareto(q, lower, upper, shape)
qtpareto(p, lower, upper, shape)
rtpareto(n, lower, upper, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Must be a single positive integer. }
  \item{lower, upper, shape}{
  the lower, upper and shape (\eqn{k}) parameters.
  If necessary, values are recycled.
  }
}
\value{
  \code{dtpareto} gives the density,
  \code{ptpareto} gives the distribution function,
  \code{qtpareto} gives the quantile function, and
  \code{rtpareto} generates random deviates.
}
\references{
  Aban, I. B., Meerschaert, M. M. and Panorska, A. K. (2006).
  Parameter estimation for the truncated Pareto distribution,
  \emph{Journal of the American Statistical Association},
  \bold{101}(473),
  270--277.

}
\author{ T. W. Yee }
\details{
  See \code{\link{tpareto1}}, the \pkg{VGAM} family function
  for estimating the parameter \eqn{k} by maximum likelihood estimation,
  for the formula of the probability density function and the
  range restrictions imposed on the parameters.

}
%%\note{
%%  The truncated Pareto distribution is 
%%}
\seealso{
  \code{\link{tpareto1}}.

}
\examples{
lower = 3; upper = 8; k = exp(0.5)
\dontrun{
x = seq(lower, upper, len=200)
plot(x, dtpareto(x, lo=lower, up=upper, shape=k), type="l",
     main="Truncated Pareto density split into 10 equal areas")
abline(h=0, col="blue", lty=2)
qq = qtpareto(seq(0.1,0.9,by=0.1),lo=lower, up=upper,shape=k)
lines(qq, dtpareto(qq, lo=lower, up=upper, shape=k),
      col="purple", lty=3, type="h")
}
pp = seq(0.1,0.9,by=0.1)
qq = qtpareto(pp, lo=lower, up=upper, shape=k)
ptpareto(qq, lo=lower, up=upper, shape=k)
qtpareto(ptpareto(qq, lo=lower, up=upper, shape=k),
         lo=lower, up=upper, shape=k) - qq # Should be all 0
}
\keyword{distribution}


