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
dtpareto(x, lower, upper, shape, log = FALSE)
ptpareto(q, lower, upper, shape)
qtpareto(p, lower, upper, shape)
rtpareto(n, lower, upper, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n, log}{Same meaning as \code{\link[stats:Uniform]{runif}}. }
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
  Aban, I. B., Meerschaert, M. M. and Panorska, A. K. (2006)
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
\examples{ lower <- 3; upper <- 8; kay <- exp(0.5)
\dontrun{ xx <- seq(lower - 0.5, upper + 0.5, len = 401)
plot(xx, dtpareto(xx, low = lower, upp = upper, shape = kay),
     main = "Truncated Pareto density split into 10 equal areas",
     type = "l", ylim = 0:1, xlab = "x")
abline(h = 0, col = "blue", lty = 2)
qq <- qtpareto(seq(0.1, 0.9, by = 0.1), low = lower, upp = upper,
               shape = kay)
lines(qq, dtpareto(qq, low = lower, upp = upper, shape = kay),
      col = "purple", lty = 3, type = "h")
lines(xx, ptpareto(xx, low = lower, upp = upper, shape = kay),
      col = "orange") }
pp <- seq(0.1, 0.9, by = 0.1)
qq <- qtpareto(pp, low = lower, upp = upper, shape = kay)

ptpareto(qq, low = lower, upp = upper, shape = kay)
qtpareto(ptpareto(qq, low = lower, upp = upper, shape = kay),
         low = lower, upp = upper, shape = kay) - qq # Should be all 0
}
\keyword{distribution}



