\name{Truncpareto}
\alias{Truncpareto}
\alias{dtruncpareto}
\alias{ptruncpareto}
\alias{qtruncpareto}
\alias{rtruncpareto}
\title{The Truncated Pareto Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the upper truncated Pareto(I) distribution with parameters
  \code{lower}, \code{upper} and \code{shape}.


}
\usage{
dtruncpareto(x, lower, upper, shape, log = FALSE)
ptruncpareto(q, lower, upper, shape, lower.tail = TRUE, log.p = FALSE)
qtruncpareto(p, lower, upper, shape)
rtruncpareto(n, lower, upper, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n, log}{Same meaning as \code{\link[stats:Uniform]{runif}}.


  }
  \item{lower, upper, shape}{
  the lower, upper and shape (\eqn{k}) parameters.
  If necessary, values are recycled.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }

}
\value{
  \code{dtruncpareto} gives the density,
  \code{ptruncpareto} gives the distribution function,
  \code{qtruncpareto} gives the quantile function, and
  \code{rtruncpareto} generates random deviates.


}
\references{

  Aban, I. B., Meerschaert, M. M. and Panorska, A. K. (2006)
  Parameter estimation for the truncated Pareto distribution,
  \emph{Journal of the American Statistical Association},
  \bold{101}(473),
  270--277.



}
\author{ T. W. Yee and Kai Huang }
\details{

  See \code{\link{truncpareto}}, the \pkg{VGAM} family function
  for estimating the parameter \eqn{k} by maximum likelihood estimation,
  for the formula of the probability density function and the
  range restrictions imposed on the parameters.



}
%%\note{
%%  The truncated Pareto distribution is
%%}
\seealso{
  \code{\link{truncpareto}}.


}
\examples{ lower <- 3; upper <- 8; kay <- exp(0.5)
\dontrun{ xx <- seq(lower - 0.5, upper + 0.5, len = 401)
plot(xx, dtruncpareto(xx, low = lower, upp = upper, shape = kay),
     main = "Truncated Pareto density split into 10 equal areas",
     type = "l", ylim = 0:1, xlab = "x")
abline(h = 0, col = "blue", lty = 2)
qq <- qtruncpareto(seq(0.1, 0.9, by = 0.1), low = lower, upp = upper,
                   shape = kay)
lines(qq, dtruncpareto(qq, low = lower, upp = upper, shape = kay),
      col = "purple", lty = 3, type = "h")
lines(xx, ptruncpareto(xx, low = lower, upp = upper, shape = kay),
      col = "orange") }
pp <- seq(0.1, 0.9, by = 0.1)
qq <- qtruncpareto(pp, lower = lower, upper = upper, shape = kay)

ptruncpareto(qq, lower = lower, upper = upper, shape = kay)
qtruncpareto(ptruncpareto(qq, lower = lower, upper = upper, shape = kay),
         lower = lower, upper = upper, shape = kay) - qq  # Should be all 0
}
\keyword{distribution}



