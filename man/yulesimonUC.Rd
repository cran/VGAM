\name{Yules}
\alias{Yules}
\alias{dyules}
\alias{pyules}
\alias{qyules}
\alias{ryules}
\title{ Yule-Simon Distribution }
\description{
  Density, distribution function, quantile function and random
  generation for the Yule-Simon distribution.

}
\usage{
dyules(x, shape, log = FALSE)
pyules(q, shape, lower.tail = TRUE, log.p = FALSE)
qyules(p, shape)
ryules(n, shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
  Same meaning as in \code{\link[stats]{Normal}}.
%  Vector of quantiles. For the density, it should be a vector
%  with positive integer values in order for the probabilities
%  to be positive.


  }
% \item{p}{vector of probabilities.}
% \item{n}{number of observations. Same as in \code{\link[stats]{runif}}. }
  \item{shape}{
   See \code{\link{yulesimon}}.


  }
% \item{log}{logical; if TRUE, the logarithm is returned. }
  \item{log, lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }
}
\details{
  See \code{\link{yulesimon}}, the \pkg{VGAM} family function
  for estimating the parameter,
  for the formula of the probability density function and other details.


}
\value{
  \code{dyules} gives the density,
  \code{pyules} gives the distribution function,
  \code{qyules} gives the quantile function, and
  \code{ryules} generates random deviates.




}
%\references{
%
%}
\author{ T. W. Yee }
\note{
  Numerical problems may occur with
  \code{qyules()} when \code{p} is very close to 1.



}

\seealso{
  \code{\link{yulesimon}}.


}
\examples{
dyules(1:20, 2.1)
ryules(20, 2.1)

round(1000 * dyules(1:8, 2))
table(ryules(1000, 2))

\dontrun{ x <- 0:6
plot(x, dyules(x, shape = 2.2), type = "h", las = 1, col = "blue") }
}
\keyword{distribution}
