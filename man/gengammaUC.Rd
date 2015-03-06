\name{gengammaUC}
\alias{gengammaUC}
\alias{dgengamma.stacy}
\alias{pgengamma.stacy}
\alias{qgengamma.stacy}
\alias{rgengamma.stacy}
\title{The Generalized Gamma Distribution }
\description{
  Density, distribution function, quantile function and random
  generation for the generalized  gamma distribution with
  scale parameter \code{scale},
  and parameters \code{d} and \code{k}.


}
\usage{
dgengamma.stacy(x, scale = 1, d = 1, k = 1, log = FALSE)
pgengamma.stacy(q, scale = 1, d = 1, k = 1,
                lower.tail = TRUE, log.p = FALSE)
qgengamma.stacy(p, scale = 1, d = 1, k = 1,
                lower.tail = TRUE, log.p = FALSE)
rgengamma.stacy(n, scale = 1, d = 1, k = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as in \code{\link[stats]{runif}}.


  }
  \item{scale}{the (positive) scale parameter \eqn{b}.}
  \item{d, k}{the (positive) parameters \eqn{d} and \eqn{k}.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }

}
\value{
  \code{dgengamma.stacy} gives the density,
  \code{pgengamma.stacy} gives the distribution function,
  \code{qgengamma.stacy} gives the quantile function, and
  \code{rgengamma.stacy} generates random deviates.


}
\references{

Stacy, E. W. and Mihram, G. A. (1965)
Parameter estimation for a generalized gamma distribution.
\emph{Technometrics}, \bold{7}, 349--358.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{gengamma.stacy}}, the \pkg{VGAM} family function
  for estimating the generalized gamma distribution
  by maximum likelihood estimation,
  for formulae and other details.
  Apart from \code{n}, all the above arguments may be vectors and
  are recyled to the appropriate length if necessary.


}
% \note{
% }
\seealso{
  \code{\link{gengamma.stacy}}.


}
\examples{
\dontrun{ x <- seq(0, 14, by = 0.01); d <- 1.5; Scale <- 2; k <- 6
plot(x, dgengamma.stacy(x, Scale, d = d, k = k), type = "l", col = "blue", ylim = 0:1,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple are 5,10,...,95 percentiles", las = 1, ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(qgengamma.stacy(seq(0.05, 0.95, by = 0.05), Scale, d = d, k = k),
      dgengamma.stacy(qgengamma.stacy(seq(0.05,0.95,by = 0.05), Scale, d = d, k = k),
                      Scale, d = d, k = k), col = "purple", lty = 3, type = "h")
lines(x, pgengamma.stacy(x, Scale, d = d, k = k), type = "l", col = "orange")
abline(h = 0, lty = 2) }
}
\keyword{distribution}


