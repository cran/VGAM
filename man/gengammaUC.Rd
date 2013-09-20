\name{gengammaUC}
\alias{gengammaUC}
\alias{dgengamma}
\alias{pgengamma}
\alias{qgengamma}
\alias{rgengamma}
\title{The Generalized Gamma Distribution }
\description{
  Density, distribution function, quantile function and random
  generation for the generalized  gamma distribution with
  scale parameter \code{scale},
  and parameters \code{d} and \code{k}.


}
\usage{
dgengamma(x, scale = 1, d = 1, k = 1, log = FALSE)
pgengamma(q, scale = 1, d = 1, k = 1)
qgengamma(p, scale = 1, d = 1, k = 1)
rgengamma(n, scale = 1, d = 1, k = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Positive integer of length 1.}
  \item{scale}{the (positive) scale parameter \eqn{b}.}
  \item{d, k}{the (positive) parameters \eqn{d} and \eqn{k}.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }

}
\value{
  \code{dgengamma} gives the density,
  \code{pgengamma} gives the distribution function,
  \code{qgengamma} gives the quantile function, and
  \code{rgengamma} generates random deviates.


}
\references{

Stacy, E. W. and Mihram, G. A. (1965)
Parameter estimation for a generalized gamma distribution.
\emph{Technometrics}, \bold{7}, 349--358.


}
\author{ T. W. Yee }
\details{
  See \code{\link{gengamma}}, the \pkg{VGAM} family function
  for estimating the generalized gamma distribution
  by maximum likelihood estimation,
  for formulae and other details.
  Apart from \code{n}, all the above arguments may be vectors and
  are recyled to the appropriate length if necessary.


}
% \note{
% }
\seealso{
  \code{\link{gengamma}}.


}
\examples{
\dontrun{ x <- seq(0, 14, by = 0.01); d <- 1.5; Scale <- 2; k <- 6
plot(x, dgengamma(x, Scale, d, k), type = "l", col = "blue", ylim = 0:1,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple are 5,10,...,95 percentiles", las = 1, ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(qgengamma(seq(0.05,0.95,by = 0.05), Scale, d, k),
      dgengamma(qgengamma(seq(0.05,0.95,by = 0.05), Scale, d, k),
                Scale, d, k), col = "purple", lty = 3, type = "h")
lines(x, pgengamma(x, Scale, d, k), type = "l", col = "orange")
abline(h = 0, lty = 2) }
}
\keyword{distribution}


