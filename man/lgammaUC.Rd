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
dlgamma(x, location = 0, scale = 1, shape = 1, log = FALSE)
plgamma(q, location = 0, scale = 1, shape = 1,
        lower.tail = TRUE, log.p = FALSE)
qlgamma(p, location = 0, scale = 1, shape = 1,
        lower.tail = TRUE, log.p = FALSE)
rlgamma(n, location = 0, scale = 1, shape = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
   Same as \code{\link[stats]{runif}}.


  }
  \item{location}{the location parameter \eqn{a}.}
  \item{scale}{the (positive) scale parameter \eqn{b}.}
  \item{shape}{the (positive) shape parameter \eqn{k}.}
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
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{lgamma1}}, the \pkg{VGAM} family function for
  estimating the one parameter standard log-gamma distribution by maximum
  likelihood estimation, for formulae and other details.  Apart from
  \code{n}, all the above arguments may be vectors and are recyled to
  the appropriate length if necessary.



}
\note{
  The \pkg{VGAM} family function \code{\link{lgamma3}} is
  for the three parameter (nonstandard) log-gamma distribution.



}
\seealso{
  \code{\link{lgamma1}},
  \code{\link{prentice74}}.


}
\examples{
\dontrun{ loc <- 1; Scale <- 1.5; shape <- 1.4
x <- seq(-3.2, 5, by = 0.01)
plot(x, dlgamma(x, loc = loc, Scale, shape = shape), type = "l",
     col = "blue", ylim = 0:1,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple are 5,10,...,95 percentiles", las = 1, ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(qlgamma(seq(0.05, 0.95, by = 0.05), loc = loc, Scale, shape = shape),
      dlgamma(qlgamma(seq(0.05, 0.95, by = 0.05), loc = loc, scale = Scale,
                      shape = shape),
      loc = loc, Scale, shape = shape), col = "purple", lty = 3, type = "h")
lines(x, plgamma(x, loc = loc, Scale, shape = shape), col = "orange")
abline(h = 0, lty = 2) }
}
\keyword{distribution}


