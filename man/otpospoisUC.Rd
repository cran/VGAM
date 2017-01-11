\name{Otpospois}
\alias{Otpospois}
\alias{dotpospois}
\alias{potpospois}
\alias{qotpospois}
\alias{rotpospois}
\title{ One-truncated Positive-Poisson Distribution }
\description{
  Density, distribution function,
  quantile function,
  and random generation
  for the one-truncated positive-Poisson distribution.



}
\usage{
dotpospois(x, lambda, log = FALSE)
potpospois(q, lambda, log.p = FALSE)
qotpospois(p, lambda)
rotpospois(n, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
   Same as \code{\link{Pospois}}.


  }
  \item{lambda, log, log.p}{
   Same as \code{\link{Pospois}}.



  }
}
\details{
  The one-truncated positive-Poisson is a Poisson distribution but with
  the probability of a one and a zero being 0.
  That is, its support is 2, 3, \ldots.
  The other probabilities are scaled to add to unity.
  Some more details are given in \code{\link{pospoisson}}.


}
\value{
  \code{dotpospois} gives the density,
  \code{potpospois} gives the distribution function,
  \code{qotpospois} gives the quantile function, and
  \code{rotpospois} generates random deviates.




}
%\references{
%}
\author{ T. W. Yee }
\note{
  Given some response data, the \pkg{VGAM} family function
  \code{\link{otpospoisson}} estimates the parameter \code{lambda}.



}

\seealso{
  \code{\link{otpospoisson}},
  \code{\link{Pospois}},
  \code{\link{Oipospois}}.


}
\examples{
dotpospois(1:20, 0.5)
rotpospois(20, 0.5)

\dontrun{ lambda <- 4; x <- 1:10
plot(x, dotpospois(x, lambda = lambda), type = "h", ylim = 0:1,
     sub = "lambda=4", las = 1, col = "blue", ylab = "Probability",
     main = "1-truncated positive-Poisson distribution: blue=PMF; orange=CDF")
lines(x + 0.1, potpospois(x, lambda = lambda), col = "orange", lty = 3, type = "h") }
}
\keyword{distribution}



