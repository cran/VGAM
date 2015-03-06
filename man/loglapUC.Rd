\name{loglapUC}
\alias{dloglap}
\alias{ploglap}
\alias{qloglap}
\alias{rloglap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Log-Laplace Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the 3-parameter log-Laplace distribution with location
  parameter \code{location.ald}, scale parameter \code{scale.ald}
  (on the log scale),
  and asymmetry parameter \code{kappa}.


}
\usage{
dloglap(x, location.ald = 0, scale.ald = 1,
        tau = 0.5, kappa = sqrt(tau/(1-tau)), log = FALSE)
ploglap(q, location.ald = 0, scale.ald = 1,
        tau = 0.5, kappa = sqrt(tau/(1-tau)), lower.tail = TRUE, log.p = FALSE)
qloglap(p, location.ald = 0, scale.ald = 1,
        tau = 0.5, kappa = sqrt(tau/(1-tau)), lower.tail = TRUE, log.p = FALSE)
rloglap(n, location.ald = 0, scale.ald = 1,
        tau = 0.5, kappa = sqrt(tau/(1-tau)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{
  number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.
  }
  \item{location.ald, scale.ald}{
    the location parameter \eqn{\xi}{xi} and
    the (positive) scale parameter \eqn{\sigma}{sigma},
    on the log scale.
  }
  \item{tau}{
  the quantile parameter \eqn{\tau}{tau}.
  Must consist of values in \eqn{(0,1)}.
  This argument is used to specify \code{kappa} and is ignored
  if \code{kappa} is assigned.

  }
  \item{kappa}{
  the asymmetry parameter \eqn{\kappa}{kappa}.
  Must consist of positive values.
  }
  \item{log}{
  if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }

}
\details{
  A positive random variable \eqn{Y} is said to have a log-Laplace
  distribution if \eqn{\log(Y)} has an asymmetric Laplace distribution
  (ALD). There are many variants of ALDs and the one used here
  is described in \code{\link{alaplace3}}.


}
\value{
  \code{dloglap} gives the density,
  \code{ploglap} gives the distribution function,
  \code{qloglap} gives the quantile function, and
  \code{rloglap} generates random deviates.


}
\references{

Kozubowski, T. J. and Podgorski, K. (2003)
Log-Laplace distributions.
\emph{International Mathematical Journal},
\bold{3}, 467--495.



}
\author{ T. W. Yee and Kai Huang }
%\note{
%  The \pkg{VGAM} family function \code{\link{loglaplace3}}
%  estimates the three parameters by maximum likelihood estimation.
%}

\seealso{
  \code{\link{dalap}},
  \code{\link{alaplace3}},
  \code{\link{loglaplace1}}.


% \code{\link{loglaplace3}}.



}
\examples{
loc <- 0; sigma <- exp(0.5); kappa <- 1
x <- seq(-0.2, 5, by = 0.01)
\dontrun{
plot(x, dloglap(x, loc, sigma, kappa = kappa), type = "l", col = "blue",
     main = "Blue is density, red is cumulative distribution function",
     ylim = c(0,1), sub = "Purple are 5,10,...,95 percentiles", las = 1, ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(qloglap(seq(0.05,0.95,by = 0.05), loc, sigma, kappa = kappa),
      dloglap(qloglap(seq(0.05,0.95,by = 0.05), loc, sigma, kappa = kappa),
                loc, sigma, kappa = kappa), col = "purple", lty = 3, type = "h")
lines(x, ploglap(x, loc, sigma, kappa = kappa), type = "l", col = "red")
abline(h = 0, lty = 2)
}
ploglap(qloglap(seq(0.05,0.95,by = 0.05), loc, sigma, kappa = kappa),
          loc, sigma, kappa = kappa)
}
\keyword{distribution}

