\name{alaplaceUC}
\alias{dalap}
\alias{palap}
\alias{qalap}
\alias{ralap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Laplace Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the 3-parameter asymmetric Laplace distribution with location
  parameter \code{location}, scale parameter \code{scale},
  and asymmetry parameter \code{kappa}.


}
\usage{
dalap(x, location = 0, scale = 1, tau = 0.5, kappa = sqrt(tau/(1-tau)),
      log = FALSE)
palap(q, location = 0, scale = 1, tau = 0.5, kappa = sqrt(tau/(1-tau)),
      lower.tail = TRUE, log.p = FALSE)
qalap(p, location = 0, scale = 1, tau = 0.5, kappa = sqrt(tau/(1-tau)),
      lower.tail = TRUE, log.p = FALSE)
ralap(n, location = 0, scale = 1, tau = 0.5, kappa = sqrt(tau/(1-tau)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{
  number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.


  }
  \item{location}{
    the location parameter \eqn{\xi}{xi}.


  }
  \item{scale}{
  the scale parameter \eqn{\sigma}{sigma}.
  Must consist of positive values.


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
  There are many variants of asymmetric Laplace distributions (ALDs) and
  this one is known as \emph{the} ALD by Kotz et al. (2001).
  See \code{\link{alaplace3}}, the \pkg{VGAM} family function
  for estimating the three parameters by maximum likelihood estimation,
  for formulae and details.


}
\value{
  \code{dalap} gives the density,
  \code{palap} gives the distribution function,
  \code{qalap} gives the quantile function, and
  \code{ralap} generates random deviates.


}
\references{

Kotz, S., Kozubowski, T. J. and Podgorski, K. (2001)
\emph{The Laplace distribution and generalizations:
a revisit with applications to communications,
economics, engineering, and finance},
Boston: Birkhauser.


}
\author{ T. W. Yee and Kai Huang }
%\note{
%  The \pkg{VGAM} family function \code{\link{alaplace3}}
%  estimates the three parameters by maximum likelihood estimation.
%}

\seealso{
  \code{\link{alaplace3}}.
% \code{\link{dloglap}}.



}
\examples{
x <- seq(-5, 5, by = 0.01)
loc <- 0; sigma <- 1.5; kappa <- 2
\dontrun{ plot(x, dalap(x, loc, sigma, kappa = kappa), type = "l", col = "blue",
     main = "Blue is density, orange is cumulative distribution function",
     ylim = c(0, 1), sub = "Purple are 5, 10, ..., 95 percentiles",
     las = 1, ylab = "", cex.main = 0.5)
abline(h = 0, col = "blue", lty = 2)
lines(qalap(seq(0.05, 0.95, by = 0.05), loc, sigma, kappa = kappa),
      dalap(qalap(seq(0.05, 0.95, by = 0.05), loc, sigma, kappa = kappa),
            loc, sigma, kappa = kappa), col = "purple", lty = 3, type = "h")
lines(x, palap(x, loc, sigma, kappa = kappa), type = "l", col = "orange")
abline(h = 0, lty = 2) }

pp <- seq(0.05, 0.95, by = 0.05)  # Test two functions
max(abs(palap(qalap(pp, loc, sigma, kappa = kappa),
              loc, sigma, kappa = kappa) - pp))  # Should be 0
}
\keyword{distribution}

