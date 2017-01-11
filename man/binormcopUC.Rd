\name{Binormcop}
\alias{Binormcop}
\alias{dbinormcop}
\alias{pbinormcop}
\alias{rbinormcop}
\title{Gaussian Copula (Bivariate) Distribution}
\description{
  Density, distribution function,
  and random generation
  for the (one parameter) bivariate
  Gaussian copula distribution.


}
\usage{
dbinormcop(x1, x2, rho = 0, log = FALSE)
pbinormcop(q1, q2, rho = 0)
rbinormcop(n, rho = 0)
}
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.
  The \code{x1} and \code{x2} should be in the interval \eqn{(0,1)}.
  Ditto for \code{q1} and \code{q2}.


  }
  \item{n}{number of observations.
    Same as \code{\link[stats]{rnorm}}.

  }
  \item{rho}{the correlation parameter.
  Should be in the interval \eqn{(-1,1)}.


  }
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm is returned.
%   Same as \code{\link[stats]{rnorm}}.


  }
}
\value{
  \code{dbinormcop} gives the density,
  \code{pbinormcop} gives the distribution function, and
  \code{rbinormcop} generates random deviates (a two-column matrix).


}
%\references{
%
%}

\author{ T. W. Yee }
\details{
  See \code{\link{binormalcop}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation,
  for the formula of the
  cumulative distribution function and other details.


}
\note{
  Yettodo: allow \code{x1} and/or \code{x2} to have values 1,
  and to allow any values for \code{x1} and/or \code{x2} to be
  outside the unit square.


}
\seealso{
  \code{\link{binormalcop}},
  \code{\link{binormal}}.


}
\examples{
\dontrun{ edge <- 0.01  # A small positive value
N <- 101; x <- seq(edge, 1.0 - edge, len = N); Rho <- 0.7
ox <- expand.grid(x, x)
zedd <- dbinormcop(ox[, 1], ox[, 2], rho = Rho, log = TRUE)
contour(x, x, matrix(zedd, N, N), col = "blue", labcex = 1.5)
zedd <- pbinormcop(ox[, 1], ox[, 2], rho = Rho)
contour(x, x, matrix(zedd, N, N), col = "blue", labcex = 1.5)
}
}
\keyword{distribution}


%plot(r <- rbinormcop(n = 3000, rho = Rho), col = "blue")
%par(mfrow = c(1, 2))
%hist(r[, 1])  # Should be uniform
%hist(r[, 2])  # Should be uniform
