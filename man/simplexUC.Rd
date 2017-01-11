\name{Simplex }
\alias{dsimplex}
%\alias{psimplex}
%\alias{qsimplex}
\alias{rsimplex}
\title{ Simplex Distribution }
\description{
  Density function,
  and random generation for the simplex distribution.

}
\usage{
dsimplex(x, mu = 0.5, dispersion = 1, log = FALSE)
rsimplex(n, mu = 0.5, dispersion = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Vector of quantiles.
  The support of the distribution is the interval \eqn{(0,1)}.


  }
  \item{mu, dispersion}{
  Mean and dispersion parameters.
  The former lies in the interval \eqn{(0,1)} and the latter is positive.


  }
  \item{n, log}{
  Same usage as \code{\link[stats:Uniform]{runif}}.


  }
}
\details{

  The \pkg{VGAM} family function \code{\link{simplex}} fits this model;
  see that online help for more information.
  For \code{rsimplex()} the rejection method is used;
  it may be very slow if the density is highly peaked,
  and will fail if the density asymptotes at the boundary.


}
\value{
  \code{dsimplex(x)} gives the density function,
  \code{rsimplex(n)} gives \eqn{n} random variates.


}
% \references{
%
%}
\author{ T. W. Yee }

\seealso{
  \code{\link{simplex}}.


}

\examples{
sigma <- c(4, 2, 1)  # Dispersion parameter
mymu  <- c(0.1, 0.5, 0.7); xxx <- seq(0, 1, len = 501)
\dontrun{ par(mfrow = c(3, 3))  # Figure 2.1 of Song (2007)
for (iii in 1:3)
  for (jjj in 1:3) {
    plot(xxx, dsimplex(xxx, mymu[jjj], sigma[iii]),
         type = "l", col = "blue", xlab = "", ylab = "", main =
         paste("mu = ", mymu[jjj], ", sigma = ", sigma[iii], sep = "")) } }
}
\keyword{distribution}

% mean(rsimplex(1000, mymu[2], sigma[2]))  # Should be mu below
%  var(rsimplex(1000, mymu[2], sigma[2]))  # Should be as below
% (mu <- mymu[2])
% lambda <- (1 / sigma[2])^2
% mu * (1 - mu) - sqrt(lambda / 2) * exp(lambda / (mu^2 * (1 - mu)^2)) *
% pgamma(lambda / (2 * mu^2 * (1 - mu)^2), 0.5, lower = FALSE) * gamma(0.5)


