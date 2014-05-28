\name{PoissonPoints}
\alias{PoissonPoints}
\alias{dpois.points}
%\alias{ppois.points}
%\alias{qpois.points}
\alias{rpois.points}
\title{Poisson Points Distribution}
\description{
  Density
  for the
  PoissonPoints distribution.


% distribution function, quantile function
% and random generation

}
\usage{
dpois.points(x, lambda, ostatistic, dimension = 2, log = FALSE)
}
%ppois.points(q, lambda, ostatistic, dimension = 2, log = FALSE)
%qpois.points(p, lambda, ostatistic, dimension = 2, log = FALSE)
%rpois.points(n, lambda, ostatistic, dimension = 2, log = FALSE)
\arguments{
  \item{x}{vector of quantiles.}
  \item{lambda}{
  the mean density of points.


  }
  \item{ostatistic}{
  positive values, usually integers.


  }
  \item{dimension}{
  Either 2 and/or 3.


  }
% \item{p}{vector of probabilities.}
% \item{n}{number of observations.
% Same as \code{\link[stats:Uniform]{runif}}.


% }
  \item{log}{
  Logical; if TRUE, the logarithm is returned.


  }
}
\value{
  \code{dpois.points} gives the density. % and


% \code{ppois.points} gives the distribution function,
% \code{qpois.points} gives the quantile function, and
% \code{rpois.points} generates random deviates.


}
%\author{ T. W. Yee }
\details{
  See \code{\link{poisson.points}}, the \pkg{VGAM} family function
  for estimating the parameters,
  for the formula of the probability density function and other details.


}
%\section{Warning }{


%}
\seealso{
  \code{\link{poisson.points}},
  \code{\link[stats:Poisson]{dpois}},
  \code{\link{Maxwell}}.


}
\examples{
\dontrun{ lambda <- 1; xvec <- seq(0, 2, length = 400)
plot(xvec, dpois.points(xvec, lambda, ostat = 1, dimension = 2),
     type = "l", las = 1, col = "blue",
     sub = "First order statistic",
     main = paste("PDF of PoissonPoints distribution with lambda = ",
                  lambda, " and on the plane", sep = "")) }
}
\keyword{distribution}


