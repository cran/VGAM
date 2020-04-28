\name{Otlog}
\alias{Otlog}
\alias{dotlog}
\alias{potlog}
\alias{qotlog}
\alias{rotlog}
\title{ One-truncated Logarithmic Distribution }
\description{
  Density, distribution function,
  quantile function,
  and random generation
  for the one-truncated logarithmic distribution.



}
\usage{
dotlog(x, shape, log = FALSE)
potlog(q, shape, log.p = FALSE)
qotlog(p, shape)
rotlog(n, shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{
   Vector of quantiles. For the density, it should be a vector with
   integer values \eqn{> 1} in order for the probabilities to be positive.


  }
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as in \code{\link[stats]{runif}}.


  }
  \item{shape}{
   The parameter value \eqn{c} described in in \code{\link{logff}}.
   Here it is called \code{shape} because \eqn{0<c<1} is the range.


%   For \code{rotlog()} this pa%arameter must be of length 1.


  }
  \item{log, log.p}{
  Logical.
  If \code{log.p = TRUE} then all probabilities \code{p} are
  given as \code{log(p)}.


  }
}
\details{
  The one-truncated logarithmic distribution is a logarithmic distribution but with
  the probability of a one being zero. The other probabilities are scaled
  to add to unity.
  Some more details are given in \code{\link{logff}}.


}
\value{
  \code{dotlog} gives the density,
  \code{potlog} gives the distribution function,
  \code{qotlog} gives the quantile function, and
  \code{rotlog} generates random deviates.




}
%\references{
%}
\author{ T. W. Yee }
\note{
  Given some response data, the \pkg{VGAM} family function
  \code{\link{otlog}} estimates the parameter \code{shape}.
  Function \code{potlog()} suffers from the problems that
  \code{\link{plog}} sometimes has.



}

\seealso{
  \code{\link{Gaitlog}},
  \code{\link{otlog}},
  \code{\link{rlog}},
  \code{\link{Oilog}}.


}
\examples{
dotlog(1:20, 0.5)
rotlog(20, 0.5)

\dontrun{ shape <- 0.8; x <- 1:10
plot(x, dotlog(x, shape = shape), type = "h", ylim = 0:1,
     sub = "shape=0.8", las = 1, col = "blue", ylab = "Probability",
     main = "1-truncated logarithmic distribution: blue=PMF; orange=CDF")
lines(x + 0.1, potlog(x, shape = shape), col = "orange", lty = 3, type = "h") }
}
\keyword{distribution}



