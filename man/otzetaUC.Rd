\name{Otzeta}
\alias{Otzeta}
\alias{dotzeta}
\alias{potzeta}
\alias{qotzeta}
\alias{rotzeta}
\title{ One-truncated Zeta Distribution }
\description{
  Density, distribution function,
  quantile function,
  and random generation
  for the one-truncated zeta distribution.



}
\usage{
dotzeta(x, shape, log = FALSE)
potzeta(q, shape, log.p = FALSE)
qotzeta(p, shape)
rotzeta(n, shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
  Same as in \code{\link[stats]{runif}}.


  }
  \item{shape}{
   The positive shape parameter described in in \code{\link{zetaff}}.
   Here it is called \code{shape} because it is positive.


%   For \code{rotzeta()} this pa%arameter must be of length 1.


  }
  \item{log, log.p}{
  Same as in \code{\link[stats]{runif}}.


  }
}
\details{
  The one-truncated zeta distribution is a zeta distribution but with
  the probability of a one being zero. The other probabilities are scaled
  to add to unity.
  Some more details are given in \code{\link{zetaff}}.


}
\value{
  \code{dotzeta} gives the density,
  \code{potzeta} gives the distribution function,
  \code{qotzeta} gives the quantile function, and
  \code{rotzeta} generates random deviates.




}
%\references{
%}
\author{ T. W. Yee }
\note{
  Given some response data, the \pkg{VGAM} family function
  \code{\link{otzeta}} estimates the parameter \code{shape}.


% Function \code{potzeta()} suffers from the problems that
% \code{\link{plog}} sometimes has.



}

\seealso{
  \code{\link{Otzeta}},
  \code{\link{zetaff}},
  \code{\link{Oizeta}}.


}
\examples{
dotzeta(1:20, 0.5)
rotzeta(20, 0.5)

\dontrun{ shape <- 0.8; x <- 1:10
plot(x, dotzeta(x, shape = shape), type = "h", ylim = 0:1,
     sub = "shape=0.8", las = 1, col = "blue", ylab = "Probability",
     main = "1-truncated zeta distribution: blue=PMF; orange=CDF")
lines(x + 0.1, potzeta(x, shape = shape), col = "orange", lty = 3, type = "h") }
}
\keyword{distribution}



