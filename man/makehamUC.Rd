\name{Makeham}
\alias{Makeham}
\alias{dmakeham}
\alias{pmakeham}
\alias{qmakeham}
\alias{rmakeham}
\title{The Makeham Distribution}
\description{
  Density,
  cumulative distribution function,
  quantile function
  and
  random generation for
  the Makeham distribution.

}
\usage{
dmakeham(x, shape, scale = 1, epsilon = 0, log = FALSE)
pmakeham(q, shape, scale = 1, epsilon = 0)
qmakeham(p, shape, scale = 1, epsilon = 0)
rmakeham(n, shape, scale = 1, epsilon = 0)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }
  \item{shape, scale}{positive shape and scale parameters. }
  \item{epsilon}{another parameter. Must be non-negative. See below. }

}
\value{
  \code{dmakeham} gives the density,
  \code{pmakeham} gives the cumulative distribution function,
  \code{qmakeham} gives the quantile function, and
  \code{rmakeham} generates random deviates.


}
\author{ T. W. Yee }
\details{
  See \code{\link{makeham}} for details.
  The default value of \code{epsilon = 0} corresponds
  to the Gompertz distribution.
  The function \code{\link{pmakeham}} uses \code{\link{lambertW}}.


}
\references{

Jodra, P. (2009) 
A closed-form expression for the quantile function of the
Gompertz-Makeham distribution.
\emph{Mathematics and Computers in Simulation},
\bold{79}, 3069--3075.



}


%\note{
%
%}
\seealso{
  \code{\link{makeham}},
  \code{\link{lambertW}}.


}
\examples{
probs <- seq(0.01, 0.99, by = 0.01)
Shape <- exp(-1); Scale <- exp(1); eps = Epsilon <- exp(-1)
max(abs(pmakeham(qmakeham(p = probs, Shape, sca = Scale, eps = Epsilon),
                 Shape, sca = Scale, eps = Epsilon) - probs)) # Should be 0

\dontrun{ x <- seq(-0.1, 2.0, by = 0.01);
plot(x, dmakeham(x, Shape, sca = Scale, eps = Epsilon), type = "l",
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles",
     col = "blue", las = 1, ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(x, pmakeham(x, Shape, sca = Scale, eps = Epsilon), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qmakeham(probs, Shape, sca = Scale, eps = Epsilon)
lines(Q, dmakeham(Q, Shape, sca = Scale, eps = Epsilon),
      col = "purple", lty = 3, type = "h")
pmakeham(Q, Shape, sca = Scale, eps = Epsilon) - probs # Should be all zero
abline(h = probs, col = "purple", lty = 3) }
}
\keyword{distribution}


