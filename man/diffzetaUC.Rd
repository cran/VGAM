\name{Diffzeta}
\alias{Diffzeta}
\alias{ddiffzeta}
\alias{pdiffzeta}
\alias{qdiffzeta}
\alias{rdiffzeta}
\title{ Differenced Zeta Distribution }
\description{
  Density, distribution function,
  quantile function,
  and random generation
  for the differenced zeta distribution.



}
\usage{
ddiffzeta(x, shape, start = 1, log = FALSE)
pdiffzeta(q, shape, start = 1, lower.tail = TRUE)
qdiffzeta(p, shape, start = 1)
rdiffzeta(n, shape, start = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
  Same as in \code{\link[stats]{runif}}.


  }
  \item{shape, start}{
   Details at \code{\link{diffzeta}}.


%   For \code{rdiffzeta()} this pa%arameter must be of length 1.


  }
  \item{log, lower.tail}{
  Same as in \code{\link[stats]{runif}}.


  }
}
\details{
  This distribution appears to work well on the distribution
  of English words in such texts.
  Some more details are given in \code{\link{diffzeta}}.


}
\value{
  \code{ddiffzeta} gives the density,
  \code{pdiffzeta} gives the distribution function,
  \code{qdiffzeta} gives the quantile function, and
  \code{rdiffzeta} generates random deviates.




}
%\references{
%}
\author{ T. W. Yee }
\note{
  Given some response data, the \pkg{VGAM} family function
  \code{\link{diffzeta}} estimates the parameter \code{shape}.


  Function \code{pdiffzeta()} suffers from the problems that
  \code{\link{plog}} sometimes has, i.e., when \code{p}
  is very close to 1.



}

\seealso{
  \code{\link{diffzeta}},
  \code{\link{zetaff}},
  \code{\link{zipf}},
  \code{\link{Oizeta}}.


}
\examples{
ddiffzeta(1:20, 0.5, start = 2)
rdiffzeta(20, 0.5)

\dontrun{ shape <- 0.8; x <- 1:10
plot(x, ddiffzeta(x, shape = shape), type = "h", ylim = 0:1,
     sub = "shape=0.8", las = 1, col = "blue", ylab = "Probability",
     main = "Differenced zeta distribution: blue=PMF; orange=CDF")
lines(x + 0.1, pdiffzeta(x, shape = shape), col = "orange", lty = 3, type = "h") }
}
\keyword{distribution}



