\name{Hzeta}
\alias{Hzeta}
\alias{dhzeta}
\alias{phzeta}
\alias{qhzeta}
\alias{rhzeta}
\title{ Haight's Zeta Distribution  }
\description{
  Density, distribution function, quantile function and random generation
  for Haight's zeta distribution with parameter \code{shape}.


}
\usage{
dhzeta(x, shape, log = FALSE)
phzeta(q, shape, log.p = FALSE)
qhzeta(p, shape)
rhzeta(n, shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
  Same meaning as \code{\link[stats]{runif}}.


  }
  \item{shape}{
   The positive shape parameter.
   Called \eqn{\alpha}{alpha} below.


  }
  \item{log,log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }
}
\details{
   The probability function is
   \deqn{f(x) = (2x-1)^{(-\alpha)} - (2x+1)^{(-\alpha)},}{%
         f(x) = (2x-1)^(-alpha) - (2x+1)^(-alpha),}
   where \eqn{\alpha>0}{alpha>0} and \eqn{x=1,2,\ldots}{x=1,2,...}.


}
\value{
  \code{dhzeta} gives the density,
  \code{phzeta} gives the distribution function,
  \code{qhzeta} gives the quantile function, and
  \code{rhzeta} generates random deviates.


}
%\references{
%
%    Pages 533--4 of
%    Johnson N. L., Kemp, A. W. and Kotz S. (2005)
%    \emph{Univariate Discrete Distributions},
%    3rd edition,
%    Hoboken, New Jersey: Wiley.
%
%
%}
\author{ T. W. Yee and Kai Huang }
\note{
  Given some response data, the \pkg{VGAM} family function
  \code{\link{hzeta}} estimates the parameter \code{shape}.


}

\seealso{
  \code{\link{hzeta}},
  \code{\link{zeta}},
  \code{\link{zetaff}},
  \code{\link{simulate.vlm}}.


}
\examples{
dhzeta(1:20, 2.1)
rhzeta(20, 2.1)

round(1000 * dhzeta(1:8, 2))
table(rhzeta(1000, 2))

\dontrun{ shape <- 1.1; x <- 1:10
plot(x, dhzeta(x, shape = shape), type = "h", ylim = 0:1, lwd = 2,
     sub = paste("shape =", shape), las = 1, col = "blue", ylab = "Probability",
     main = "Haight's zeta: blue = density; orange = distribution function")
lines(x+0.1, phzeta(x, shape = shape), col = "orange", lty = 3, lwd = 2,
      type = "h")
}
}
\keyword{distribution}
