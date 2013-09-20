\name{Hzeta}
\alias{Hzeta}
\alias{dhzeta}
\alias{phzeta}
\alias{qhzeta}
\alias{rhzeta}
\title{ Haight's Zeta Distribution  }
\description{
  Density, distribution function, quantile function and random generation
  for Haight's zeta distribution with parameter \code{alpha}.


}
\usage{
dhzeta(x, alpha, log = FALSE)
phzeta(q, alpha)
qhzeta(p, alpha)
rhzeta(n, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{
   Vector of quantiles. For the density, it should be a vector with
   positive integer values in order for the probabilities to be positive.


  }
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as \code{\link[stats]{runif}}.


  }
  \item{alpha}{ 
   The parameter value. Must contain positive values and is recycled to
   the length of \code{x} or \code{p} or \code{q} if necessary.


  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


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
\references{ 

    Page 533 of
    Johnson N. L., Kemp, A. W. and Kotz S. (2005)
    \emph{Univariate Discrete Distributions},
    3rd edition,
    Hoboken, New Jersey: Wiley.


}
\author{ T. W. Yee }
\note{ 
  Given some response data, the \pkg{VGAM} family function
  \code{\link{hzeta}} estimates the parameter \code{alpha}.


}

\seealso{
  \code{\link{hzeta}},
  \code{\link{zeta}},
  \code{\link{zetaff}}.


}
\examples{
dhzeta(1:20, 2.1)
rhzeta(20, 2.1)

round(1000 * dhzeta(1:8, 2))
table(rhzeta(1000, 2))

\dontrun{ alpha <- 1.1; x <- 1:10
plot(x, dhzeta(x, alpha = alpha), type = "h", ylim = 0:1, lwd = 2,
     sub = paste("alpha =", alpha), las = 1, col = "blue", ylab = "Probability",
     main = "Haight's zeta: blue = density; orange = distribution function")
lines(x+0.1, phzeta(x, alpha = alpha), col = "orange", lty = 3, lwd = 2, type = "h")
}
}
\keyword{distribution}
