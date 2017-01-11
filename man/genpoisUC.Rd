\name{dgenpois}
\alias{dgenpois}
%\alias{pgenpois}
%\alias{qgenpois}
%\alias{rgenpois}
\title{The Generalized Poisson Distribution}
\description{
  Density for the Generalized Poisson Distribution.

}
\usage{
dgenpois(x, lambda = 0, theta, log = FALSE)
}
\arguments{
  \item{x,}{vector of quantiles.}

  \item{lambda, theta}{
    See \code{\link{genpoisson}}.
    The default value of \code{lambda} corresponds to an
    ordinary Poisson distribution.


  }
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm of the density is returned.


  }
}
\value{
  \code{dgenpois} gives the density.
  The value \code{NaN} is returned for elements not satisfying
  the parameter restrictions, e.g., if \eqn{\lambda > 1}{lambda > 1}.



% \code{pgenpois} gives the distribution function, and
% \code{qgenpois} gives the quantile function, and
% \code{rgenpois} generates random deviates.


}
\author{ T. W. Yee }
\details{
  Most of the background to this function is given
  in \code{\link{genpoisson}}.
  Some warnings relevant to this distribution are given there,
  especially relating to the complicated range of the
  parameter \code{lambda} about or near \eqn{-1}.


  Note that numerical round off errors etc. can occur; see
  below for an example.


}
%\note{
%}
\seealso{
  \code{\link{genpoisson}},
  \code{\link[stats:Poisson]{dpois}}.


}
\examples{
sum(dgenpois(0:1000, lambda = -0.5, theta = 2))  # Not perfect...
\dontrun{
lambda <- -0.2; theta <- 2; y <- 0:10
proby <- dgenpois(y, lambda = lambda, theta = theta, log = FALSE)
plot(y, proby, type = "h", col = "blue", lwd = 2, ylab = "P[Y=y]",
     main = paste("Y ~ Generalized Poisson(lambda=", lambda,
                  ", theta=", theta, ")", sep = ""), las = 1,
     sub = "Orange is the Poisson probability function")
sum(proby)
lines(y + 0.1, dpois(y, theta), type = "h", lwd = 2, col = "orange")
}
}
\keyword{distribution}


