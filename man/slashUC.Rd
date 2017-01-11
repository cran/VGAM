\name{Slash}
\alias{Slash}
\alias{dslash}
\alias{pslash}
\alias{rslash}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Slash Distribution }
\description{
  Density function, distribution function, and
  random generation for the slash distribution.

}
\usage{
dslash(x, mu = 0, sigma = 1, log = FALSE, smallno = .Machine$double.eps*1000)
pslash(q, mu = 0, sigma = 1, very.negative = -10000,
       lower.tail = TRUE, log.p = FALSE)
rslash(n, mu = 0, sigma = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{n}{
  Same as \code{\link[stats]{runif}}.


% number of observations. Must be a single positive integer.


  }
  \item{mu, sigma}{the mean and standard deviation of
  the univariate normal distribution.


  }
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm of the density is returned.


  }
  \item{very.negative}{
  Numeric, of length 1.
  A large negative value.
  For \code{(q-mu)/sigma} values less than this,
  the value 0 is returned because
  \code{\link[stats]{integrate}} tends to fail.
  A warning is issued.
  Similarly, if \code{(q-mu)/sigma} is greater than
  \code{abs(very.negative)} then 1 is returned
  with a warning.


  }
  \item{smallno}{
  See \code{\link{slash}}.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }
}
\details{
  See \code{\link{slash}}, the \pkg{VGAM} family function
  for estimating the two parameters by maximum likelihood estimation,
  for the formula of the probability density function and other details.


  Function \code{\link{pslash}} uses a \code{for ()} loop and
  \code{\link[stats]{integrate}}, meaning it's very slow.
  It may also be inaccurate for extreme values of \code{q},
  and returns with 1 or 0 values when too extreme compared
  to \code{very.negative}.



}
\value{
  \code{dslash} gives the density, and
  \code{pslash} gives the distribution function,
  \code{rslash} generates random deviates.


}
%\references{ }
\author{ Thomas W. Yee and C. S. Chee}
\note{
  \code{pslash} is very slow.

}

\seealso{
  \code{\link{slash}}.


}
\examples{
\dontrun{
curve(dslash, col = "blue", ylab = "f(x)", -5, 5, ylim = c(0, 0.4), las = 1,
     main = "Standard slash, normal and Cauchy densities", lwd = 2)
curve(dnorm, col = "black", lty = 2, lwd = 2, add = TRUE)
curve(dcauchy, col = "orange", lty = 3, lwd = 2, add = TRUE)
legend("topleft", c("slash", "normal", "Cauchy"), lty = 1:3,
       col = c("blue","black","orange"), lwd = 2)

curve(pslash, col = "blue", -5, 5, ylim = 0:1)
pslash(c(-Inf, -20000, 20000, Inf))  # Gives a warning
}
}
\keyword{distribution}
