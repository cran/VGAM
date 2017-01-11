\name{powerlink}
\alias{powerlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Power Link Function }
\description{
  Computes the power transformation, including its inverse and the
  first two derivatives.

}
\usage{
powerlink(theta, power = 1, inverse = FALSE, deriv = 0,
          short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{power}{
  This denotes the power or exponent.


  }

  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The power link function raises a parameter by a certain value of
  \code{power}.
  Care is needed because it is very easy to get numerical
  problems, e.g., if \code{power=0.5} and \code{theta} is
  negative.



}
\value{
  For \code{powerlink} with \code{deriv = 0}, then \code{theta} raised
  to the power of \code{power}.
  And if \code{inverse = TRUE} then
  \code{theta} raised to the power of \code{1/power}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


}
%\references{
%  McCullagh, P. and Nelder, J. A. (1989)
%  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.
%
%}
\author{ Thomas W. Yee }

\note{
  Numerical problems may occur for certain combinations of
  \code{theta} and \code{power}.
  Consequently this link function should be used with caution.


}

\seealso{
  \code{\link{Links}},
  \code{\link{loge}}.


}
\examples{
powerlink("a", power = 2, short = FALSE, tag = TRUE)
powerlink(x <- 1:5)
powerlink(x, power = 2)
max(abs(powerlink(powerlink(x, power = 2),
                  power = 2, inverse = TRUE) - x))  # Should be 0
powerlink(x <- (-5):5, power = 0.5)  # Has NAs

# 1/2 = 0.5
pdata <- data.frame(y = rbeta(n = 1000, shape1 = 2^2, shape2 = 3^2))
fit <- vglm(y ~ 1, betaR(lshape1 = powerlink(power = 0.5), i1 = 3,
                         lshape2 = powerlink(power = 0.5), i2 = 7), data = pdata)
t(coef(fit, matrix = TRUE))
Coef(fit)  # Useful for intercept-only models
vcov(fit, untransform = TRUE)
}
\keyword{math}
\keyword{models}
\keyword{regression}

