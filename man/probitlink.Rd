\name{probitlink}
\alias{probitlink}
%\alias{probit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Probit Link Function }
\description{
  Computes the probit transformation, including its inverse and the
  first two derivatives.

}
\usage{
probitlink(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{bvalue}{
  See \code{\link{Links}}.


  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The probit link function is commonly used for parameters that
  lie in the unit interval.
  It is the inverse CDF of the standard normal distribution.
  Numerical values of \code{theta} close to 0 or 1 or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.



}
\value{
  For \code{deriv = 0}, the probit of \code{theta}, i.e.,
  \code{qnorm(theta)} when \code{inverse = FALSE}, and if \code{inverse =
  TRUE} then \code{pnorm(theta)}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


}
\references{
  McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
\author{ Thomas W. Yee }

\note{
  Numerical instability may occur when \code{theta} is close to 1 or 0.
  One way of overcoming this is to use \code{bvalue}.


  In terms of the threshold approach with cumulative probabilities for
  an ordinal response this link function corresponds to the univariate
  normal distribution (see \code{\link{uninormal}}).


}
\seealso{
  \code{\link{Links}},
  \code{\link{logitlink}},
  \code{\link{clogloglink}},
  \code{\link{cauchitlink}},
  \code{\link[stats]{Normal}}.


}
\examples{
p <- seq(0.01, 0.99, by = 0.01)
probitlink(p)
max(abs(probitlink(probitlink(p), inverse = TRUE) - p))  # Should be 0

p <- c(seq(-0.02, 0.02, by = 0.01), seq(0.97, 1.02, by = 0.01))
probitlink(p)  # Has NAs
probitlink(p, bvalue = .Machine$double.eps)  # Has no NAs

\dontrun{p <- seq(0.01, 0.99, by = 0.01); par(lwd = (mylwd <- 2))
plot(p, logitlink(p), type = "l", col = "limegreen", ylab = "transformation",
     las = 1, main = "Some probability link functions")
lines(p,  probitlink(p), col = "purple")
lines(p, clogloglink(p), col = "chocolate")
lines(p, cauchitlink(p), col = "tan")
abline(v = 0.5, h = 0, lty = "dashed")
legend(0.1, 4, c("logitlink", "probitlink", "clogloglink", "cauchitlink"),
       col = c("limegreen", "purple", "chocolate", "tan"), lwd = mylwd)
par(lwd = 1) }
}
\keyword{math}
\keyword{models}
\keyword{regression}
