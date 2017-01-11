\name{cauchit}
\alias{cauchit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cauchit Link Function }
\description{
  Computes the cauchit (tangent) link transformation, including its
  inverse and the first two derivatives.

}
\usage{
cauchit(theta, bvalue = .Machine$double.eps,
        inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
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
  This link function is an alternative link function for parameters that
  lie in the unit interval.
  This type of link bears the same relation to the Cauchy distribution
  as the probit link bears to the Gaussian. One characteristic of this
  link function is that the tail is heavier relative to the other links
  (see examples below).


  Numerical values of \code{theta} close to 0 or 1 or out of range result
  in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.


}
\value{
  For \code{deriv = 0}, the tangent of \code{theta}, i.e.,
  \code{tan(pi * (theta-0.5))} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{0.5 + atan(theta)/pi}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of
  \code{theta} if \code{inverse = FALSE}, else if \code{inverse = TRUE}
  then it returns the reciprocal.


}
\references{
  McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
\author{ Thomas W. Yee }

\note{
  Numerical instability may occur when \code{theta} is close to 1 or 0.
  One way of overcoming this is to use \code{bvalue}.


  As mentioned above,
  in terms of the threshold approach with cumulative probabilities for
  an ordinal response this link function corresponds to the
  Cauchy distribution (see \code{\link{cauchy1}}).


}

\seealso{
    \code{\link{logit}},
    \code{\link{probit}},
    \code{\link{cloglog}},
    \code{\link{loge}},
    \code{\link{cauchy}},
    \code{\link{cauchy1}}.


 }
\examples{
p <- seq(0.01, 0.99, by=0.01)
cauchit(p)
max(abs(cauchit(cauchit(p), inverse = TRUE) - p))  # Should be 0

p <- c(seq(-0.02, 0.02, by=0.01), seq(0.97, 1.02, by = 0.01))
cauchit(p)  # Has no NAs

\dontrun{
par(mfrow = c(2, 2), lwd = (mylwd <- 2))
y <- seq(-4, 4, length = 100)
p <- seq(0.01, 0.99, by = 0.01)

for (d in 0:1) {
  matplot(p, cbind(logit(p, deriv = d), probit(p, deriv = d)),
          type = "n", col = "purple", ylab = "transformation",
          las = 1, main = if (d == 0) "Some probability link functions"
          else "First derivative")
  lines(p,   logit(p, deriv = d), col = "limegreen")
  lines(p,  probit(p, deriv = d), col = "purple")
  lines(p, cloglog(p, deriv = d), col = "chocolate")
  lines(p, cauchit(p, deriv = d), col = "tan")
  if (d == 0) {
    abline(v = 0.5, h = 0, lty = "dashed")
    legend(0, 4.5, c("logit", "probit", "cloglog", "cauchit"), lwd = mylwd,
           col = c("limegreen","purple","chocolate", "tan"))
  } else
    abline(v = 0.5, lty = "dashed")
}

for (d in 0) {
  matplot(y, cbind( logit(y, deriv = d, inverse = TRUE),
                   probit(y, deriv = d, inverse = TRUE)),
          type  = "n", col = "purple", xlab = "transformation", ylab = "p",
          main = if (d == 0) "Some inverse probability link functions"
          else "First derivative", las=1)
  lines(y,   logit(y, deriv = d, inverse = TRUE), col = "limegreen")
  lines(y,  probit(y, deriv = d, inverse = TRUE), col = "purple")
  lines(y, cloglog(y, deriv = d, inverse = TRUE), col = "chocolate")
  lines(y, cauchit(y, deriv = d, inverse = TRUE), col = "tan")
  if (d == 0) {
      abline(h = 0.5, v = 0, lty = "dashed")
      legend(-4, 1, c("logit", "probit", "cloglog", "cauchit"), lwd = mylwd,
             col = c("limegreen", "purple", "chocolate", "tan"))
  }
}
par(lwd = 1)
}
}
\keyword{math}
\keyword{models}
\keyword{regression}

%plot(y, logit(y, inverse = TRUE), type = "l", col = "limegreen",
%     xlab = "transformation", ylab = "p",
%     lwd=2, las=1, main = "Some inverse probability link functions")
%lines(y, probit(y, inverse = TRUE), col = "purple", lwd=2)
%lines(y, cloglog(y, inverse = TRUE), col = "chocolate", lwd=2)
%abline(h=0.5, v = 0, lty = "dashed")




