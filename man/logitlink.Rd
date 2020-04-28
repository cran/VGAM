\name{logitlink}
\alias{logitlink}
%\alias{logit}
\alias{extlogitlink}
%\alias{extlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Logit Link Function }
\description{
  Computes the logit transformation, including its inverse and the
  first two derivatives.

}
\usage{
logitlink(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
          short = TRUE, tag = FALSE)
extlogitlink(theta, min = 0, max = 1, bminvalue = NULL, bmaxvalue = NULL,
             inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.


  }
  \item{bvalue, bminvalue, bmaxvalue}{
  See \code{\link{Links}}.
  These are boundary values.
  For \code{extlogitlink}, values of \code{theta} less than or equal
  to \eqn{A} or greater than or equal to \eqn{B} can be replaced
  by \code{bminvalue} and \code{bmaxvalue}.


  }

% Extra argument for passing in additional information.
% For \code{logitlink}, values of \code{theta} which are equal to 0 or 1 are
% replaced by \code{earg} or \code{1-earg}
% (respectively, and if given) before computing the logit.

  \item{min, max}{
  For \code{extlogitlink},
  \code{min} gives \eqn{A},
  \code{max} gives \eqn{B}, and for out of range values,
  \code{bminvalue} and \code{bmaxvalue}.

  }

  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The logit link function is very commonly used for parameters that
  lie in the unit interval.
  It is the inverse CDF of the logistic distribution.
  Numerical values of \code{theta} close to 0 or 1 or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.

  

  The \emph{extended} logit link function \code{extlogitlink} should be used
  more generally for parameters that lie in the interval \eqn{(A,B)},
  say.
  The formula is
  \deqn{\log((\theta-A)/(B-\theta))}{%
         log((theta-A)/(B-theta))}
  and the default values for \eqn{A} and \eqn{B} correspond to
  the ordinary
  logit function.
  Numerical values of \code{theta} close to \eqn{A} or \eqn{B} or
  out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.
  However these can be replaced by values \eqn{bminvalue} and
  \eqn{bmaxvalue} first before computing the link function.



}
\value{
  For \code{logitlink} with \code{deriv = 0}, the logit of \code{theta}, i.e.,
  \code{log(theta/(1-theta))} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{exp(theta)/(1+exp(theta))}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


  Here, all logarithms are natural logarithms, i.e., to base \emph{e}.


}
\references{
  McCullagh, P. and Nelder, J. A. (1989).
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
\author{ Thomas W. Yee }

\note{
  Numerical instability may occur when \code{theta} is
  close to 1 or 0 (for \code{logitlink}), or close to \eqn{A} or \eqn{B} for
  \code{extlogitlink}.
  One way of overcoming this is to use, e.g., \code{bvalue}.


  In terms of the threshold approach with cumulative probabilities for
  an ordinal response this link function corresponds to the univariate
  logistic distribution (see \code{\link{logistic}}).


}

\seealso{
    \code{\link{Links}},
    \code{\link{logitoffsetlink}},
    \code{\link{probitlink}},
    \code{\link{clogloglink}},
    \code{\link{cauchitlink}},
    \code{\link{logistic1}},
    \code{\link{loglink}},
    \code{\link[stats]{Logistic}},
    \code{\link{multilogitlink}}.


}
\examples{
p <- seq(0.01, 0.99, by = 0.01)
logitlink(p)
max(abs(logitlink(logitlink(p), inverse = TRUE) - p))  # Should be 0

p <- c(seq(-0.02, 0.02, by = 0.01), seq(0.97, 1.02, by = 0.01))
logitlink(p)  # Has NAs
logitlink(p, bvalue = .Machine$double.eps)  # Has no NAs

p <- seq(0.9, 2.2, by = 0.1)
extlogitlink(p, min = 1, max = 2,
             bminvalue = 1 + .Machine$double.eps,
             bmaxvalue = 2 - .Machine$double.eps)  # Has no NAs

\dontrun{ par(mfrow = c(2,2), lwd = (mylwd <- 2))
y <- seq(-4, 4, length = 100)
p <- seq(0.01, 0.99, by = 0.01)
for (d in 0:1) {
  myinv <- (d > 0)
  matplot(p, cbind( logitlink(p, deriv = d, inverse = myinv),
                   probitlink(p, deriv = d, inverse = myinv)),
          type = "n", col = "purple", ylab = "transformation", las = 1,
          main = if (d ==  0) "Some probability link functions"
          else "1 / first derivative")
  lines(p,   logitlink(p, deriv = d, inverse = myinv), col = "limegreen")
  lines(p,  probitlink(p, deriv = d, inverse = myinv), col = "purple")
  lines(p, clogloglink(p, deriv = d, inverse = myinv), col = "chocolate")
  lines(p, cauchitlink(p, deriv = d, inverse = myinv), col = "tan")
  if (d ==  0) {
    abline(v = 0.5, h = 0, lty = "dashed")
    legend(0, 4.5, c("logitlink", "probitlink", "clogloglink", "cauchitlink"),
           col = c("limegreen", "purple", "chocolate", "tan"), lwd = mylwd)
  } else
    abline(v = 0.5, lty = "dashed")
}

for (d in 0) {
  matplot(y, cbind(logitlink(y, deriv = d, inverse = TRUE),
                   probitlink(y, deriv = d, inverse = TRUE)), las = 1,
          type = "n", col = "purple", xlab = "transformation", ylab = "p",
          main = if (d ==  0) "Some inverse probability link functions"
          else "First derivative")
  lines(y,   logitlink(y, deriv = d, inverse = TRUE), col = "limegreen")
  lines(y,  probitlink(y, deriv = d, inverse = TRUE), col = "purple")
  lines(y, clogloglink(y, deriv = d, inverse = TRUE), col = "chocolate")
  lines(y, cauchitlink(y, deriv = d, inverse = TRUE), col = "tan")
  if (d ==  0) {
    abline(h = 0.5, v = 0, lty = "dashed")
    legend(-4, 1, c("logitlink", "probitlink", "clogloglink", "cauchitlink"),
           col = c("limegreen", "purple", "chocolate", "tan"), lwd = mylwd)
  }
}

p <- seq(0.21, 0.59, by = 0.01)
plot(p, extlogitlink(p, min = 0.2, max = 0.6),
     type = "l", col = "black", ylab = "transformation", xlim = c(0, 1),
     las = 1, main = "extlogitlink(p, min = 0.2, max = 0.6)")
par(lwd = 1)
}
}
\keyword{math}
\keyword{models}
\keyword{regression}

%plot(y, logitlink(y, inverse = TRUE), type = "l", col = "limegreen",
%     xlab = "transformation", ylab = "p",
%     lwd = 2, las = 1, main = "Some inverse probability link functions")
%lines(y, probitlink(y, inverse = TRUE), col = "purple", lwd = 2)
%lines(y, clogloglink(y, inverse = TRUE), col = "chocolate", lwd = 2)
%abline(h = 0.5, v = 0, lty = "dashed")




