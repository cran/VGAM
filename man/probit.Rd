\name{probit}
\alias{probit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Probit Link Function }
\description{
  Computes the probit transformation, including its inverse and the
  first two derivatives.
}
\usage{
probit(theta, earg = list(), inverse = FALSE, deriv = 0,
       short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{earg}{
  Optional list. Extra argument for passing in additional information.
  Values of \code{theta} which are less than or equal to 0 can be
  replaced by the \code{bvalue} component of the list \code{earg}
  before computing the link function value.
  Values of \code{theta} which are greater than or equal to 1 can be
  replaced by 1 minus the \code{bvalue} component of the list \code{earg}
  before computing the link function value.
  The component name \code{bvalue} stands for ``boundary value''.
  See \code{\link{Links}} for general information about \code{earg}.

  }
  \item{inverse}{
   Logical. If \code{TRUE} the inverse function is computed.

  }
  \item{deriv}{
   Order of the derivative. Integer with value 0, 1 or 2.

  }
  \item{short}{
  Used for labelling the \code{blurb} slot of a \code{\link{vglmff-class}}
  object.

  }
  \item{tag}{
  Used for labelling the linear/additive predictor in the
  \code{initialize} slot of a \code{\link{vglmff-class}} object.
  Contains a little more information if \code{TRUE}.

  }
}
\details{
  The probit link function is commonly used for parameters that
  lie in the unit interval. 
  Numerical values of \code{theta} close to 0 or 1 or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.
  The arguments \code{short} and \code{tag} are used only if
  \code{theta} is character.

}
\value{
  For \code{deriv = 0}, the probit of \code{theta}, i.e.,
  \code{qnorm(theta)} when \code{inverse = FALSE}, and if \code{inverse =
  TRUE} then \code{pnorm(theta)}.

  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function of \code{theta}
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
  One way of overcoming this is to use \code{earg}.

  In terms of the threshold approach with cumulative probabilities for
  an ordinal response this link function corresponds to the univariate
  normal distribution (see \code{\link{normal1}}).

}
\seealso{ 
  \code{\link{Links}},
  \code{\link{logit}},
  \code{\link{cloglog}},
  \code{\link{cauchit}}.
}
\examples{
p = seq(0.01, 0.99, by=0.01)
probit(p)
max(abs(probit(probit(p), inverse=TRUE) - p)) # Should be 0

p = c(seq(-0.02, 0.02, by=0.01), seq(0.97, 1.02, by=0.01))
probit(p)  # Has NAs
probit(p, earg=list(bvalue= .Machine$double.eps))  # Has no NAs

\dontrun{
plot(p, logit(p), type="l", col="limegreen", ylab="transformation",
     lwd=2, las=1, main="Some probability link functions")
lines(p, probit(p), col="purple", lwd=2)
lines(p, cloglog(p), col="chocolate", lwd=2)
lines(p, cauchit(p), col="tan", lwd=2)
abline(v=0.5, h=0, lty="dashed")
legend(0.1, 4.0, c("logit", "probit", "cloglog", "cauchit"), 
       col=c("limegreen","purple","chocolate","tan"), lwd=2) }
}
\keyword{math}
\keyword{models}
\keyword{regression}
