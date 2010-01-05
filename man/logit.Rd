\name{logit}
\alias{logit}
\alias{elogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Logit Link Function }
\description{
  Computes the logit transformation, including its inverse and the
  first two derivatives.

}
\usage{
logit(theta, earg = list(), inverse = FALSE, deriv = 0,
      short = TRUE, tag = FALSE)
elogit(theta, earg = list(min=0, max=1), inverse = FALSE, deriv = 0,
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
  Similarly, for \code{elogit}, values of \code{theta} less than or equal
  to \eqn{A} or greater than or equal to \eqn{B} can be replaced 
  by the \code{bminvalue} and \code{bmaxvalue} components of the list \code{earg}. 

% Extra argument for passing in additional information.
% For \code{logit}, values of \code{theta} which are equal to 0 or 1 are
% replaced by \code{earg} or \code{1-earg}
% (respectively, and if given) before computing the logit.

  For \code{elogit}, \code{earg} should be a list with components
  \code{min} giving \eqn{A}, 
  \code{max} giving \eqn{B}, and for out of range values, 
  \code{bminvalue} and \code{bmaxvalue}.
  If \code{earg} is used, these
  component names should not be abbreviated.

  }
  \item{inverse}{
  Logical. If \code{TRUE} the inverse function is computed.
  The inverse logit function is known as the \emph{expit} function.

  }
  \item{deriv}{
  Order of the derivative. Integer with value 0, 1 or 2.

  }
  \item{short}{
  Used for labelling the \code{blurb} slot of a
  \code{\link{vglmff-class}} object.

  }
  \item{tag}{
  Used for labelling the linear/additive predictor in the
  \code{initialize} slot of a \code{\link{vglmff-class}} object.
  Contains a little more information if \code{TRUE}.

  }
}
\details{
  The logit link function is very commonly used for parameters that
  lie in the unit interval. 
  Numerical values of \code{theta} close to 0 or 1 or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.

  The \emph{extended} logit link function \code{elogit} should be used
  more generally for parameters that lie in the interval \eqn{(A,B)}, say.
  The formula is
  \deqn{\log((\theta-A)/(B-\theta))}{%
        log((theta-A)/(B-theta))}
  and the default values for \eqn{A} and \eqn{B} correspond to the ordinary
  logit function.
  Numerical values of \code{theta} close to \eqn{A} or \eqn{B} or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}. 
  However these can be replaced by values \eqn{bminvalue} and
  \eqn{bmaxvalue} first before computing the link function.

  The arguments \code{short} and \code{tag} are used only if
  \code{theta} is character.

}
\value{
  For \code{logit} with \code{deriv = 0}, the logit of \code{theta}, i.e.,
  \code{log(theta/(1-theta))} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{exp(theta)/(1+exp(theta))}.

  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.

  Here, all logarithms are natural logarithms, i.e., to base \emph{e}.

}
\references{
  McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.

}
\author{ Thomas W. Yee }

\note{
  Numerical instability may occur when \code{theta} is
  close to 1 or 0 (for \code{logit}), or close to \eqn{A} or \eqn{B} for
  \code{elogit}.
  One way of overcoming this is to use \code{earg}.

  In terms of the threshold approach with cumulative probabilities for
  an ordinal response this link function corresponds to the univariate
  logistic distribution (see \code{\link{logistic}}).

}

\seealso{ 
    \code{\link{Links}},
    \code{\link{probit}},
    \code{\link{cloglog}},
    \code{\link{cauchit}},
    \code{\link{logistic1}},
    \code{\link{loge}}.

 }
\examples{
p = seq(0.01, 0.99, by=0.01)
logit(p)
max(abs(logit(logit(p), inverse=TRUE) - p)) # Should be 0

p = c(seq(-0.02, 0.02, by=0.01), seq(0.97, 1.02, by=0.01))
logit(p)  # Has NAs
logit(p, earg=list(bvalue= .Machine$double.eps))  # Has no NAs

p = seq(0.9, 2.2, by=0.1)
elogit(p, earg=list(min=1, max=2,
                    bminvalue = 1 + .Machine$double.eps,
                    bmaxvalue = 2 - .Machine$double.eps))  # Has no NAs

\dontrun{
par(mfrow=c(2,2))
y = seq(-4, 4, length=100)
for(d in 0:1) {
    matplot(p, cbind(logit(p, deriv=d), probit(p, deriv=d)),
            type="n", col="purple", ylab="transformation",
            lwd=2, las=1,
            main = if (d == 0) "Some probability link functions"
            else "First derivative")
    lines(p, logit(p, deriv=d), col="limegreen", lwd=2)
    lines(p, probit(p, deriv=d), col="purple", lwd=2)
    lines(p, cloglog(p, deriv=d), col="chocolate", lwd=2)
    lines(p, cauchit(p, deriv=d), col="tan", lwd=2)
    if (d == 0) {
        abline(v=0.5, h=0, lty="dashed")
        legend(0, 4.5, c("logit", "probit", "cloglog", "cauchit"),
               col=c("limegreen","purple","chocolate", "tan"), lwd=2)
    } else
        abline(v=0.5, lty="dashed")
}

for(d in 0) {
    matplot(y, cbind(logit(y, deriv=d, inverse=TRUE),
                     probit(y, deriv=d, inverse=TRUE)),
            type="n", col="purple", xlab="transformation", ylab="p",
            lwd=2, las=1,
            main = if (d == 0) "Some inverse probability link functions"
            else "First derivative")
    lines(y, logit(y, deriv=d, inverse=TRUE), col="limegreen", lwd=2)
    lines(y, probit(y, deriv=d, inverse=TRUE), col="purple", lwd=2)
    lines(y, cloglog(y, deriv=d, inverse=TRUE), col="chocolate", lwd=2)
    lines(y, cauchit(y, deriv=d, inverse=TRUE), col="tan", lwd=2)
    if (d == 0) {
        abline(h=0.5, v=0, lty="dashed")
        legend(-4, 1, c("logit", "probit", "cloglog", "cauchit"),
               col=c("limegreen","purple","chocolate", "tan"), lwd=2)
    }
}

p = seq(0.21, 0.59, by=0.01)
plot(p, elogit(p, earg=list(min=0.2, max=0.6)), lwd=2, 
     type="l", col="black", ylab="transformation", xlim=c(0,1),
     las=1, main="elogit(p, earg=list(min=0.2, max=0.6)")
}
}
\keyword{math}
\keyword{models}
\keyword{regression}

%plot(y, logit(y, inverse=TRUE), type="l", col="limegreen",
%     xlab="transformation", ylab="p",
%     lwd=2, las=1, main="Some inverse probability link functions")
%lines(y, probit(y, inverse=TRUE), col="purple", lwd=2)
%lines(y, cloglog(y, inverse=TRUE), col="chocolate", lwd=2)
%abline(h=0.5, v=0, lty="dashed")




