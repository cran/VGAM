\name{multilogitlink}
\alias{multilogitlink}
% \alias{multilogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multi-logit Link Function }
\description{
  Computes the multilogit transformation, including its inverse and the
  first two derivatives.

}
\usage{
multilogitlink(theta, refLevel = "(Last)", M = NULL, whitespace = FALSE,
               bvalue = NULL, inverse = FALSE, deriv = 0, all.derivs = FALSE,
               short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{refLevel, M, whitespace}{
  See \code{\link{multinomial}}.


  }
  \item{bvalue}{
  See \code{\link{Links}}.


  }

  \item{all.derivs}{
  Logical.
  This is currently experimental only.



% If \code{TRUE} then more partial derivatives are returned;
% these is needed by, e.g.,
% \code{\link{hdeff.vglm}} for \code{\link{multinomial}} fits.
% This argument might work for only some combinations of
% the arguments, e.g., it should work at least for
% \code{inverse = TRUE} and \code{deriv = 1}.



  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The \code{multilogitlink()} link function is a generalization of the
  \code{\link{logitlink}} link to \eqn{M} levels/classes.
  It forms the basis of the \code{\link{multinomial}} logit model.
  It is sometimes called the \emph{multi-logit} link
  or the \emph{multinomial logit} link.
  When its inverse function is computed it returns values which
  are positive and add to unity.


}
\value{
  For \code{multilogitlink} with \code{deriv = 0}, the multilogit of \code{theta},
  i.e.,
  \code{log(theta[, j]/theta[, M+1])} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{exp(theta[, j])/(1+rowSums(exp(theta)))}.


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
  close to 1 or 0 (for \code{multilogitlink}).
  One way of overcoming this is to use, e.g., \code{bvalue}.
  Currently \code{care.exp()} is used to avoid \code{NA}s being
  returned if the probability is too close to 1.


}

\seealso{
    \code{\link{Links}},
    \code{\link{multinomial}},
    \code{\link{logitlink}},
    \code{\link{gaitpoisson.mlm}},
    \code{\link{normal.vcm}},
    \code{\link{CommonVGAMffArguments}}.


 }
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let,
            multinomial, trace = TRUE, data = pneumo)  # For illustration only!
fitted(fit)
predict(fit)

multilogitlink(fitted(fit))
multilogitlink(fitted(fit)) - predict(fit)  # Should be all 0s

multilogitlink(predict(fit), inverse = TRUE)  # rowSums() add to unity
multilogitlink(predict(fit), inverse = TRUE, refLevel = 1)  # For illustration only
multilogitlink(predict(fit), inverse = TRUE) - fitted(fit)  # Should be all 0s

multilogitlink(fitted(fit), deriv = 1)
multilogitlink(fitted(fit), deriv = 2)
}
\keyword{math}
\keyword{models}
\keyword{regression}

