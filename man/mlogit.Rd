\name{mlogit}
\alias{mlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multi-logit Link Function }
\description{
  Computes the mlogit transformation, including its inverse and the
  first two derivatives.

}
\usage{
mlogit(theta, refLevel = "last", M = NULL, whitespace = FALSE,
       bvalue = NULL, inverse = FALSE, deriv = 0,
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

  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The \code{mlogit()} link function is a generalization of the
  \code{\link{logit}} link to \eqn{M} levels/classes.
  It forms the basis of the \code{\link{multinomial}} logit model.
  It is sometimes called the \emph{multi-logit} link
  or the \emph{multinomial logit} link.
  When its inverse function is computed it returns values which
  are positive and add to unity.


}
\value{
  For \code{mlogit} with \code{deriv = 0}, the mlogit of \code{theta},
  i.e.,
  \code{log(theta[, j]/theta[, M+1])} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{exp(theta[, j])/(1+rowSums(exp(theta)))}.


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
  close to 1 or 0 (for \code{mlogit}).
  One way of overcoming this is to use, e.g., \code{bvalue}.
  Currently \code{care.exp()} is used to avoid \code{NA}s being
  returned if the probability is too close to 1.


}

\seealso{ 
    \code{\link{Links}},
    \code{\link{multinomial}},
    \code{\link{logit}},
    \code{\link{normal.vcm}},
    \code{\link{CommonVGAMffArguments}}.


 }
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let,
            multinomial, trace = TRUE, pneumo)  # For illustration only!
fitted(fit)
predict(fit)

mlogit(fitted(fit))
mlogit(fitted(fit)) - predict(fit)  # Should be all 0s

mlogit(predict(fit), inverse = TRUE)  # rowSums() add to unity
mlogit(predict(fit), inverse = TRUE, refLevel = 1)  # For illustration only
mlogit(predict(fit), inverse = TRUE) - fitted(fit)  # Should be all 0s

mlogit(fitted(fit), deriv = 1)
mlogit(fitted(fit), deriv = 2)
}
\keyword{math}
\keyword{models}
\keyword{regression}

