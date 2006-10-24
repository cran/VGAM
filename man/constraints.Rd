\name{constraints}
\alias{constraints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Constraint Matrices }
\description{
  Returns the \emph{constraint matrices} of objects in the
  \pkg{VGAM} package.

}
\usage{
constraints(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Some \pkg{VGAM} object, for example, having
  class \code{\link{vglmff-class}}.
  }
  \item{\dots}{ Other possible arguments. }
}
\details{
  Constraint matrices describe the relationship of
  coefficients/component functions of a particular explanatory
  variable between
  the linear/additive predictors in VGLM/VGAM
  models. For example, they may be all different (constraint
  matrix is the identity matrix) or all the same (constraint
  matrix has one column and has unit values).

  VGLMs and VGAMs have constraint matrices which are \emph{known}.
  The class of RR-VGLMs have constraint matrices which are
  \emph{unknown} and are to be estimated.
}
\value{
  This extractor function returns a list comprising of
  constraint matrices---one for each column of the
  LM model matrix, and in that order.
  The list is labelled with the variable names.
  Each constraint matrix has \eqn{M} rows, where
  \eqn{M} is the number of linear/additive predictors,
  and whose rank is equal to the number of columns.
  A model with no constraints at all has an order
  \eqn{M} identity matrix as each variable's
  constraint matrix.
}
\author{T. W. Yee }
\note{
  In all \pkg{VGAM} family functions \code{zero=NULL} means
  none of the linear/additive predictors are modelled as
  intercepts-only.
  Other arguments found in certain \pkg{VGAM} family functions
  which affect constraint matrices include
  \code{parallel} and \code{exchangeable}.

  The \code{constraints} argument in \code{\link{vglm}}
  and \code{\link{vgam}} allows constraint matrices to
  be inputted. If so, then \code{constraints(fit)} should
  return the same as the input.
}

\references{
Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.

Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.

\url{http://www.stat.auckland.ac.nz/~yee} contains additional
information.
}

\section{Warning }{The \code{xij} argument changes things,
  and this has not been fully resolved yet.
}

\seealso{
  VGLMs are described in \code{\link{vglm-class}};
  RR-VGLMs are described in \code{\link{rrvglm-class}}.

  Arguments such as \code{zero} and \code{parallel}
  found in many \pkg{VGAM}
  family functions are a way of creating/modifying constraint
  matrices conveniently, e.g., see \code{\link{zero}}.
}
\examples{
# Fit the proportional odds model
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
fit = vglm(cbind(normal, mild, severe) ~ let,
           cumulative(parallel=TRUE, reverse=TRUE), pneumo)
fit
coef(fit, matrix=TRUE)
constraints(fit)  # Parallel assumption results in this




# Fit a rank-1 stereotype (RR-multinomial logit) model 
data(car.all)
fit = rrvglm(Country ~ Width + Height + HP, multinomial, car.all, Rank=1)
constraints(fit) # All except the first are the A matrix 
}
\keyword{models}
\keyword{regression}

