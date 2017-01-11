\name{constraints}
\alias{constraints}
\alias{constraints.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Constraint Matrices }
\description{
  Extractor function for the \emph{constraint matrices} of objects
  in the \pkg{VGAM} package.

}
\usage{
constraints(object, ...)
constraints.vlm(object, type = c("lm", "term"), all = TRUE, which,
                matrix.out = FALSE, colnames.arg = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Some \pkg{VGAM} object, for example, having
  class \code{\link{vglmff-class}}.


  }
  \item{type}{
  Character. Whether LM- or term-type constraints are to be returned.
  The number of such matrices returned is equal to
  \code{nvar(object, type = "lm")} and
  the number of terms, respectively.


  }
  \item{all, which}{
  If \code{all = FALSE} then \code{which} gives the integer index or a
  vector of logicals specifying the selection.


  }
  \item{matrix.out}{
  Logical. If \code{TRUE} then the constraint matrices
  are \code{\link[base]{cbind}()ed} together.
  The result is usually more compact because the default
  is a list of constraint matrices.


  }
  \item{colnames.arg}{
  Logical. If \code{TRUE} then column names are assigned
  corresponding to the variables.


  }
  \item{\dots}{
  Other possible arguments such as \code{type}.


  }

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
  The extractor function
  \code{constraints()}
  returns a list comprising of
  constraint matrices---usually one for each column of the
  VLM model matrix, and in that order.
  The list is labelled with the variable names.
  Each constraint matrix has \eqn{M} rows, where
  \eqn{M} is the number of linear/additive predictors,
  and whose rank is equal to the number of columns.
  A model with no constraints at all has an order
  \eqn{M} identity matrix as each variable's
  constraint matrix.


  For \code{\link{vglm}} and \code{\link{vgam}} objects,
  feeding in \code{type = "term"} constraint matrices back
  into the same model should work and give an identical model.
  The default are the \code{"lm"}-type constraint matrices;
  this is a list with one constraint matrix per column of
  the LM matrix.
  See the \code{constraints} argument of \code{\link{vglm}},
  and the example below.


}
\author{T. W. Yee }
\note{
  In all \pkg{VGAM} family functions \code{zero = NULL} means
  none of the linear/additive predictors are modelled as
  intercepts-only.
  Other arguments found in certain \pkg{VGAM} family functions
  which affect constraint matrices include
  \code{parallel} and \code{exchangeable}.


  The \code{constraints} argument in \code{\link{vglm}}
  and \code{\link{vgam}} allows constraint matrices to
  be inputted. If so, then \code{constraints(fit, type = "lm")} can
  be fed into the \code{constraints} argument of the same object
  to get the same model.


The \code{xij} argument does not affect constraint matrices; rather,
it allows each row of the constraint matrix to be multiplied by a
specified vector.


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


%\url{http://www.stat.auckland.ac.nz/~yee} contains additional
%information.


}


\seealso{
  \code{\link{is.parallel}},
  \code{\link{is.zero}}.
  VGLMs are described in \code{\link{vglm-class}};
  RR-VGLMs are described in \code{\link{rrvglm-class}}.


  Arguments such as \code{zero} and \code{parallel}
  found in many \pkg{VGAM}
  family functions are a way of creating/modifying constraint
  matrices conveniently, e.g., see \code{\link{zero}}.
  See \code{\link{CommonVGAMffArguments}} for more information.


}
\examples{
# Fit the proportional odds model:
pneumo <- transform(pneumo, let = log(exposure.time))
(fit1 <- vglm(cbind(normal, mild, severe) ~ sm.bs(let, 3),
              cumulative(parallel = TRUE, reverse = TRUE), data = pneumo))
coef(fit1, matrix = TRUE)
constraints(fit1)  # Parallel assumption results in this
constraints(fit1, type = "term")  # Same as the default ("vlm"-type)
is.parallel(fit1)

# An equivalent model to fit1 (needs the type "term" constraints):
clist.term <- constraints(fit1, type = "term")  # "term"-type constraints
(fit2 <- vglm(cbind(normal, mild, severe) ~ sm.bs(let, 3), data = pneumo,
              cumulative(reverse = TRUE), constraints = clist.term))
abs(max(coef(fit1, matrix = TRUE) -
        coef(fit2, matrix = TRUE)))  # Should be zero

# Fit a rank-1 stereotype (RR-multinomial logit) model:
fit <- rrvglm(Country ~ Width + Height + HP, multinomial, data = car.all)
constraints(fit)  # All except the first are the estimated A matrix
}
\keyword{models}
\keyword{regression}

