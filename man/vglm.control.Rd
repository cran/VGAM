\name{vglm.control}
\alias{vglm.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Control function for vglm }
\description{
  Algorithmic constants and parameters for running \code{vglm} are set
  using this function.
}
\usage{
vglm.control(backchat = if (is.R()) FALSE else TRUE,
             checkwz=TRUE, criterion = names(.min.criterion.VGAM),
             epsilon = 1e-07, half.stepsizing = TRUE,
             maxit = 30, stepsize = 1, save.weight = FALSE,
             trace = FALSE, wzepsilon = .Machine$double.eps^0.75, 
             xij = NULL, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{backchat}{ logical indicating if a backchat is to be used (not
    applicable in \R). }
  \item{checkwz}{ logical indicating whether the diagonal elements of
    the working weight matrices should be checked whether they are
    sufficiently positive, i.e., greater than \code{wzepsilon}. If not,
    any values less than \code{wzepsilon} are replaced with this value.

  }
  \item{criterion}{ character variable describing what criterion is to
    be used to test for convergence. 
    The possibilities are listed in \code{.min.criterion.VGAM}, but
    most family functions only implement a few of these.
  }
  \item{epsilon}{ positive convergence tolerance epsilon. Roughly
    speaking, the Newton-Raphson/Fisher-scoring iterations 
    are assumed to have   
    converged when two successive \code{criterion} values are within
    \code{epsilon} of each other.
    }
    \item{half.stepsizing}{ logical indicating if half-stepsizing is
      allowed. For example, in maximizing a log-likelihood, if the
      next iteration has a log-likelihood that is less than the current
      value of the log-likelihood, then a half step will be taken.
      If the log-likelihood is still less than at the current position,
      a quarter-step will be taken etc. Eventually a step will be taken
      so that an improvement is made to the convergence criterion. 
      \code{half.stepsizing} is ignored if
      \code{criterion=="coefficients"}.
      }
      \item{maxit}{ maximum number of 
        Newton-Raphson/Fisher-scoring iterations allowed. }
      \item{stepsize}{ usual step size to be taken between each
	Newton-Raphson/Fisher-scoring iteration. It should be a value
	between 0 and 1, where 
	a value of unity corresponds to an ordinary step. 
        A value of 0.5 means half-steps are taken.
	Setting a value near zero will cause convergence 
        to be generally slow but may help increase the chances of
        successful convergence for some family functions. 
	}
  \item{save.weight}{ logical indicating whether the \code{weights} slot
  of a \code{"vglm"} object will be saved on the object. If not, it will
  be reconstructed when needed, e.g., \code{summary}.
  Some family functions have \code{save.weight=TRUE} and others have 
  \code{save.weight=FALSE} in their control functions.
}
  \item{trace}{ logical indicating if output should be produced for
  each iteration. }
  \item{wzepsilon}{
  Small positive number used to test whether the diagonals of the working
  weight matrices are sufficiently positive.

  }
  \item{xij}{
    formula giving terms making up a covariate-dependent term (a variable
    that takes on different values for each linear/additive predictor.
    For example, the ocular pressure of each eye).
    There should be \eqn{M} unique terms; use
    \code{\link{fill1}}, \code{fill2}, \code{fill3}, etc. if necessary.
    Each formula should have a response which is taken as the name of
    that variable, and the terms are enumerated in sequential order.
    With more than one formula, use a list of formulas.
    See Example 2 below.

}
  \item{\dots}{ other parameters that may be picked up from control
    functions that are specific to the \pkg{VGAM} family function. }
}
\details{
  Most of the control parameters are used within \code{vglm.fit} and you
  will have to look at that to understand the full details.

  Setting \code{save.weight=FALSE} is useful for some models because the
  \code{weights} slot of the object is the largest and so less memory
  is used to store the object. However, for some \pkg{VGAM} family function,
  it is necessary to set \code{save.weight=TRUE} because the \code{weights}
  slot cannot be reconstructed later.
  }
\value{
  A list with components matching the input names. A little error
  checking is done, but not much.
  The list is assigned to the \code{control} slot of \code{vglm} objects.
  }
  \references{
Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.

  }
\author{ Thomas W. Yee}
\note{ 
In Example 2 below there are two covariates that have linear/additive
predictor specific values.
These are handled using the \code{xij} argument.

}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{vglm}},
  \code{\link{fill}}.
}

\examples{
# Example 1.
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
vglm(cbind(normal,mild,severe) ~ let, multinomial, pneumo,
     crit="coef", step=0.5, trace=TRUE, eps=1e-8, maxit=40)


# Example 2. The use of the xij argument
set.seed(111)
n = 1000
ymat = rdiric(n, shape=c(4,7,3,1))
mydat = data.frame(x1=runif(n), x2=runif(n), x3=runif(n), x4=runif(n),
                   z1=runif(n), z2=runif(n), z3=runif(n), z4=runif(n))
mydat = round(mydat, dig=2)
fit = vglm(ymat ~ x1 + x2 + x3 + x4 + z1 + z2 + z3 + z4,
           fam = dirichlet, data=mydat, crit="c",
           xij = list(z ~ z1 + z2 + z3 + z4,
                      x ~ x1 + x2 + x3 + x4))
model.matrix(fit, type="lm")[1:7,]   # LM model matrix
model.matrix(fit, type="vlm")[1:7,]  # Big VLM model matrix
coef(fit)
coef(fit, matrix=TRUE)
coef(fit, matrix=TRUE, compress=FALSE)
max(abs(predict(fit)-predict(fit, new=mydat))) # Predicts correctly
summary(fit)
}
\keyword{models}
\keyword{regression}

