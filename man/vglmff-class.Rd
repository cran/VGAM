\name{vglmff-class}
\docType{class}
\alias{vglmff-class}
\title{Class ``vglmff'' }
\description{  Family functions for the \pkg{VGAM} package }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("vglmff", ...)}.

}
\section{Slots}{
  In the following, \eqn{M} is the number of linear/additive
  predictors.

\describe{
  \item{\code{blurb}:}{
  Object of class \code{"character"} giving
  a small description of the model. Important arguments such as
  parameter link functions can be expressed here.

  }
  \item{\code{constraints}:}{
  Object of class \code{"expression"}
  which sets up any constraint matrices defined by arguments in the
  family function. A \code{zero} argument is always fed into
  \code{cm.zero.vgam}, whereas other constraints are fed into
  \code{cm.vgam}.

  }
  \item{\code{deviance}:}{
  Object of class \code{"function"}
  returning the deviance of the model. This slot is optional.
  If present, the function must have arguments
  \code{function(mu, y, w, residuals = FALSE, eta, extra = NULL)}.
  Deviance residuals are returned if \code{residuals = TRUE}.

  }
  \item{\code{fini}:}{
  Object of class \code{"expression"} to insert
  code at a special position in \code{vglm.fit} or
  \code{vgam.fit}.
  This code is evaluated immediately after the fitting.

  }
  \item{\code{first}:}{
  Object of class \code{"expression"} to insert
  code at a special position in \code{\link{vglm}} or
  \code{\link{vgam}}.

  }
  \item{\code{infos}:}{
  Object of class \code{"function"} which
  returns a list with components such as \code{M1}.
  At present only a very few \pkg{VGAM} family functions have this
  feature implemented.
  Those that do do not require specifying the \code{M1}
  argument when used with \code{\link{rcim}}.


  }
  \item{\code{initialize}:}{
  Object of class \code{"expression"} used
  to perform error checking (especially for the variable \code{y})
  and obtain starting values for the model.
  In general, \code{etastart} or
  \code{mustart} are assigned values based on the variables \code{y},
  \code{x} and \code{w}.

  }
  \item{\code{linkinv}:}{
  Object of class \code{"function"} which
  returns the fitted values, given the linear/additive predictors.
  The function must have arguments
  \code{function(eta, extra = NULL)}.

  }
  \item{\code{last}:}{
  Object of class \code{"expression"} to insert code at a
  special position (at the very end) of \code{vglm.fit()}
  or \code{vgam.fit()}.
  This code is evaluated after the fitting.
  The list \code{misc} is often assigned components in this slot,
  which becomes the \code{misc} slot on the fitted object.

  }
  \item{\code{linkfun}:}{
  Object of class \code{"function"} which,
  given the fitted values, returns the linear/additive predictors.
  If present, the function must have arguments
  \code{function(mu, extra = NULL)}.
  Most \pkg{VGAM} family functions do not have
  a \code{linkfun} function. They largely are for
  classical exponential families, i.e., GLMs.


  }
  \item{\code{loglikelihood}:}{
  Object of class \code{"function"}
  returning the log-likelihood of the model. This slot is optional.
  If present, the function must have arguments
  \code{function(mu, y, w, residuals = FALSE, eta, extra = NULL)}.
  The argument \code{residuals} can be ignored because
  log-likelihood residuals aren't defined.


  }
  \item{\code{middle}:}{
  Object of class \code{"expression"} to insert
  code at a special position in \code{vglm.fit} or
  \code{vgam.fit}.


  }
  \item{\code{middle2}:}{
  Object of class \code{"expression"} to insert
  code at a special position in \code{vglm.fit} or
  \code{vgam.fit}.


  }
  \item{\code{simslot}:}{
  Object of class \code{"function"} to allow
  \code{\link[stats]{simulate}} to work.



  }
  \item{\code{summary.dispersion}:}{
  Object of class \code{"logical"}
  indicating whether the general VGLM formula (based on a residual
  sum of squares) can be used for computing the scaling/dispersion
  parameter. It is \code{TRUE} for most models except for nonlinear
  regression models.


  }
  \item{\code{vfamily}:}{
  Object of class \code{"character"}
  giving class information about the family function. Although
  not developed at this stage, more flexible classes are planned
  in the future. For example, family functions
  \code{sratio}, \code{cratio}, \code{cumulative}, and \code{acat}
  all operate on categorical data, therefore will have a special class
  called \code{"VGAMcat"}, say. Then if \code{fit} was
  a \code{vglm} object, then \code{coef(fit)} would print
  out the \code{vglm} coefficients plus \code{"VGAMcat"}
  information as well.


  }
  \item{\code{deriv}:}{
  Object of class \code{"expression"} which
  returns a \eqn{M}-column matrix of first derivatives of the
  log-likelihood function
  with respect to the linear/additive predictors, i.e., the
  score vector. In Yee and Wild (1996) this is the
  \eqn{\bold{d}_i}{\bold{d}i} vector. Thus each row of the
  matrix returned by this slot is such a vector.


  }
  \item{\code{weight}:}{
  Object of class \code{"expression"} which
  returns the second derivatives of the log-likelihood function
  with respect to the linear/additive predictors.
  This can be either the observed or expected information matrix, i.e.,
  Newton-Raphson or Fisher-scoring respectively.
  In Yee and Wild (1996) this is the
  \eqn{\bold{W}_i}{\bold{W}i} matrix. Thus each row of the
  matrix returned by this slot is such a matrix.
  Like the \code{weights} slot of \code{vglm}/\code{vgam}, it is
  stored in
  \emph{matrix-band} form, whereby the first \eqn{M}
  columns of the matrix are the
  diagonals, followed by the upper-diagonal band, followed by the
  band above that, etc. In this case, there can be up to \eqn{M(M+1)}
  columns, with the last column corresponding to the (1,\eqn{M}) elements
  of the weight matrices.

  }

  \item{\code{validfitted, validparams}:}{
  Functions that test that the fitted values and
  all parameters are within range.
  These functions can issue a warning if violations are detected.


  }


}
}

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "vglmff")}:
      short summary of the family function.
    }
  }
}
\references{

Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.


%\url{http://www.stat.auckland.ac.nz/~yee} contains further
%information on how to write \pkg{VGAM} family functions.
%The file is amongst other \pkg{VGAM} PDF documentation.


}
\author{ Thomas W. Yee }
\note{
  With link functions etc., one must use \code{substitute} to
  embed the options into the code. There are two different forms:
  \code{eval(substitute(expression({...}), list(...)))}
  for expressions, and
  \code{eval(substitute( function(...) { ... }, list(...) )) }
  for functions.


% 20130322; this is obsolete, and can delete it:
% A unified method of handling arguments is to use
% \code{match.arg}. This allows, for example,
% \code{vglm(..., family = cratio(link = logit))}
% and
% \code{vglm(..., family = cratio(link = "logi"))}
% to be equivalent (Nb. there is a \code{logit} function).


  The \code{extra} argument in
  \code{linkinv}, \code{linkfun}, \code{deviance},
  \code{loglikelihood}, etc.
  matches with the argument \code{extra}
  in \code{\link{vglm}}, \code{\link{vgam}} and \code{\link{rrvglm}}.
  This allows input to be fed into all slots of a \pkg{VGAM}
  family function.


  The expression \code{derivative} is evaluated immediately
  prior to \code{weight}, so there is provision for re-use
  of variables etc.  Programmers must be careful to choose
  variable names that do not interfere with \code{vglm.fit},
  \code{vgam.fit()} etc.


  Programmers of \pkg{VGAM} family functions are encouraged
  to keep to previous conventions regarding the naming of arguments,
  e.g.,
  \code{link} is the argument for parameter link functions,
  \code{zero} for allowing some of the
  linear/additive predictors to be an intercept term only, etc.


  In general, Fisher-scoring is recommended over
  Newton-Raphson where tractable. Although usually slightly
  slower in convergence, the weight matrices from using the
  expected information are positive-definite over a larger
  parameter space.


  }

\section{Warning }{
  \pkg{VGAM} family functions are not compatible with
  \code{\link[stats]{glm}}, nor \code{gam} (from either
  \pkg{gam} or \pkg{mgcv} packages).


}

\seealso{
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{rrvglm}},
  \code{\link{rcim}}.


}
\examples{
cratio()
cratio(link = "cloglog")
cratio(link = "cloglog", reverse = TRUE)
}
\keyword{classes}






