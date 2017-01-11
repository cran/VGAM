\name{rrvglm-class}
\docType{class}
\alias{rrvglm-class}
\title{Class ``rrvglm'' }
\description{
Reduced-rank vector generalized linear models.
}
\section{Objects from the Class}{
Objects can be created by calls to \code{\link{rrvglm}}.
}
\section{Slots}{
  \describe{
    \item{\code{extra}:}{
  Object of class \code{"list"};
  the \code{extra} argument on entry to \code{vglm}. This
  contains any extra information that might be needed
  by the family function.
  }
  \item{\code{family}:}{
  Object of class \code{"vglmff"}.
  The family function.  }
    \item{\code{iter}:}{
  Object of class \code{"numeric"}.
  The number of IRLS iterations used.
  }
  \item{\code{predictors}:}{
  Object of class \code{"matrix"}
  with \eqn{M} columns which holds the \eqn{M} linear predictors.
  }
  \item{\code{assign}:}{
  Object of class \code{"list"},
  from class \code{ "vlm"}.
  This named list gives information matching the columns and the
  (LM) model matrix terms.
  }
  \item{\code{call}:}{
  Object of class \code{"call"}, from class \code{ "vlm"}.
  The matched call.
  }
  \item{\code{coefficients}:}{
  Object of class
  \code{"numeric"}, from class \code{ "vlm"}.
  A named vector of coefficients.
  }
  \item{\code{constraints}:}{
  Object of class \code{"list"}, from
  class \code{ "vlm"}.
  A named list of constraint matrices used in the fitting.
  }
  \item{\code{contrasts}:}{
  Object of class \code{"list"}, from
  class \code{ "vlm"}.
  The contrasts used (if any).
  }
  \item{\code{control}:}{
  Object of class \code{"list"}, from class
  \code{ "vlm"}.
  A list of parameters for controlling the fitting process.
          See \code{\link{vglm.control}} for details.
  }
  \item{\code{criterion}:}{
  Object of class \code{"list"}, from
  class \code{ "vlm"}.
  List of convergence criterion evaluated at the
                   final IRLS iteration.
  }
  \item{\code{df.residual}:}{
  Object of class
  \code{"numeric"}, from class \code{ "vlm"}.
  The residual degrees of freedom.
  }
  \item{\code{df.total}:}{
  Object of class \code{"numeric"},
  from class \code{ "vlm"}.
  The total degrees of freedom.
  }
  \item{\code{dispersion}:}{
  Object of class \code{"numeric"},
  from class \code{ "vlm"}.
  The scaling parameter.
  }
  \item{\code{effects}:}{
  Object of class \code{"numeric"},
  from class \code{ "vlm"}.
  The effects.
  }
  \item{\code{fitted.values}:}{
  Object of class
  \code{"matrix"}, from class \code{ "vlm"}.
  The fitted values. This is usually the mean but may be
  quantiles, or the location parameter, e.g., in the Cauchy model.

  }
  \item{\code{misc}:}{
  Object of class \code{"list"},
  from class \code{ "vlm"}.
  A named list to hold miscellaneous parameters.
  }
  \item{\code{model}:}{
  Object of class \code{"data.frame"},
  from class \code{ "vlm"}.
  The model frame.
  }
  \item{\code{na.action}:}{
  Object of class \code{"list"},
  from class \code{ "vlm"}.
  A list holding information about missing values.
  }
  \item{\code{offset}:}{
  Object of class \code{"matrix"},
  from class \code{ "vlm"}.
  If non-zero, a \eqn{M}-column matrix of offsets.
  }
  \item{\code{post}:}{
  Object of class \code{"list"},
  from class \code{ "vlm"}
  where post-analysis results may be put.
  }
  \item{\code{preplot}:}{
  Object of class \code{"list"},
  from class \code{ "vlm"}
  used by \code{\link{plotvgam}}; the plotting parameters
  may be put here.
  }
  \item{\code{prior.weights}:}{
  Object of class
  \code{"matrix"}, from class \code{ "vlm"}
  holding the initially supplied weights.
  }
  \item{\code{qr}:}{
  Object of class \code{"list"},
  from class \code{ "vlm"}.
  QR decomposition at the final iteration.
  }
  \item{\code{R}:}{
  Object of class \code{"matrix"},
  from class \code{ "vlm"}.
  The \bold{R} matrix in the QR decomposition used in the fitting.
  }
  \item{\code{rank}:}{
  Object of class \code{"integer"},
  from class \code{ "vlm"}.
  Numerical rank of the fitted model.
  }
  \item{\code{residuals}:}{
  Object of class \code{"matrix"},
  from class \code{ "vlm"}.
  The \emph{working} residuals at the final IRLS iteration.
  }
  \item{\code{ResSS}:}{
  Object of class \code{"numeric"},
  from class \code{ "vlm"}.
  Residual sum of squares at the final IRLS iteration with
  the adjusted dependent vectors and weight matrices.
  }
  \item{\code{smart.prediction}:}{
  Object of class
  \code{"list"}, from class \code{ "vlm"}.
  A list of data-dependent parameters (if any)
  that are used by smart prediction.
  }
  \item{\code{terms}:}{
  Object of class \code{"list"},
  from class \code{ "vlm"}.
  The \code{\link[stats]{terms}} object used.
  }
  \item{\code{weights}:}{
  Object of class \code{"matrix"},
  from class \code{ "vlm"}.
  The weight matrices at the final IRLS iteration.
  This is in matrix-band form.
  }
  \item{\code{x}:}{
  Object of class \code{"matrix"},
  from class \code{ "vlm"}.
  The model matrix (LM, not VGLM).
  }
  \item{\code{xlevels}:}{
  Object of class \code{"list"},
  from class \code{ "vlm"}.
  The levels of the factors, if any, used in fitting.
  }
  \item{\code{y}:}{
  Object of class \code{"matrix"},
  from class \code{ "vlm"}.
  The response, in matrix form.
  }
    \item{\code{Xm2}:}{Object of class \code{"matrix"},
 from class \code{ "vlm"}.
 See \code{\link{vglm-class}}).
 }
    \item{\code{Ym2}:}{Object of class \code{"matrix"},
 from class \code{ "vlm"}.
 See \code{\link{vglm-class}}).
 }
  \item{\code{callXm2}:}{
  Object of class \code{"call"}, from class \code{ "vlm"}.
  The matched call for argument \code{form2}.
  }
  }
}

\section{Extends}{
Class \code{"vglm"}, directly.
Class \code{"vlm"}, by class "vglm".
}
\section{Methods}{
  \describe{
    \item{biplot}{\code{signature(x = "rrvglm")}: biplot. }
    \item{Coef}{\code{signature(object = "rrvglm")}: more detailed
       coefficients giving \bold{A},
       \eqn{\bold{B}_1}{\bold{B}1}, \bold{C}, etc.
    }
  \item{biplot}{\code{signature(object = "rrvglm")}:
    biplot. }
  \item{print}{\code{signature(x = "rrvglm")}:
    short summary of the object. }
  \item{summary}{\code{signature(object = "rrvglm")}:
    a more detailed summary of the object. }
  }
}
\references{
Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.


%\url{http://www.stat.auckland.ac.nz/~yee}


}
\author{ Thomas W. Yee }

\note{
  The slots of \code{"rrvglm"} objects are currently identical to \code{"vglm"}
  objects.


}

% ~Make other sections like Warning with \section{Warning }{....} ~
% zzz need to make sure this function matches \code{\link{vglm-class}},
%where \code{\link{vglm-class}} is definitive.

\seealso{
  \code{\link{rrvglm}},
  \code{\link{lvplot.rrvglm}},
  \code{\link{vglmff-class}}.


}

\examples{
\dontrun{ # Rank-1 stereotype model of Anderson (1984)
pneumo <- transform(pneumo, let = log(exposure.time),
                            x3  = runif(nrow(pneumo)))  # x3 is unrelated
fit <- rrvglm(cbind(normal, mild, severe) ~ let + x3,
              multinomial, data = pneumo, Rank = 1)
Coef(fit)
}
}
\keyword{classes}
% set.seed(111)
