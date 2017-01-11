\name{vglm-class}
\docType{class}
\alias{vglm-class}
\title{Class ``vglm'' }
\description{  Vector generalized linear models. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{vglm(...)}.
%   ~~ describe objects here ~~
}
\section{Slots}{

In the following, \eqn{M} is the number of linear predictors.

  \describe{
    \item{\code{extra}:}{Object of class \code{"list"};
the \code{extra} argument on entry to \code{vglm}. This
contains any extra information that might be needed
    by the family function. }
    \item{\code{family}:}{Object of class \code{"vglmff"}.
The family function.  }
    \item{\code{iter}:}{Object of class \code{"numeric"}.
The number of IRLS iterations used. }
    \item{\code{predictors}:}{Object of class \code{"matrix"}
with \eqn{M} columns which holds the \eqn{M} linear predictors. }
    \item{\code{assign}:}{Object of class \code{"list"},
  from class \code{ "vlm"}.
  This named list gives information matching the columns and the
  (LM) model matrix terms.
    }
    \item{\code{call}:}{Object of class \code{"call"}, from class
  \code{ "vlm"}.
  The matched call.
    }
    \item{\code{coefficients}:}{Object of class
\code{"numeric"}, from class \code{ "vlm"}.
A named vector of coefficients.
    }
    \item{\code{constraints}:}{Object of class \code{"list"}, from
  class \code{ "vlm"}.
  A named list of constraint matrices used in the fitting.
    }
    \item{\code{contrasts}:}{Object of class \code{"list"}, from
  class \code{ "vlm"}.
  The contrasts used (if any).
    }
    \item{\code{control}:}{Object of class \code{"list"}, from class
  \code{ "vlm"}.
  A list of parameters for controlling the fitting process.
          See \code{\link{vglm.control}} for details.
    }
    \item{\code{criterion}:}{Object of class \code{"list"}, from
  class \code{ "vlm"}.
  List of convergence criterion evaluated at the
                   final IRLS iteration.
    }
    \item{\code{df.residual}:}{Object of class
\code{"numeric"}, from class \code{ "vlm"}.
The residual degrees of freedom.
    }
    \item{\code{df.total}:}{Object of class \code{"numeric"},
 from class \code{ "vlm"}.
The total degrees of freedom.
    }
    \item{\code{dispersion}:}{Object of class \code{"numeric"},
 from class \code{ "vlm"}.
The scaling parameter.
 }
    \item{\code{effects}:}{Object of class \code{"numeric"},
 from class \code{ "vlm"}.
The effects.
 }
    \item{\code{fitted.values}:}{Object of class
\code{"matrix"}, from class \code{ "vlm"}.
The fitted values.
%This may be missing or consist entirely
%of \code{NA}s, e.g., the Cauchy model.
 }
    \item{\code{misc}:}{Object of class \code{"list"},
 from class \code{ "vlm"}.
A named list to hold miscellaneous parameters.
 }
    \item{\code{model}:}{Object of class \code{"data.frame"},
 from class \code{ "vlm"}.
The model frame.
 }
    \item{\code{na.action}:}{Object of class \code{"list"},
 from class \code{ "vlm"}.
A list holding information about missing values.
 }
    \item{\code{offset}:}{Object of class \code{"matrix"},
 from class \code{ "vlm"}.
If non-zero, a \eqn{M}-column matrix of offsets.
 }
    \item{\code{post}:}{Object of class \code{"list"},
 from class \code{ "vlm"}
where post-analysis results may be put.
 }
    \item{\code{preplot}:}{Object of class \code{"list"},
 from class \code{ "vlm"}
used by \code{\link{plotvgam}}; the plotting parameters
        may be put here.
 }
    \item{\code{prior.weights}:}{Object of class
\code{"matrix"}, from class \code{ "vlm"}
holding the initially supplied weights.
 }
    \item{\code{qr}:}{Object of class \code{"list"},
 from class \code{ "vlm"}.
QR decomposition at the final iteration.
 }
    \item{\code{R}:}{Object of class \code{"matrix"},
 from class \code{ "vlm"}.
The \bold{R} matrix in the QR decomposition used in the fitting.
 }
    \item{\code{rank}:}{Object of class \code{"integer"},
 from class \code{ "vlm"}.
Numerical rank of the fitted model.
 }
    \item{\code{residuals}:}{Object of class \code{"matrix"},
 from class \code{ "vlm"}.
The \emph{working} residuals at the final IRLS iteration.
 }
    \item{\code{ResSS}:}{Object of class \code{"numeric"},
 from class \code{ "vlm"}.
Residual sum of squares at the final IRLS iteration with
the adjusted dependent vectors and weight matrices.
 }
    \item{\code{smart.prediction}:}{Object of class
\code{"list"}, from class \code{ "vlm"}.
A list of data-dependent parameters (if any)
    that are used by smart prediction.
 }
    \item{\code{terms}:}{Object of class \code{"list"},
 from class \code{ "vlm"}.
The \code{\link[stats]{terms}} object used.
 }
    \item{\code{weights}:}{Object of class \code{"matrix"},
 from class \code{ "vlm"}.
The weight matrices at the final IRLS iteration.
    This is in matrix-band form.
 }
    \item{\code{x}:}{Object of class \code{"matrix"},
 from class \code{ "vlm"}.
The model matrix (LM, not VGLM).
 }
    \item{\code{xlevels}:}{Object of class \code{"list"},
 from class \code{ "vlm"}.
The levels of the factors, if any, used in fitting.
 }
    \item{\code{y}:}{Object of class \code{"matrix"},
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
Class \code{"vlm"}, directly.
}
\section{Methods}{
  \describe{
    \item{cdf}{\code{signature(object = "vglm")}:
cumulative distribution function.
Applicable to, e.g., quantile regression and extreme value data models.}
    \item{deplot}{\code{signature(object = "vglm")}:
Applicable to, e.g., quantile regression.}
    \item{deviance}{\code{signature(object = "vglm")}:
deviance of the model (where applicable). }
    \item{plot}{\code{signature(x = "vglm")}:
diagnostic plots. }
    \item{predict}{\code{signature(object = "vglm")}:
extract the linear predictors or
predict the linear predictors at a new data frame.}
    \item{print}{\code{signature(x = "vglm")}:
short summary of the object. }
    \item{qtplot}{\code{signature(object = "vglm")}:
quantile plot (only applicable to some models). }
    \item{resid}{\code{signature(object = "vglm")}:
residuals. There are various types of these. }
    \item{residuals}{\code{signature(object = "vglm")}:
residuals. Shorthand for \code{resid}. }
    \item{rlplot}{\code{signature(object = "vglm")}: return level plot.
Useful for extreme value data models.}
    \item{summary}{\code{signature(object = "vglm")}:
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
%\note{ ~~further notes~~ }

%~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{
  \code{\link{vglm}},
  \code{\link{vglmff-class}},
  \code{\link{vgam-class}}.



}

\examples{
# Multinomial logit model
pneumo <- transform(pneumo, let = log(exposure.time))
vglm(cbind(normal, mild, severe) ~ let, multinomial, data = pneumo)
}
\keyword{classes}
