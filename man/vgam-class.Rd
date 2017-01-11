\name{vgam-class}
\docType{class}
\alias{vgam-class}
\title{Class ``vgam'' }
\description{ Vector generalized additive models. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{vgam(...)}.
%   ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
    \item{\code{nl.chisq}:}{Object of class \code{"numeric"}.
Nonlinear chi-squared values. }
    \item{\code{nl.df}:}{Object of class \code{"numeric"}.
Nonlinear chi-squared degrees of freedom values. }
    \item{\code{spar}:}{Object of class \code{"numeric"}
containing the (scaled) smoothing parameters. }
    \item{\code{s.xargument}:}{Object of class \code{"character"}
holding the variable name of any \code{s()} terms. }
    \item{\code{var}:}{Object of class \code{"matrix"} holding
 approximate pointwise standard error information. }
    \item{\code{Bspline}:}{Object of class \code{"list"}
      holding the scaled (internal and boundary) knots, and the
      fitted B-spline coefficients. These are used
for prediction. }
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
    \item{\code{fitted.values}:}{
  Object of class
  \code{"matrix"}, from class \code{ "vlm"}.
  The fitted values. This is usually the mean but may be
  quantiles, or the location parameter, e.g., in the Cauchy model.

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
Class \code{"vglm"}, directly.
Class \code{"vlm"}, by class \code{"vglm"}.
}
\section{Methods}{
  \describe{
    \item{cdf}{\code{signature(object = "vglm")}:
cumulative distribution function.
Useful for quantile regression and extreme value data models.}
    \item{deplot}{\code{signature(object = "vglm")}:
      density plot.
      Useful for quantile regression models.}
    \item{deviance}{\code{signature(object = "vglm")}:
deviance of the model (where applicable). }
    \item{plot}{\code{signature(x = "vglm")}:
diagnostic plots. }
    \item{predict}{\code{signature(object = "vglm")}:
extract the additive predictors or
predict the additive predictors at a new data frame.}
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

Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.


%\url{http://www.stat.auckland.ac.nz/~yee}

}
\author{ Thomas W. Yee }
\note{
  VGAMs have all the slots that \code{\link{vglm}} objects
  have (\code{\link{vglm-class}}), plus the first few slots
  described in the section above.


}

%~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{
\code{\link{vgam.control}},
\code{\link{vglm}},
\code{\link[VGAM]{s}},
\code{\link{vglm-class}},
\code{\link{vglmff-class}}.


}

\examples{
# Fit a nonparametric proportional odds model
pneumo <- transform(pneumo, let = log(exposure.time))
vgam(cbind(normal, mild, severe) ~ s(let),
     cumulative(parallel = TRUE), data = pneumo)
}
\keyword{classes}
\keyword{models}
\keyword{regression}
\keyword{smooth}
