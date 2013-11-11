% 20120312
% Modified from file src/library/stats/man/influence.measures.Rd

\name{hatvalues}
%\title{Regression Deletion Diagnostics}
\title{Hat Values and Regression Deletion Diagnostics}

%\concept{studentized residuals}
%\concept{standardized residuals}
%\concept{Cook's distances}
%\concept{Covariance ratios}
\concept{DFBETAs}
%\concept{DFFITs}

%\alias{influence.measures}
%\alias{print.infl}
%\alias{summary.infl}
%\alias{hat}
\alias{hatvalues}
%\alias{hatvalues.lm}
\alias{hatvaluesvlm}
\alias{hatplot}
\alias{hatplot.vlm}
%\alias{rstandard}
%\alias{rstandard.lm}
%\alias{rstandard.glm}
%\alias{rstudent}
%\alias{rstudent.lm}
%\alias{rstudent.glm}
\alias{dfbeta}
\alias{dfbetavlm}
%\alias{dfbetas}
%\alias{dfbetas.lm}
%\alias{dffits}
%\alias{covratio}
%\alias{cooks.distance}
%\alias{cooks.distance.lm}
%\alias{cooks.distance.glm}
\usage{
hatvalues(model, \dots)
hatvaluesvlm(model, type = c("diagonal", "matrix", "centralBlocks"), \dots)
hatplot(model, \dots)
hatplot.vlm(model, multiplier = c(2, 3), lty = "dashed",
            xlab = "Observation", ylab = "Hat values", ylim = NULL, \dots)
dfbetavlm(model, maxit.new = 1,
          trace.new = FALSE,
          smallno = 1.0e-8, ...)
}
\arguments{
  \item{model}{an \R object, typically returned by \code{\link{vglm}}.
%or \code{\link{glm}}.


  }
  \item{type}{Character.
  The default is the first choice, which is
  a \eqn{nM \times nM}{nM x nM} matrix.
  If \code{type = "matrix"} then the \emph{entire} hat matrix is
  returned.
  If \code{type = "centralBlocks"} then \eqn{n} central
  \eqn{M \times M}{M x M} block matrices, in matrix-band format.

 }

% \item{diag}{Logical. If \code{TRUE} then the diagonal elements
% of the hat matrix are returned, else the \emph{entire} hat matrix is
% returned.
% In the latter case, it is a \eqn{nM \times nM}{nM x nM} matrix.


% }


  \item{multiplier}{Numeric, the multiplier.
  The usual rule-of-thumb is that values greater than two or three
  times the average leverage (at least for the linear model) should
  be checked.


  }
  \item{lty, xlab, ylab, ylim}{Graphical parameters, see
  \code{\link[graphics]{par}} etc.
  The default of \code{ylim} is \code{c(0, max(hatvalues(model)))}
  which means that if the horizontal dashed lines cannot be seen
  then there are no particularly influential observations.


  }
  \item{maxit.new, trace.new, smallno}{
  Having \code{maxit.new = 1} will give a one IRLS step approximation
  from the ordinary solution (and no warnings!).
  Else having \code{maxit.new = 10}, say, should usually mean
  convergence will occur for all observations when they are
  removed one-at-a-time.
  Else having \code{maxit.new = 2}, say, should usually mean
  some lack of convergence will occur when observations are
  removed one-at-a-time.
  Setting \code{trace.new = TRUE} will produce some running output
  at each IRLS iteration and for each individual row of the model matrix.
  The argument \code{smallno} multiplies each value of the
  original prior weight (often unity); setting it identically
  to zero will result in an error, but setting a very small value
  effectively removes that observation.


  }

% \item{infl}{influence structure as returned by
%   \code{\link{lm.influence}} or \code{\link{influence}} (the latter
%   only for the \code{glm} method of \code{rstudent} and
%   \code{cooks.distance}).}
% \item{res}{(possibly weighted) residuals, with proper default.}
% \item{sd}{standard deviation to use, see default.}
% \item{dispersion}{dispersion (for \code{\link{glm}} objects) to use,
%   see default.}
% \item{hat}{hat values \eqn{H_{ii}}{H[i,i]}, see default.}
% \item{type}{type of residuals for \code{glm} method for \code{rstandard.}}

% \item{x}{the \eqn{X} or design matrix.}
% \item{intercept}{should an intercept column be prepended to \code{x}?}
  \item{\dots}{further arguments,
  for example, graphical parameters for \code{hatplot.vlm()}.
% passed to or from other methods.


  }

}
\description{
  When complete, a
  suite of functions that can be used to compute some of the
  regression (leave-one-out deletion) diagnostics,
  for the VGLM class.


% This suite of functions can be used to compute some of the
% regression (leave-one-out deletion) diagnostics for linear and
% generalized linear models discussed in Belsley, Kuh and Welsch
% (1980), Cook and Weisberg (1982), etc.


}
\details{
  The invocation \code{hatvalues(vglmObject)} should return a
  \eqn{n \times M}{n x M} matrix of the diagonal elements of the
  hat (projection) matrix of a \code{\link{vglm}} object.
  To do this,
  the QR decomposition of the object is retrieved or
  reconstructed, and then straightforward calculations
  are performed.


  The invocation \code{hatplot(vglmObject)} should plot
  the diagonal of the hat matrix for each of the \eqn{M}
  linear/additive predictors.
  By default, two horizontal dashed lines are added;
  hat values higher than these ought to be checked.




% The primary high-level function is \code{influence.measures}
% which produces a class \code{"infl"} object tabular display
% showing the DFBETAS for each model variable, DFFITS, covariance
% ratios, Cook's distances and the diagonal elements of the
% hat matrix.  Cases which are influential with respect to any
% of these measures are marked with an asterisk.


% The functions \code{dfbetas}, \code{dffits}, \code{covratio}
% and \code{cooks.distance} provide direct access to the
% corresponding diagnostic quantities.  Functions \code{rstandard}
% and \code{rstudent} give the standardized and Studentized
% residuals respectively. (These re-normalize the residuals to
% have unit variance, using an overall and leave-one-out measure
% of the error variance respectively.)


% Values for generalized linear models are approximations, as
% described in Williams (1987) (except that Cook's distances
% are scaled as \eqn{F} rather than as chi-square values).  The
% approximations can be poor when some cases have large influence.


% The optional \code{infl}, \code{res} and \code{sd} arguments are
% there to encourage the use of these direct access functions,
% in situations where, e.g., the underlying basic influence
% measures (from \code{\link{lm.influence}} or the generic
% \code{\link{influence}}) are already available.


% Note that cases with \code{weights == 0} are \emph{dropped} from all
% these functions, but that if a linear model has been fitted with
% \code{na.action = na.exclude}, suitable values are filled in for the
% cases excluded during fitting.


% The function \code{hat()} exists mainly for S (version 2)
% compatibility; we recommend using \code{hatvalues()} instead.



}
\note{
  It is hoped, soon, that the full suite of functions described at
  \code{\link[stats]{influence.measures}} will be written for VGLMs.
  This will enable general regression deletion diagnostics to be
  available for the entire VGLM class.


%  For \code{hatvalues}, \code{dfbeta}, and \code{dfbetas}, the method
%  for linear models also works for generalized linear models.


}
\author{
  T. W. Yee.


}
%\references{
%  Belsley, D. A., Kuh, E. and Welsch, R. E. (1980)
%  \emph{Regression Diagnostics}.
%  New York: Wiley.
%
%  Cook, R. D. and Weisberg, S. (1982)
%  \emph{Residuals and Influence in Regression}.
%  London: Chapman and Hall.
%
%  Williams, D. A. (1987)
%  Generalized linear model diagnostics using the deviance and single
%  case deletions. \emph{Applied Statistics} \bold{36}, 181--191.
%
%  Fox, J. (1997)
%  \emph{Applied Regression, Linear Models, and Related Methods}. Sage.
%
%  Fox, J. (2002)
%  \emph{An R and S-Plus Companion to Applied Regression}.
%  Sage Publ.; \url{http://www.socsci.mcmaster.ca/jfox/Books/Companion/}.
%
%
%}

\seealso{
  \code{\link{vglm}},
  \code{\link{cumulative}},
  \code{\link[stats]{influence.measures}}.


}
\examples{
# Proportional odds model, p.179, in McCullagh and Nelder (1989)
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let, cumulative, data = pneumo)
hatvalues(fit)  # n x M matrix, with positive values
all.equal(sum(hatvalues(fit)), fit@rank)  # Should be TRUE
\dontrun{ par(mfrow = c(1, 2))
hatplot(fit, ylim = c(0, 1), las = 1, col = "blue") }
}
\keyword{regression}



