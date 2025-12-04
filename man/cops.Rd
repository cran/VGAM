\name{cops}
\alias{cops}
\alias{copsvglm}
\alias{cops,vglm-method}
%- Also NEED an '\alias' for EACH other topic
%- documented here.
\title{ Centre of the Parameter Space }
\description{
  Returns a vector similar to \code{coef(object)},
  comprising the centre of the parameter space (COPS)
  values,
  given a fitted VGLM.


}
\usage{
cops(object, ...)
copsvglm(object, level = 0.999,
    beta.range = confint(object, level = level),
    muxrange = 4, iter.max = 8, tol = 1e-5, subset = NULL,
    do1 = TRUE, slowtrain = FALSE, ...)
}
%- maybe also 'usage' for other objects
%- documented here.
%muxrange = 4, iter.max = 20,
%    tol = 1e-5, nmax = NULL, subset = NULL, ...
\arguments{
  \item{object}{
    A \code{\link{vglm}} object,
    e.g., representing a logistic regression.



% However, this function will not work for all
% such objects.



  }
  \item{level}{
    Fed into \code{\link[stats]{confint}} to obtain
    an interval for each coefficient.
    Then the centered confidence limits are multiplied
    by \code{muxrange} and \code{iter} to widen
    the interval, so
    as to allow a gridsearch.
    Hopefully each maximum is bracketed by
    two endpoints.
    So a value of \code{level} closer to 1 will mean
    the initial confidence intervals will be wider.


  }

  \item{beta.range}{
    Numeric.
    Interval for the numerical search.
    It is common for the COPS to be relatively
    near the estimated regression coefficients.
    Argument \code{level} provides an initial
    wide range about these estimates.
    Alternatively, a 2-column matrix can be inputted
    to replace the default, where the first
    column is the lower limit, etc.


%   , however they should not
%   be excessively far because of numerical
%   problems.
%   If so, a warning will be given.


  }
  \item{muxrange, iter.max}{
    Numeric.
    Used to obtain the interval endpoints that
    cover each COPS value.
    The initial confidence intervals are linearly
    widened by multiplication by \code{muxrange}
    and \code{1:iter.max}.
    This provides the positions for a gridsearch
    whereby the maximum can be covered.
    If the highest information value of the is at
    the most extreme gridpoints (i.e., the first
    or the last) then a warning is given for
    the user to increase
    at least one of these arguments.
    Note that \code{muxrange > 1} in order
    to expand the gridsearch outward.


  }
  \item{tol}{
    Numeric.
    Fed into \code{tol} in \code{\link[stats]{optimize}}.
    Probably \code{tol} should only occasionally used,
    e.g.,
    if \code{coef(object)} contains some very large
    or small values.


  }
  \item{do1}{
  Logical. Include the intercepts?
  This operates in conjunction with \code{subset}.



  }
  \item{slowtrain}{
  Logical. Used internally, not important. 



  }
  \item{subset}{
  Same as in \code{\link{wsdm}}.
  To control whether intercepts are included,
  use \code{do1}.


% Specify a subset of the regression
% coefficients?
% May be numeric, logical or character.
% Should work if \code{coef(object)[subset]}
% works.
% The default means \code{subset <- TRUE}.


  }



  \item{\dots}{
    currently unused but may be used in the
    future for further arguments passed into
    other methods.
% e.g., \code{subset}.


  }
}
\details{
  For many models, some COPS values will be
  \code{Inf} or \code{-Inf}
  so that manual checking is needed,
  for example, \code{\link{poissonff}}.
  The answers returned by this function only
  make sense if the COPSs are in the
  interior of the parameter space.
  This function was written specifically for
  logistic regression but has much wider
  applicability.
  Currently the result returned depends critically
  on \code{beta.range},
  \code{muxrange} and \code{iter.max}
  so that the answer ought to be checked.


% after several values are fed into that argument.


%  Each value returned may be effectively
%  that of \code{beta.range}
%  or \code{NA}.



  Yee (2025) shows that not all VGLMs have a
  finite COPS. For example, an intercept-only
  \code{\link{propodds}} with second
  cumulative probability less than 0.25 will
  not have a finite COPS for the first
  intercept.
  For now, this function should probably
  only applied to
  \code{\link{binomialff}} models
  with the four most common link functions.



% Using \code{subset} (with one selection)
% with \code{beta.range} gives greater
% control for computing one COPS value.
% In contrast, computing them all with
% a single \code{beta.range} can easily
% lead to problems.




}
\value{
  A named vector, similar to \code{\link{coefvlm}}.


%  If \code{trace.} then a list is returned,
%  having a component comprising a
%  matrix of function evaluations used by
%  \code{\link[stats]{optimize}}.


}
\references{



%Yee, T. W. (2021).
%Some new results concerning the Hauck-Donner
%effect.
%\emph{Manuscript in preparation}.



%Yee, T. W. (2022).
%Some new results concerning the Wald tests and
%the parameter space.
%\emph{In preparation}.



%Yee, T. W. (2024).
%Musings and new results on the parameter space.
%\emph{Under review}.




Yee, T. W. (2025).
The centre of the parameter space and its
mapping with the
WSDM function:
A diagnostic for logistic regression and
beyond.
\emph{In preparation}.



}
\author{ Thomas W. Yee.  }

%\section{Warning }{



%}

\note{


  This function could be made more robust
  in the future.
  Of course, the function will not return finite
  results if the COPS is not finite, e.g.,
  a Poisson regression.




% may be fragile
% experimental and hopefully will
% be generalized to other regression models


}
\seealso{
  \code{\link{hdeff}},
  \code{\link{wsdm}},
  \code{\link{cops3.glm}},
  \code{\link{binomialff}},
  \code{\link{coefvlm}},
  \code{\link[stats]{coef}}.



}

\examples{
\dontrun{data("xs.nz", package = "VGAMdata")
data1 <- na.omit(xs.nz[, c("age", "cancer", "sex")])
fit1 <- vglm(cancer ~ age + sex, binomialff, data1)
cops(fit1)
}}

% Add >=1 standard keywords, see file 'KEYWORDS'
% in the R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{htest}
%\concept{Hauck--Donner effect}




