% File src/library/stats/man/step.Rd
% Part of the R package, https://www.R-project.org
% Copyright 1995-2014 R Core Team
% Distributed under GPL 2 or later

\name{step4}
\alias{step4}
\alias{step4vglm}
\title{
Choose a model by AIC in a Stepwise Algorithm

}
\description{
  Select a formula-based model by AIC.

}
\usage{
step4(object, \dots)
step4vglm(object, scope, direction = c("both", "backward", "forward"),
          trace = 1, keep = NULL, steps = 1000, k = 2, \dots)
}
%#           constraints = NULL,
\arguments{
  \item{object}{
    an object of class \code{"vglm"}.
    This is used as the initial model in the stepwise search.


%    It is strongly recommended that this be the full model
%    because a backward direction is taken first.


  }
  \item{scope}{
  See \code{\link[stats]{step}}.


%    defines the range of models examined in the stepwise search.
%    This should be either a single formula, or a list containing
%    components \code{upper} and \code{lower}, both formulae.  See the
%    details for how to specify the formulae and how they are used.



  }

  \item{direction}{
  See \code{\link[stats]{step}}.


%   the mode of stepwise search, can be one of \code{"both"},
%   \code{"backward"}, or \code{"forward"},
%   with a default being the first value.
%   If the \code{scope} argument is missing the default for
%   \code{direction} is also \code{"backward"}.
%   Values can be abbreviated.


  }


%  \item{trace}{
%    if positive, information is printed during the running
%    of \code{step}.
%    Larger values may give more detailed information.



%  }
  \item{trace, keep}{
  See \code{\link[stats]{step}}.


%   a filter function whose input is a fitted model object and the
%    associated \code{AIC} statistic, and whose output is arbitrary.
%    Typically \code{keep} will select a subset of the components of
%    the object and return them. The default is not to keep anything.



  }
  \item{steps, k}{
  See \code{\link[stats]{step}}.


%    the maximum number of steps to be considered.  The default is 1000
%    (essentially as many as required).  It is typically used to stop the
%    process early.




  }
%  \item{k}{
%  See \code{\link[stats]{step}}.


%    the multiple of the number of degrees of freedom used for the
%    penalty.  Only \code{k = 2} gives the genuine AIC: \code{k = log(n)}
%    is sometimes referred to as BIC or SBC.



%  }
  \item{\dots}{
    any additional arguments to
    \code{\link{extractAIC.vglm}},
    \code{\link{drop1.vglm}} and
    \code{\link{add1.vglm}}.


  }
}
\value{
  The results are placed in the \code{post} slot of the
  stepwise-selected model that is returned.
  There are up to two additional components.
  There is an \code{"anova"}
  component corresponding to the steps taken in the search,
  as well as a \code{"keep"} component if the \code{keep=} argument
  was supplied in the call.



% the stepwise-selected model is returned, with up to two additional
% components.  There is an \code{"anova"} component corresponding to the
% steps taken in the search, as well as a \code{"keep"} component if the
% \code{keep=} argument was supplied in the call. The
% \code{"Resid. Dev"} column of the analysis of deviance table refers to
% a constant minus twice the maximized log likelihood: it will be a
% deviance only in cases where a saturated model is well-defined (thus
% excluding \code{lm}, \code{aov} and \code{survreg} fits, for example).




}
\details{
  This function is a direct adaptation of
  \code{\link[stats]{step}}
  for \code{\link{vglm-class}} objects.
  Since \code{\link[stats]{step}} is not generic,
  the name \code{step4()} was adopted
  and it \emph{is} generic, as well as being S4 rather than S3.
  It is the intent that this function should work as similar as
  possible to \code{\link[stats]{step}}.



  Internally, the methods function for \code{\link{vglm-class}}
  objects calls \code{\link{add1.vglm}} and
  \code{\link{drop1.vglm}}
  repeatedly.




% ; it will work for any method for which they work, and that
% is determined by having a valid method for \code{\link{extractAIC}}.
% When the additive constant can be chosen so that AIC is equal to
% Mallows' \eqn{C_p}{Cp}, this is done and the tables are labelled
% appropriately.



% The set of models searched is determined by the \code{scope} argument.
% The right-hand-side of its \code{lower} component is always included
% in the model, and right-hand-side of the model is included in the
% \code{upper} component.  If \code{scope} is a single formula, it
% specifies the \code{upper} component, and the \code{lower} model is
% empty.  If \code{scope} is missing, the initial model is used as the
% \code{upper} model.



% Models specified by \code{scope} can be templates to update
% \code{object} as used by \code{\link{update.formula}}.  So using
% \code{.} in a \code{scope} formula means \sQuote{what is
% already there}, with \code{.^2} indicating all interactions of
% existing terms.



% There is a potential problem in using \code{\link{glm}} fits with a
% variable \code{scale}, as in that case the deviance is not simply
% related to the maximized log-likelihood.  The \code{"glm"} method for
% function \code{\link{extractAIC}} makes the
% appropriate adjustment for a \code{gaussian} family, but may need to be
% amended for other cases.  (The \code{binomial} and \code{poisson}
% families have fixed \code{scale} by default and do not correspond
% to a particular maximum-likelihood problem for variable \code{scale}.)



}
%\note{
%  This is a minimal implementation.  Use \code{\link[MASS]{stepAIC}}
%  in package \CRANpkg{MASS} for a wider range of object classes.



%}
\section{Warning}{
  In general,
  the same warnings in
  \code{\link[stats]{drop1.glm}} and
  \code{\link{drop1.vglm}}
  apply here.


  
  This function

  
  

%  The model fitting must apply the models to the same dataset. This
%  may be a problem if there are missing values and \R's default of
%  \code{na.action = na.omit} is used.  We suggest you remove the
%  missing values first.



%  Calls to the function \code{\link{nobs}} are used to check that the
%  number of observations involved in the fitting process remains
%  unchanged.



}
\seealso{
  \code{\link{add1.vglm}},
  \code{\link{drop1.vglm}},
  \code{\link{vglm}},
  \code{\link{trim.constraints}},
  \code{\link[stats]{add1.glm}},
  \code{\link[stats]{drop1.glm}},
  \code{\link{backPain2}},
  \code{\link[stats]{step}},
  \code{\link[stats]{update}}.



  
%  \code{\link[MASS]{stepAIC}} in \CRANpkg{MASS},



}
%\references{
%  Hastie, T. J. and Pregibon, D. (1992).
%  \emph{Generalized linear models.}
%  Chapter 6 of \emph{Statistical Models in S}
%  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

%  Venables, W. N. and Ripley, B. D. (2002).
%  \emph{Modern Applied Statistics with S.}
%  New York: Springer (4th ed).
%}



%\author{
%  B. D. Ripley: \code{step} is a slightly simplified version of
%  \code{\link[MASS]{stepAIC}} in package \CRANpkg{MASS}
%  (Venables & Ripley, 2002 and earlier editions).

%  The idea of a \code{step} function follows that described in
%  Hastie & Pregibon (1992); but the implementation
%  in \R is more general.
%}




\examples{
data("backPain2", package = "VGAM")
summary(backPain2)
fit1 <- vglm(pain ~ x2 + x3 + x4 + x2:x3 + x2:x4 + x3:x4,
             propodds, data = backPain2)
spom1 <- step4(fit1)
summary(spom1)
spom1@post$anova
}
\keyword{models}


%\donttest{}
%\dontshow{utils::example("lm", echo = FALSE)}
