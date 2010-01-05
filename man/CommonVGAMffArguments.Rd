\name{CommonVGAMffArguments}
\alias{CommonVGAMffArguments}
\alias{TypicalVGAMfamilyFunction}
\title{Common VGAM family function Arguments }
\description{
  Here is a description of some common and typical arguments found
  in many \pkg{VGAM} family functions, e.g.,
  \code{lsigma}, 
  \code{isigma}, \code{nsimEI}, \code{parallel} and \code{zero}.

}
\usage{
TypicalVGAMfamilyFunction(lsigma = "loge", esigma = list(), isigma = NULL,
                          parallel = TRUE, shrinkage.init = 0.95,
                          nointercept = NULL, method.init = 1,
                          nsimEIM = 100, zero = NULL)
}
\arguments{

  \item{lsigma}{
  Character.
  Link function applied to a parameter and not necessarily a mean.
  See \code{\link{Links}} for a selection of choices.
  If there is only one parameter then this argument is often called
  \code{link}.

  }
  \item{esigma}{
  List.
  Extra argument allowing for additional information, specific to the
  link function.
  See \code{\link{Links}} for more information.
  If there is only one parameter then this argument is often called
  \code{earg}.

  }
  \item{isigma}{
  Optional initial values can often be inputted using an argument
  beginning with \code{"i"}.
  For example, \code{"isigma"} and \code{"ilocation"}, or just
  \code{"init"} if there is one parameter.
  A value of \code{NULL} means a value is computed internally, i.e.,
  a \emph{self-starting} \pkg{VGAM} family function.
  If a failure to converge occurs make use of these types of arguments.

  }
  \item{parallel}{
  A logical, or formula specifying which terms have equal/unequal
  coefficients.
  This argument is common in \pkg{VGAM} family functions for categorical
  responses, e.g., \code{\link{cumulative}},  \code{\link{acat}}, 
  \code{\link{cratio}}, \code{\link{sratio}}.
  For the proportional odds model (\code{\link{cumulative}}) having
  parallel constraints applied to each explanatory variable (except for
  the intercepts) means the fitted probabilities do not become negative
  or greater than 1. However this parallelism or proportional-odds
  assumption ought to be checked.

  }
  \item{nsimEIM}{
  Some \pkg{VGAM} family functions use simulation to obtain an approximate
  expected information matrix (EIM).
  For those that do, the \code{nsimEIM} argument specifies the number
  of random variates used per observation; the mean of \code{nsimEIM}
  random variates is taken.
  Thus \code{nsimEIM} controls the accuracy and a larger value may be
  necessary if the EIMs are not positive-definite.
  For intercept-only models (\code{y ~ 1)} the value of \code{nsimEIM}
  can be smaller (since the common value used is also then taken as the
  mean over the observations), especially if the number of observations
  is large.

  Some \pkg{VGAM} family functions provide two algorithms for estimating
  the EIM.
  If applicable, set \code{nsimEIM=NULL} to choose the other algorithm.

  }
  \item{method.init}{
  An integer with value \code{1} or \code{2} or \code{3} or ... which
  specifies the initialization method for some parameters or a specific
  parameter.
  If failure to converge occurs try the next higher value, and
  continue until success.
  For example, \code{method.init=1} might be the method of moments,
  and \code{method.init=2} might be another method.
  If no value of \code{method.init} works then it will be necessary
  to use arguments such as \code{isigma}.
  For many \pkg{VGAM} family functions it is advisable to try
  this argument with all possible values to safeguard against
  problems such as converging to a local solution.
  \pkg{VGAM} family functions with this argument usually correspond
  to a model or distribution that is relatively hard to fit successfully,
  therefore care is needed to ensure the global solution is obtained.
  So using all possible values that this argument supplies is a good idea.

  }
  \item{zero}{
  An integer specifying which linear/additive predictor is modelled as
  intercepts-only. That is,
  the regression coefficients are set to zero for all covariates
  except for the intercept.
  If \code{zero} is specified then it is a vector with values from the
  set \eqn{\{1,2,\ldots,M\}}.
  The value \code{zero=NULL} means model all linear/additive predictors
  as functions of the explanatory variables.
  Here, \eqn{M} is the number of linear/additive predictors.


  Some \pkg{VGAM} family functions allow the \code{zero} argument to
  accept negative values; if so then its absolute value is recycled
  over each response. For example, \code{zero=-2} would mean,
  for each response, the second linear/additive predictor is modelled
  as intercepts-only. A specific example is the \eqn{k} parameter in
  \code{\link{negbinomial}}; this \pkg{VGAM} family function can handle
  a matrix of responses.

  }
  \item{shrinkage.init}{
  Shrinkage factor \eqn{s} used for obtaining initial values.
  Numeric, between 0 and 1.
  In general, the formula used is something like
  \eqn{s \mu + (1-s) y}{s*mu + (1-s)*y}
  where \eqn{\mu}{mu} is a measure of central tendency such as a weighted
  mean or median, and \eqn{y} is the response vector.
  The the initial values are slight perturbations of the mean towards the
  actual data.
  For many type of models this method seems to work well and is often
  reasonably robust to outliers in the response.
  Often this argument is only used if
  the argument \code{method.init} is assigned a certain value.

  }
  \item{nointercept}{
  An integer-valued vector specifying which
  linear/additive predictors have no intercepts.
  Any values must be from the set \{1,2,\ldots,\eqn{M}\}.
  A value of \code{NULL} means no such constraints.

  }

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\section{Warning }{
  The \code{zero} argument is supplied for convenience but conflicts
  can arise with other arguments, e.g., the \code{constraints}
  argument of \code{\link{vglm}} and \code{\link{vgam}}.
  See Example 5 below for an example.
  If not sure, use, e.g., \code{constraints(fit)} and
  \code{coef(fit, matrix=TRUE)} to check the result of a fit \code{fit}.


  The arguments \code{zero} and \code{nointercept} can be inputted
  with values that fail. For example,
  \code{multinomial(zero = 2, nointercept = 1:3)} means the
  second linear/additive predictor is identically zero,
  which will cause a failure.


  Be careful about the use of other potentially contradictory constraints,
  e.g., \code{multinomial(zero = 2, parallel = TRUE ~ x3)}. If in doubt,
  apply \code{constraints()} to the fitted object to check.



  \pkg{VGAM} family functions with the \code{nsimEIM} may have inaccurate
  working weight matrices. If so, then the standard errors of the
  regression coefficients may be inaccurate. Thus output from
  \code{summary(fit)},
  \code{vcov(fit)},
  etc. may be misleading.

}

\details{
  Full details will be given in documentation yet to be written,
  at a later date!

}
%\references{
%}

\seealso{
  \code{\link{Links}},
  \code{\link{vglmff-class}}.

}
\author{T. W. Yee}
%\note{
%
%}
\examples{
# Example 1
cumulative()
cumulative(link="probit", reverse=TRUE, parallel=TRUE)

# Example 2
x = runif(n <- 1000)
y = rweibull(n, shape=2+exp(1+x), scale = exp(-0.5))
fit = vglm(y ~ x, weibull(lshape="logoff", eshape=list(offset=-2), zero=2))
coef(fit, mat=TRUE)
Coef(fit)

# Example 3
x = runif(n <- 500)
y1 = rnbinom(n, mu=exp(3+x), size=exp(1)) # k is size
y2 = rnbinom(n, mu=exp(2-x), size=exp(0))
fit = vglm(cbind(y1,y2) ~ x, negbinomial(zero=-2)) # multivariate response
coef(fit, matrix=TRUE)

# Example 4
\dontrun{
# fit1 and fit2 are equivalent
fit1 = vglm(ymatrix ~ x2 + x3 + x4 + x5,
            cumulative(parallel=FALSE ~ 1 + x3 + x5), mydataframe)
fit2 = vglm(ymatrix ~ x2 + x3 + x4 + x5,
            cumulative(parallel=TRUE ~ x2 + x4), mydataframe)
}

# Example 5
x = rnorm(n <- 200)
y1 = rnorm(n, mean=1-3*x, sd=exp(1+0.2*x))
y2 = rnorm(n, mean=1-3*x, sd=exp(1))
args(normal1)
fit1 = vglm(y1 ~ x, normal1)     # This is ok
fit2 = vglm(y2 ~ x, normal1(zero=2))     # This is ok

# This creates potential conflict
clist = list("(Intercept)"=diag(2), "x"=diag(2))
fit3 = vglm(y2 ~ x, normal1(zero=2), constraints=clist) # Conflict!
coef(fit3, matrix=TRUE)   # Shows that clist[["x"]] was overwritten,
constraints(fit3) # i.e., 'zero' seems to override the 'constraints' arg
}
\keyword{models}

