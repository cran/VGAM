\name{CommonVGAMffArguments}
\alias{CommonVGAMffArguments}
\alias{TypicalVGAMfamilyFunction}
\title{Common VGAM family function Arguments }
\description{
  Here is a description of some common and typical arguments found
  in many \pkg{VGAM} family functions, e.g.,
  \code{lsigma}, 
  \code{isigma},
  \code{nsimEI},
  \code{parallel} and
  \code{zero}.

}
\usage{
TypicalVGAMfamilyFunction(lsigma = "loge", esigma = list(), isigma = NULL,
                          parallel = TRUE, shrinkage.init = 0.95,
                          nointercept = NULL, imethod = 1,
                          prob.x = c(0.15, 0.85), mv = FALSE,
                          oim = FALSE, nsimEIM = 100, zero = NULL)
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
  If applicable, set \code{nsimEIM = NULL} to choose the other algorithm.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} or \code{3} or ... which
  specifies the initialization method for some parameters or a specific
  parameter.
  If failure to converge occurs try the next higher value, and
  continue until success.
  For example, \code{imethod = 1} might be the method of moments,
  and \code{imethod = 2} might be another method.
  If no value of \code{imethod} works then it will be necessary
  to use arguments such as \code{isigma}.
  For many \pkg{VGAM} family functions it is advisable to try
  this argument with all possible values to safeguard against
  problems such as converging to a local solution.
  \pkg{VGAM} family functions with this argument usually correspond
  to a model or distribution that is relatively hard to fit successfully,
  therefore care is needed to ensure the global solution is obtained.
  So using all possible values that this argument supplies is a good idea.


  }
  \item{prob.x}{
  Numeric, of length two.
  The probabilites that define quantiles with respect to some vector,
  usually an \code{x} of some sort.
  This is used to create two subsets of data corresponding to `low' and
  `high' values of x.
  Each value is separately fed into the \code{probs} argument
  of \code{\link[stats:quantile]{quantile}}.
  If the data set size is small then it may be necessary to
  increase/decrease slightly the first/second values respectively.


  }
  \item{oim}{
  Logical.
  Should the observed information matrices (OIMs) be used for
  the working weights?
  In general, setting \code{oim = TRUE} means the Newton-Raphson
  algorithm, and \code{oim = FALSE} means Fisher-scoring.
  The latter uses the EIM, and is usually recommended.
  If \code{oim = TRUE} then \code{nsimEIM} is ignored.

  }
  \item{zero}{
  An integer specifying which linear/additive predictor is modelled
  as intercepts-only. That is, the regression coefficients are
  set to zero for all covariates except for the intercept.
  If \code{zero} is specified then it is a vector with values
  from the set \eqn{\{1,2,\ldots,M\}}.
  The value \code{zero = NULL} means model all linear/additive
  predictors as functions of the explanatory variables.
  Here, \eqn{M} is the number of linear/additive predictors.


  Some \pkg{VGAM} family functions allow the \code{zero} argument to
  accept negative values; if so then its absolute value is recycled
  over each (usual) response. For example, \code{zero = -2} for the
  two-parameter negative binomial distribution would mean,
  for each response, the second linear/additive predictor is modelled
  as intercepts-only. That is, for all the \eqn{k} parameters in
  \code{\link{negbinomial}} (this \pkg{VGAM} family function can handle
  a matrix of responses).


  Suppose \code{zero = zerovec} where \code{zerovec} is a vector
  of negative values. If \eqn{G} is the usual \eqn{M} value for
  a univariate response then the actual values for argument \code{zero}
  are all values in
  \code{c(abs(zerovec), G + abs(zerovec), 2*G + abs(zerovec), ... )}
  lying in the integer range \eqn{1} to \eqn{M}.
  For example, setting \code{zero = -c(2, 3)} for a matrix response
  of 4 columns with \code{\link{zinegbinomial}} (which
  usually has \eqn{G = M = 3} for a univariate response)
  would be equivalent to \code{zero = c(2, 3, 5, 6, 8, 9, 11, 12)}.
  This example has \eqn{M = 12}.
  Note that if \code{zerovec} contains negative values then their
  absolute values should be elements from the set \code{1:G}.


  Note: \code{zero} may have positive and negative values,
  for example, setting \code{zero = c(-2, 3)} in the above example
  would be equivalent to \code{zero = c(2, 3, 5, 8, 11)}.


  }
  \item{shrinkage.init}{
  Shrinkage factor \eqn{s} used for obtaining initial values.
  Numeric, between 0 and 1.
  In general, the formula used is something like
  \eqn{s \mu + (1-s) y}{s*mu + (1-s)*y}
  where \eqn{\mu}{mu} is a measure of central tendency such as a weighted
  mean or median, and \eqn{y} is the response vector.
  For example, the initial values are slight perturbations of
  the mean towards the actual data.
  For many types of models this method seems to work well and is often
  reasonably robust to outliers in the response.
  Often this argument is only used if
  the argument \code{imethod} is assigned a certain value.

  }
  \item{nointercept}{
  An integer-valued vector specifying which
  linear/additive predictors have no intercepts.
  Any values must be from the set \{1,2,\ldots,\eqn{M}\}.
  A value of \code{NULL} means no such constraints.

  }

  \item{mv}{
  Logical.
  Some \pkg{VGAM} family functions allow a multivariate or vector response.
  If so, then usually the response is a matrix with columns
  corresponding to the individual response variables.
  They are all fitted simultaneously.
  Arguments such as \code{parallel} may then be useful to allow
  for relationships between the regressions of each response variable.
  If \code{mv = TRUE} then sometimes the response is interpreted
  differently, e.g., \code{\link{posbinomial}} chooses the first
  column of a matrix response as success and combines the other
  columns as failure, but when \code{mv = TRUE} then each column
  of the response matrix is the number of successes and the
  \code{weights} argument is of the same dimension as the
  response and contains the number of trials.

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
  \code{coef(fit, matrix = TRUE)} to check the result of a fit \code{fit}.


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
cumulative(link = "probit", reverse = TRUE, parallel = TRUE)

# Example 2
wdata <- data.frame(x = runif(nn <- 1000))
wdata <- transform(wdata,
           y = rweibull(nn, shape = 2 + exp(1+x), scale = exp(-0.5)))
fit = vglm(y ~ x, weibull(lshape = "logoff", eshape = list(offset = -2),
                          zero = 2), wdata)
coef(fit, mat = TRUE)

# Example 3; multivariate (multiple) response
ndata <- data.frame(x = runif(nn <- 500))
ndata <- transform(ndata,
           y1 = rnbinom(nn, mu = exp(3+x), size = exp(1)), # k is size
           y2 = rnbinom(nn, mu = exp(2-x), size = exp(0)))
fit <- vglm(cbind(y1, y2) ~ x, negbinomial(zero = -2), ndata)
coef(fit, matrix = TRUE)

# Example 4
\dontrun{
# fit1 and fit2 are equivalent
fit1 <- vglm(ymatrix ~ x2 + x3 + x4 + x5,
             cumulative(parallel = FALSE ~ 1 + x3 + x5), mydataframe)
fit2 <- vglm(ymatrix ~ x2 + x3 + x4 + x5,
             cumulative(parallel = TRUE ~ x2 + x4), mydataframe)
}

# Example 5
gdata <- data.frame(x = rnorm(nn <- 200))
gdata <- transform(gdata,
           y1 = rnorm(nn, mean = 1 - 3*x, sd = exp(1 + 0.2*x)),
           y2 = rnorm(nn, mean = 1 - 3*x, sd = exp(1)))
args(normal1)
fit1 <- vglm(y1 ~ x, normal1, gdata) # This is ok
fit2 <- vglm(y2 ~ x, normal1(zero = 2), gdata) # This is ok

# This creates potential conflict
clist <- list("(Intercept)" = diag(2), "x" = diag(2))
fit3 <- vglm(y2 ~ x, normal1(zero = 2), gdata,
             constraints = clist) # Conflict!
coef(fit3, matrix = TRUE)   # Shows that clist[["x"]] was overwritten,
constraints(fit3) # i.e., 'zero' seems to override the 'constraints' arg
}
\keyword{models}

