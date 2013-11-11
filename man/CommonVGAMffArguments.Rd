\name{CommonVGAMffArguments}
\alias{CommonVGAMffArguments}
\alias{TypicalVGAMfamilyFunction}
\title{Common VGAM family function Arguments }
\description{
  Here is a description of some common and typical arguments found
  in many \pkg{VGAM} family functions, e.g.,
  \code{lsigma}, 
  \code{isigma},
  \code{gsigma},
  \code{nsimEI},
  \code{parallel} and
  \code{zero}.

}
\usage{
TypicalVGAMfamilyFunction(lsigma = "loge",
                          isigma = NULL,
                          link.list = list("(Default)" = "identity",
                                           x2          = "loge",
                                           x3          = "logoff",
                                           x4          = "mlogit",
                                           x5          = "mlogit"),
                          earg.list = list("(Default)" = list(),
                                           x2          = list(),
                                           x3          = list(offset = -1),
                                           x4          = list(),
                                           x5          = list()),
                          gsigma = exp(-5:5),
                          parallel = TRUE,
                          shrinkage.init = 0.95,
                          nointercept = NULL, imethod = 1,
                          type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
                          probs.x = c(0.15, 0.85),
                          probs.y = c(0.25, 0.50, 0.75),
                          mv = FALSE, earg.link = FALSE,
                          whitespace = FALSE, bred = FALSE,
                          oim = FALSE, nsimEIM = 100, zero = NULL)
}
\arguments{
%                         apply.parint = FALSE,

  \item{lsigma}{
  Character.
  Link function applied to a parameter and not necessarily a mean.
  See \code{\link{Links}} for a selection of choices.
  If there is only one parameter then this argument is often called
  \code{link}.

  }
% \item{esigma}{
% List.
% Extra argument allowing for additional information, specific to the
% link function.
% See \code{\link{Links}} for more information.
% If there is only one parameter then this argument is often called
% \code{earg}.

% }


  \item{link.list, earg.list}{
  Some \pkg{VGAM} family functions
  (such as \code{\link{normal.vcm}})
  implement models with
  potentially lots of parameter link functions.
  These two arguments allow many such links and extra arguments
  to be inputted more easily.
  One has something like
  \code{link.list = list("(Default)" = "identity", x2 = "loge", x3 = "logoff")}
  and 
  \code{earg.list = list("(Default)" = list(), x2 = list(), x3 = "list(offset = -1)")}.
  Then any unnamed terms will have the default link with its
  corresponding extra argument.
  Note: the \code{\link{mlogit}} link is also possible, and if so,
  at least two instances of it are necessary.
  Then the last term is the baseline/reference group.


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
  \item{gsigma}{
  Grid-search initial values can be inputted using an argument
  beginning with \code{"g"},
  e.g., \code{"gsigma"}, \code{"gshape"} and \code{"gscale"}.
  If argument \code{isigma} is inputted then that has precedence over
  \code{gsigma}, etc.
% The actual search values will be \code{unique(sort(c(gshape)))}, etc.
  If the grid search is 2-dimensional then it is advisable not to
  make the vectors too long as a nested \code{for} loop may be used.
  Ditto for 3-dimensions.


% Then the actual search values will be \code{unique(sort(c(gshape, 1/gshape)))}, etc.


  }
  \item{parallel}{
  A logical, or a simple formula specifying which terms have equal/unequal
  coefficients.
  The formula must be simple, i.e., additive with simple main effects terms.
  Interactions and nesting etc. are not handled.
  To handle complex formulas use the \code{constraints} argument
  (of \code{\link{vglm}} etc.);
  however, there is a lot more setting up involved and things will
  not be as convenient.


  Here are some examples.
1. \code{parallel = TRUE ~ x2 + x5} means the parallelism assumption
is only applied to \eqn{X_2}, \eqn{X_5} and the intercept.
2.  \code{parallel = TRUE ~ -1}
and \code{parallel = TRUE ~ 0}
mean the parallelism assumption
is applied to \emph{no} variables at all.
Similarly,
\code{parallel = FALSE ~ -1} and
\code{parallel = FALSE ~ 0}
mean the parallelism assumption
is applied to \emph{all} the variables including the intercept.
3.  \code{parallel = FALSE ~ x2 - 1}
and \code{parallel = FALSE ~ x2 + 0}
applies the
parallelism constraint to all terms (including the intercept)
except for \eqn{X_2}.


  This argument is common in \pkg{VGAM} family functions for categorical
  responses, e.g., \code{\link{cumulative}},  \code{\link{acat}}, 
  \code{\link{cratio}}, \code{\link{sratio}}.
  For the proportional odds model (\code{\link{cumulative}}) having
  parallel constraints applied to each explanatory variable (except for
  the intercepts) means the fitted probabilities do not become negative
  or greater than 1. However this parallelism or proportional-odds
  assumption ought to be checked.


  }


% \item{apply.parint}{
% \emph{This variable will be depreciated shortly}.
% Logical.
% It refers to whether the parallelism constraint is
% applied to the intercept too.
% By default, in some models it does, in other models it does not.
% Used only if \code{parallel = TRUE} (fully or partially with
% respect to all the explanatory variables).


% }


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
  \item{type.fitted}{
  Character.
  Type of fitted value returned by the \code{fitted()} methods function.
  The first choice is always the default.
  The available choices depends on what kind of family function it is.
  Using the first few letters of the chosen choice is okay.
  See \code{\link{fittedvlm}} for more details.


  }
  \item{probs.x, probs.y}{
  Numeric, with values in (0, 1).
  The probabilites that define quantiles with respect to some vector,
  usually an \code{x} or \code{y} of some sort.
  This is used to create two subsets of data corresponding to `low' and
  `high' values of x or y.
  Each value is separately fed into the \code{probs} argument
  of \code{\link[stats:quantile]{quantile}}.
  If the data set size is small then it may be necessary to
  increase/decrease slightly the first/second values respectively.


  }
  \item{whitespace}{
  Logical.
  Should white spaces (\code{" "}) be used in the
  labelling of the linear/additive predictors?
  Setting \code{TRUE} usually results in more readability but
  it occupies more columns of the output.


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
  as intercept-only. That is, the regression coefficients are
  set to zero for all covariates except for the intercept.
  If \code{zero} is specified then it may be a vector with values
  from the set \eqn{\{1,2,\ldots,M\}}.
  The value \code{zero = NULL} means model \emph{all} linear/additive
  predictors as functions of the explanatory variables.
  Here, \eqn{M} is the number of linear/additive predictors.
  Technically, if \code{zero} contains the value \eqn{j} then
  the \eqn{j}th row of every constraint matrix (except for the
  intercept) consists of all 0 values.


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
  \item{earg.link}{
  Sometimes the link argument can receive \code{earg}-type input,
  such as \code{\link{quasibinomial}} calling \code{\link{binomial}}.
  This argument should be generally ignored.


  }
  \item{bred}{
  Logical.
  Some \pkg{VGAM} family functions will allow bias-reduction based
  on the work by Kosmidis and Firth. 
  Currently none are working yet!



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

\references{

Kosmidis, I. and Firth, D. (2009)
Bias reduction in exponential family nonlinear models.
\emph{Biometrika},
\bold{96}(4), 793--804.


%Kosmidis, I. and Firth, D. (2010)
%A generic algorithm for reducing bias in parametric estimation.
%\emph{Electronic Journal of Statistics},
%\bold{4}, 1097--1112.


}

\seealso{
  \code{\link{Links}},
  \code{\link{vglmff-class}},
  \code{\link{normal.vcm}},
  \code{\link{mlogit}}.


}
\author{T. W. Yee}

\note{
  See \code{\link{Links}} regarding a major change in
  link functions, for version 0.9-0 and higher
  (released during the 2nd half of 2012).



}

\examples{
# Example 1
cumulative()
cumulative(link = "probit", reverse = TRUE, parallel = TRUE)

# Example 2
wdata <- data.frame(x2 = runif(nn <- 1000))
wdata <- transform(wdata,
         y = rweibull(nn, shape = 2 + exp(1 + x2), scale = exp(-0.5)))
fit <- vglm(y ~ x2, weibull(lshape = logoff(offset = -2), zero = 2), wdata)
coef(fit, mat = TRUE)

# Example 3; multivariate (multiple) response
\dontrun{
ndata <- data.frame(x = runif(nn <- 500))
ndata <- transform(ndata,
           y1 = rnbinom(nn, mu = exp(3+x), size = exp(1)),  # k is size
           y2 = rnbinom(nn, mu = exp(2-x), size = exp(0)))
fit <- vglm(cbind(y1, y2) ~ x, negbinomial(zero = -2), ndata)
coef(fit, matrix = TRUE)
}
# Example 4
\dontrun{
# fit1 and fit2 are equivalent
fit1 <- vglm(ymatrix ~ x2 + x3 + x4 + x5,
             cumulative(parallel = FALSE ~ 1 + x3 + x5), mydataframe)
fit2 <- vglm(ymatrix ~ x2 + x3 + x4 + x5,
             cumulative(parallel = TRUE ~ x2 + x4), mydataframe)
}

# Example 5
gdata <- data.frame(x2 = rnorm(nn <- 200))
gdata <- transform(gdata,
           y1 = rnorm(nn, mean = 1 - 3*x2, sd = exp(1 + 0.2*x2)),
           y2 = rnorm(nn, mean = 1 - 3*x2, sd = exp(1)))
args(uninormal)
fit1 <- vglm(y1 ~ x2, uninormal, gdata)            # This is okay
fit2 <- vglm(y2 ~ x2, uninormal(zero = 2), gdata)  # This is okay

# This creates potential conflict
clist <- list("(Intercept)" = diag(2), "x2" = diag(2))
fit3 <- vglm(y2 ~ x2, uninormal(zero = 2), gdata,
             constraints = clist)  # Conflict!
coef(fit3, matrix = TRUE)  # Shows that clist[["x2"]] was overwritten,
constraints(fit3)  # i.e., 'zero' seems to override the 'constraints' arg

# Example 6 ('whitespace' argument)
pneumo <- transform(pneumo, let = log(exposure.time))
fit1 <- vglm(cbind(normal, mild, severe) ~ let,
             sratio(whitespace = FALSE, parallel = TRUE), pneumo)
fit2 <- vglm(cbind(normal, mild, severe) ~ let,
             sratio(whitespace = TRUE,  parallel = TRUE), pneumo)
head(predict(fit1), 2)  # No white spaces
head(predict(fit2), 2)  # Uses white spaces
}

\keyword{models}

