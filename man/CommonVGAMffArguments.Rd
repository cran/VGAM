\name{CommonVGAMffArguments}
\alias{CommonVGAMffArguments}
\alias{TypicalVGAMfamilyFunction}
\title{Common VGAM Family Function Arguments }
\description{
  Here is a description of some common and
  typical arguments found in many \pkg{VGAM}
  family functions, e.g.,
  \code{zero},
  \code{lsigma},
  \code{isigma},
  \code{gsigma},
  \code{eq.mean},
  \code{nsimEI} and
  \code{parallel}.

}
\usage{
TypicalVGAMfamilyFunction(lsigma = "loglink",
                          isigma = NULL,
                          zero = NULL, gsigma = exp(-5:5),
                          eq.mean = FALSE,
                          parallel = TRUE,
                          imethod = 1,
                          vfl = FALSE, Form2 = NULL, 
                          type.fitted = c("mean", "quantiles", "Qlink",
                                          "pobs0", "pstr0", "onempstr0"),
                          percentiles = c(25, 50, 75),
                          probs.x = c(0.15, 0.85),
                          probs.y = c(0.25, 0.50, 0.75),
                          multiple.responses = FALSE, earg.link = FALSE,
                          ishrinkage = 0.95, nointercept = NULL,
                          whitespace = FALSE, bred = FALSE, lss = TRUE,
                          oim = FALSE, nsimEIM = 100, byrow.arg = FALSE,
                          link.list = list("(Default)" = "identitylink",
                                           x2          = "loglink",
                                           x3          = "logofflink",
                                           x4          = "multilogitlink",
                                           x5          = "multilogitlink"),
                          earg.list = list("(Default)" = list(),
                                           x2          = list(),
                                           x3          = list(offset = -1),
                                           x4          = list(),
                                           x5          = list()),
                          Thresh = NULL, nrfs = 1)
}
\arguments{
%                         apply.parint = FALSE,

  \item{lsigma}{
  Character.
  Link function applied to a parameter and not
  necessarily a mean.  See \code{\link{Links}}
  for a selection of choices.  If there is only
  one parameter then this argument is often
  called \code{link}.

  }
% \item{esigma}{
% List.
% Extra argument allowing for additional information, specific to the
% link function.
% See \code{\link{Links}} for more information.
% If there is only one parameter then this argument is often called
% \code{earg}.

% }




  \item{isigma}{
  Optional initial values can often be
  inputted using an argument beginning with
  \code{"i"}.  For example, \code{"isigma"} and
  \code{"ilocation"}, or just \code{"init"}
  if there is one parameter.  A value of
  \code{NULL} means a value is computed
  internally, i.e., a \emph{self-starting}
  \pkg{VGAM} family function.  If a failure
  to converge occurs make use of these types
  of arguments.


  }
  \item{zero}{
    An important argument,
  either an integer vector, or a vector of
  character strings.


  If an integer, then it specifies which
  linear/additive predictor is modelled
  as \emph{intercept-only}.  That is,
  the regression coefficients are set to
  zero for all covariates except for the
  intercept.  If \code{zero} is specified
  then it may be a vector with values from
  the set \eqn{\{1,2,\ldots,M\}}.  The value
  \code{zero = NULL} means model \emph{all}
  linear/additive predictors as functions of
  the explanatory variables.  Here, \eqn{M}
  is the number of linear/additive predictors.
  Technically, if \code{zero} contains the
  value \eqn{j} then the \eqn{j}th row of every
  constraint matrix (except for the intercept)
  consists of all 0 values.


  
  Some \pkg{VGAM} family functions allow the
  \code{zero} argument to accept negative
  values; if so then its absolute value is
  recycled over each (usual) response. For
  example, \code{zero = -2} for the
  two-parameter negative binomial distribution
  would mean, for each response, the second
  linear/additive predictor is modelled as
  intercepts-only. That is, for all the \eqn{k}
  parameters in \code{\link{negbinomial}}
  (this \pkg{VGAM} family function can handle
  a matrix of responses).



  
  Suppose \code{zero = zerovec} where
  \code{zerovec} is a vector of negative
  values. If \eqn{G} is the usual \eqn{M}
  value for a univariate response then the
  actual values for argument \code{zero}
  are all values in \code{c(abs(zerovec), G +
  abs(zerovec), 2*G + abs(zerovec), ... )} lying
  in the integer range \eqn{1} to \eqn{M}.
  For example, setting \code{zero = -c(2,
  3)} for a matrix response of 4 columns with
  \code{\link{zinegbinomial}} (which usually
  has \eqn{G = M = 3} for a univariate response)
  would be equivalent to
  \code{zero = c(2, 3, 5, 6, 8, 9, 11, 12)}.
  This example has \eqn{M = 12}.
  Note that if \code{zerovec} contains
  negative values then their absolute values
  should be elements from the set \code{1:G}.


  
  Note: \code{zero} may have positive and
  negative values, for example, setting
  \code{zero = c(-2, 3)} in the above example
  would be equivalent to
  \code{zero = c(2, 3, 5, 8, 11)}.



  The argument \code{zero} also accepts
  a character vector (for \pkg{VGAM}
  1.0-1 onwards).  Each value is fed into
  \code{\link[base]{grep}} with \code{fixed =
  TRUE}, meaning that wildcards \code{"*"}
  are not useful.  See the example below---all
  the variants work; those with \code{LOCAT}
  issue a warning that that value is unmatched.
Importantly, the parameter names are
\code{c("location1", "scale1", "location2",
"scale2")} because there are 2 responses.  Yee
(2015) described \code{zero} for only numerical
input.  Allowing character input is particularly
important when the number of parameters cannot
be determined without having the actual data
first. For example, with time series data,
an ARMA(\eqn{p},\eqn{q}) process might have
parameters \eqn{\theta_1,\ldots,\theta_p} which
should be intercept-only by default. Then
specifying a numerical default value for
\code{zero} would be too difficult (there
are the drift and scale parameters too).
However, it is possible with the character
representation: \code{zero = "theta"} would
achieve this.  In the future, most \pkg{VGAM}
family functions might be converted
  to the character representation---the
  advantage being that it is more readable.
  When programming a \pkg{VGAM} family function
  that allows character input, the variable
  \code{predictors.names} must be assigned
  correctly.

  

%Note that \code{zero} accepts wildcards (cf. the Linux operating system):
%\code{"location*"} means that \emph{all} location parameters
%are intercept-only.
% When programming a \pkg{VGAM} family function that allows character
% input, the variables \code{parameters.names}
% and \code{Q1}




% 20230209; the following paragraph comes from Comm*Rd:


  If the \code{constraints} argument is used
  then the family function's \code{zero}
  argument (if it exists) needs to be set to
  \code{NULL}. This avoids what could be a
  probable contradiction.  Sometimes setting
  other arguments related to constraint matrices
  to \code{FALSE} is also a good idea, e.g.,
  \code{parallel = FALSE},
  \code{exchangeable = FALSE}.



  }
  \item{gsigma}{
  Grid-search initial values can be inputted
  using an argument beginning with \code{"g"},
  e.g., \code{"gsigma"}, \code{"gshape"} and
  \code{"gscale"}.  If argument \code{isigma}
  is inputted then that has precedence over
  \code{gsigma}, etc.
  If the grid search is 2-dimensional then it
  is advisable not to make the vectors too long
  as a nested \code{for} loop may be used.
  Ditto for 3-dimensions etc.  Sometimes a
  \code{".mux"} is added as a suffix, e.g.,
  \code{gshape.mux}; this means that the grid
  is created relatively and not absolutely,
  e.g., its values are multipled by some single
  initial estimate of the parameter in order
  to create the grid on an absolute scale.


% The actual search values will be
% \code{unique(sort(c(gshape)))}, etc.



  Some family functions have an argument
  called \code{gprobs.y}.  This is fed
  into the \code{probs} argument of
  \code{\link[stats:quantile]{quantile}}
  in order to obtain some values of central
  tendency of the response, i.e., some spread
  of values in the middle.  when \code{imethod = 1}
  to obtain an initial value for the mean
  Some family functions have an argument called
  \code{iprobs.y}, and if so, then these values
  can overwrite \code{gprobs.y}.





% Then the actual search values will be
% \code{unique(sort(c(gshape, 1/gshape)))}, etc.


  }
  \item{eq.mean}{
  Logical.
  Constrain all the means to be equal?
  This type of argument is simpler than
  \code{parallel} because only a single
  \code{TRUE} or \code{FALSE} can be assigned
  and 
  not a formula.  Thus if \code{TRUE} then it
  will be enforced over all variables.


  }
  \item{parallel}{
  A logical, or a simple formula specifying
  which terms have equal/unequal coefficients.
  The formula must be simple, i.e., additive
  with simple main effects terms.  Interactions
  and nesting etc. are not handled.  To handle
  complex formulas use the \code{constraints}
  argument (of \code{\link{vglm}} etc.);
  however, there is a lot more setting up
  involved and things will not be as convenient.



Here are some examples.
1. \code{parallel = TRUE ~ x2 + x5} means
the parallelism assumption
is only applied to \eqn{X_2}, \eqn{X_5} and
the intercept.
2.  \code{parallel = TRUE ~ -1}
and \code{parallel = TRUE ~ 0}
mean the parallelism assumption
is applied to \emph{no} variables at all.
Similarly,
\code{parallel = FALSE ~ -1} and
\code{parallel = FALSE ~ 0}
mean the parallelism assumption
is applied to \emph{all} the variables
including the intercept.
3.  \code{parallel = FALSE ~ x2 - 1}
and \code{parallel = FALSE ~ x2 + 0}
applies the parallelism constraint to all terms
(including the intercept) except for \eqn{X_2}.
4.  \code{parallel = FALSE ~ x2 * x3}
probably will not work. Instead, expand it
out manually to get 
    \code{parallel = FALSE ~ x2 + x3 + x2:x3},
and that should work.
That's because the main formula processes or
expands the \code{"*"}
operator but \code{parallel} does not.
5. To check whether \code{parallel} has done
what was expected, type
\code{coef(fit, matrix = TRUE)} or
\code{constraints(fit)} for confirmation.









  This argument is common in \pkg{VGAM} family
  functions for categorical responses, e.g.,
  \code{\link{cumulative}},  \code{\link{acat}},
  \code{\link{cratio}}, \code{\link{sratio}}.
  For the proportional odds model
  (\code{\link{cumulative}}) having parallel
  constraints applied to each explanatory
  variable (except for the intercepts) means the
  fitted probabilities do not become negative
  or greater than 1. However this parallelism
  or proportional-odds assumption ought to
  be checked.


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
  Some \pkg{VGAM} family functions use
  simulation to obtain an approximate expected
  information matrix (EIM).  For those that
  do, the \code{nsimEIM} argument specifies
  the number of random variates used per
  observation; the mean of \code{nsimEIM} random
  variates is taken.  Thus \code{nsimEIM}
  controls the accuracy and a larger value
  may be necessary if the EIMs are not
  positive-definite.  For intercept-only models
  (\code{y ~ 1)} the value of \code{nsimEIM}
  can be smaller (since the common value used
  is also then taken as the mean over the
  observations), especially if the number of
  observations is large.


  Some \pkg{VGAM} family functions provide
  two algorithms for estimating the EIM.
  If applicable, set \code{nsimEIM = NULL}
  to choose the other algorithm.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2}
  or \code{3} or ... which specifies the
  initialization method for some parameters or
  a specific parameter.  If failure to converge
  occurs try the next higher value, and continue
  until success.
  For example, \code{imethod = 1} might
  be the method of moments, and
  \code{imethod = 2} might be another method.
  If no value of \code{imethod} works then
  it will be necessary to use arguments such
  as \code{isigma}.  For many \pkg{VGAM}
  family functions it is advisable to try this
  argument with all possible values to safeguard
  against problems such as converging to a local
  solution.  \pkg{VGAM} family functions with
  this argument usually correspond to a model
  or distribution that is relatively hard to
  fit successfully, therefore care is needed
  to ensure the global solution is obtained.
  So using all possible values that this
  argument supplies is a good idea.



  \pkg{VGAM} family functions such
  \code{\link{genpoisson2}} recycle
  \code{imethod} to be of length 2 corresponding
  to the 2 parameters. In the future, this
  feature will be extended to other family
  functions to confer more flexibility.



  }
  \item{Form2}{
    Formula.
    Using applied to models with \eqn{M=2}.
    Specifies the terms for \eqn{\eta_2}
    and the other terms belong to \eqn{\eta_1}.
    It is a way to partition the \bold{X}
    matrix into two sets of covariates,
    where they are assigned to each \eqn{\eta_j}
    separately.
    This argument sets up constraint matrices
    \code{rbind(0, 1)}
    for terms in \code{Form2}
    and \code{rbind(1, 0)} for
    \code{setdiff(formula, Form2)}
    so to speak.
    Note that 
    sometimes this argument is only accessed
    if \code{vfl = TRUE}.
    Arguments such as \code{Form1} and
    \code{Form3} are also possible in
    \pkg{VGAM} family functions because
    the \eqn{\eta_j} which is likely to be
    modelled more simply is
    chosen for convenience.


  }
  \item{vfl}{
    A single logical.
    This stands for
    \emph{variance--variance factored loglinear}
    (VFL)
    model.
    If \code{TRUE} then usually
    some other argument such
    as \code{Form2} or \code{parallel}
    is used to partition the main
    \code{\link{vglm}} \code{formula}
    into two sets of covariates.
    For some families
    such as \code{\link{negbinomial}}
    this enables overdispersion to be modelled
    conveniently via a loglinear model,
    given the mean.
    It is necessary to read the online help
    regarding each \pkg{VGAM} family function
    because each one may different from others.
    To fit some VFL models it is necessary to
    make a copy of existing covariates by
    creating new variable names and then
    adding it to the main formula.



    A good question is:
    why is \code{vfl} necessary?
    Wouldn't \code{Form2} be sufficient?
    Setting \code{vfl = TRUE} enables
    some internal checking such
    as avoiding conflicts.
    For example, it is often necessary to set
    \code{zero = NULL} and
    \code{parallel = FALSE}, otherwise
    there would be contradictions.



% Other VFL families include
% \code{\link{uninormal}},
% \code{\link{betabinomial}},
% \code{\link{genpoisson1}},
% \code{\link{genpoisson2}},
% \code{\link{}},
% \code{\link{}},
% \code{\link{}},
% \code{\link{}},
% \code{\link{}},





  }
  \item{type.fitted}{
  Character.
  Type of fitted value returned by
  the \code{fitted()} methods function.
  The first choice is always the default.
  The available choices depends on what kind
  of family function it is.  Using the first
  few letters of the chosen choice is okay.
  See \code{\link{fittedvlm}} for more details.



  The choice \code{"Qlink"} refers to
  quantile-links, which was introduced in
  December 2018 in \pkg{VGAMextra} 0.0-2
  for several 1-parameter distributions.
  Here, either the \code{\link{loglink}}
  or \code{\link{logitlink}} or
  \code{\link{identitylink}} of the quantile
  is the link function (and the choice is
  dependent on the support of the distribution),
  and link functions end in \code{"Qlink"}.
  A limited amount of support is provided
  for such links, e.g., \code{fitted(fit)}
  are the fitted quantiles, which is the same
  as \code{predict(fit, type = "response")}.
  However, \code{fitted(fit, percentiles = 77)}
  will not work.


  }
  \item{percentiles}{
  Numeric vector, with values between 0 and
  100 (although it is not recommended that
  exactly 0 or 100 be inputted).  Used only
  if \code{type.fitted = "quantiles"}
  or \code{type.fitted = "percentiles"},
  then this argument specifies the values of
  these quantiles.  The argument name tries
  to reinforce that the values lie between 0
  and 100.  See \code{\link{fittedvlm}} for
  more details.


  }
  \item{probs.x, probs.y}{
  Numeric, with values in (0, 1).
  The probabilites that define quantiles with
  respect to some vector, usually an \code{x} or
  \code{y} of some sort.  This is used to create
  two subsets of data corresponding to `low'
  and `high' values of x or y.  Each value is
  separately fed into the \code{probs} argument
  of \code{\link[stats:quantile]{quantile}}.
  If the data set size is small then it may be
  necessary to increase/decrease slightly the
  first/second values respectively.


  }
  \item{lss}{
  Logical.
  This stands for the ordering: location,
  scale and shape.  Should the ordering
  of the parameters be in this order?
  Almost all \pkg{VGAM} family functions
  have this order by default, but in order
  to match the arguments of existing R
  functions, one might need to set
  \code{lss = FALSE}.  For example, the arguments of
  \code{\link{weibullR}} are scale and shape,
  whereas \code{\link[stats]{rweibull}} are
  shape and scale.  As a temporary measure
  (from \pkg{VGAM} 0.9-7 onwards but prior to
  version 1.0-0), some family functions such
  as \code{\link{sinmad}} have an \code{lss}
  argument without a default. For these,
  setting \code{lss = FALSE} will work.
  Later, \code{lss = TRUE} will be the default.
  Be careful for the \code{dpqr}-type functions,
  e.g., \code{\link{rsinmad}}.



  }
  \item{Thresh}{
  \emph{Thresholds} is another name for the
  intercepts, e.g., for categorical models.
  They may be constrained by functions such as
  \code{\link{CM.equid}} and
  \code{\link{CM.qnorm}}.
  The string \code{"CM."} and
  the \code{Thresh} argument is pasted and
  then that function is called to obtain the
  constraint matrix for the intercept term.
  So
  \code{Thresh = "free"},
  \code{Thresh = "equid"},
  \code{Thresh = "qnorm"},
  \code{Thresh = "symm0"},
  \code{Thresh = "symm1"} etc.
  are possibilities.
  Families that use this include 
  \code{\link{multinomial}},
  \code{\link{cratio}},
  \code{\link{sratio}},
  \code{\link{cumulative}},
  \code{\link{acat}}.




  }
  \item{whitespace}{
  Logical.
  Should white spaces (\code{" "}) be used in the
  labelling of the linear/additive predictors?
  Setting \code{TRUE} usually results in more
  readability but it occupies more columns of
  the output.


  }
  \item{oim}{
  Logical.
  Should the observed information matrices
  (OIMs) be used for the working weights?
  In general, setting \code{oim = TRUE} means
  the Newton-Raphson algorithm,
  and \code{oim = FALSE} means Fisher-scoring.
  The latter uses the EIM, and is usually recommended.
  If \code{oim = TRUE} then \code{nsimEIM}
  is ignored.


  }
  \item{nrfs}{
    Numeric, a value in \eqn{[0, 1]}.
    Experimental argument for allowing a mixture
    of Newton-Raphson and Fisher scoring.
    The working weights are taken as a linear
    combination of the two.  If \code{nrfs = 0}
    then Newton-Raphson is used, i.e., the OIM
    is wholly used.  If \code{nrfs = 1} then
    Fisher scoring is used, i.e., the EIM is
    wholly used.  If convergence is successful
    then one might expect Newton-Raphson to
    be faster than Fisher scoring because the
    former has an order-2 convergence rate
    while the latter has a linear rate.


  },

  \item{multiple.responses}{
  Logical.
  Some \pkg{VGAM} family functions allow
  a multivariate or vector response.
  If so, then usually the response is a
  matrix with columns corresponding to the
  individual response variables.  They are
  all fitted simultaneously.  Arguments such
  as \code{parallel} may then be useful
  to allow for relationships between the
  regressions of each response variable.
  If \code{multiple.responses = TRUE} then
  sometimes the response is interpreted
  differently, e.g., \code{\link{posbinomial}}
  chooses the first column of a matrix response
  as success and combines the other columns as
  failure, but when
  \code{multiple.responses = TRUE}
  then each column of the response
  matrix is the number of successes and the
  \code{weights} argument is of the same
  dimension as the response and contains the
  number of trials.



  }
  \item{earg.link}{
  This argument should be generally ignored.


% Sometimes the link argument can receive \code{earg}-type input,
% such as \code{\link{quasibinomial}} calling \code{\link{binomial}}.



  }
  \item{byrow.arg}{
  Logical.
  Some \pkg{VGAM} family functions that handle
  multiple responses have arguments that allow
  input to be fed in which affect all the
  responses, e.g., \code{imu} for initalizing
  a \code{mu} parameter.  In such cases it is
  sometime more convenient to input one value
  per response by setting \code{byrow.arg = TRUE};
  then values are recycled in order to
  form a matrix of the appropriate dimension.
  This argument matches \code{byrow} in
  \code{\link[base]{matrix}}; in fact it is
  fed into such using
  \code{matrix(..., byrow = byrow.arg)}.
  This argument has no effect when there is
  one response.


  }
  \item{ishrinkage}{
  Shrinkage factor \eqn{s} used for obtaining
  initial values.
  Numeric, between 0 and 1.
  In general, the formula used is something like
  \eqn{s \mu + (1-s) y}{s*mu + (1-s)*y}
  where \eqn{\mu}{mu} is a measure of
  central tendency such as a weighted mean or
  median, and \eqn{y} is the response vector.
  For example, the initial values are slight
  perturbations of the mean towards the actual
  data.  For many types of models this method
  seems to work well and is often reasonably
  robust to outliers in the response.  Often
  this argument is only used if the argument
  \code{imethod} is assigned a certain value.



  }
  \item{nointercept}{
  An integer-valued vector specifying
  which linear/additive predictors have no
  intercepts.  Any values must be from the
  set \{1,2,\ldots,\eqn{M}\}.  A value of
  \code{NULL} means no such constraints.



  }
  \item{bred}{
  Logical.
  Some \pkg{VGAM} family functions will allow
  bias-reduction based on the work by Kosmidis
  and Firth.  Sometimes half-stepping is a good
  idea; set \code{stepsize = 0.5} and monitor
  convergence by setting \code{trace = TRUE}.



  }
  \item{link.list, earg.list}{
  Some \pkg{VGAM} family functions
  (such as \code{\link{normal.vcm}})
  implement models with
  potentially lots of parameter link functions.
  These two arguments allow many such links and extra arguments
  to be inputted more easily.
  One has something like
  \code{link.list = list}
\code{("(Default)" = "identitylink", x2 = "loglink", x3 = "logofflink")}
  and
  \code{earg.list = list}
  \code{("(Default)" = list(), x2 = list(), x3 = "list(offset = -1)")}.
  Then any unnamed terms will have the default
  link with its corresponding extra argument.
  Note: the \code{\link{multilogitlink}} link
  is also possible, and if so, at least two
  instances of it are necessary.  Then the last
  term is the baseline/reference group.


  }

}
\value{
  An object of class \code{"vglmff"} (see
  \code{\link{vglmff-class}}).  The object
  is used by modelling functions such as
  \code{\link{vglm}} and \code{\link{vgam}}.

  

}
\section{Warning }{
  The \code{zero} argument is supplied for
  convenience but conflicts can arise with other
  arguments, e.g., the \code{constraints}
  argument of \code{\link{vglm}} and
  \code{\link{vgam}}.  See Example 5 below
  for an example.  If not sure, use, e.g.,
  \code{constraints(fit)} and
  \code{coef(fit, matrix = TRUE)}
  to check the result of a fit \code{fit}.


  The arguments \code{zero} and
  \code{nointercept} can be inputted with values
  that fail. For example,
  \code{multinomial(zero = 2, nointercept = 1:3)}
  means the second linear/additive predictor is
  identically zero, which will cause a failure.


  Be careful about the use of other
  potentially contradictory constraints, e.g.,
  \code{multinomial(zero = 2, parallel = TRUE ~ x3)}.
  If in doubt, apply \code{constraints()}
  to the fitted object to check.



  \pkg{VGAM} family functions with the
  \code{nsimEIM} may have inaccurate working
  weight matrices. If so, then the standard
  errors of the regression coefficients
  may be inaccurate. Thus output from
  \code{summary(fit)}, \code{vcov(fit)},
  etc. may be misleading.


  Changes relating to the \code{lss} argument
  have very important consequences and users
  must beware.  Good programming style is
  to rely on the argument names and not on
  the order.



}

\details{
  Full details will be given in documentation
  yet to be written, at a later date!
  A general recommendation is to set
  \code{trace = TRUE} whenever any model fitting
  is done, since monitoring convergence is
  usually very informative.


}

\references{



Yee, T. W. (2015).
Vector Generalized Linear and Additive Models:
With an Implementation in R.
New York, USA: \emph{Springer}.




Kosmidis, I. and Firth, D. (2009).
Bias reduction in exponential family nonlinear models.
\emph{Biometrika},
\bold{96}, 793--804.


%\bold{96}(4), 793--804.


%Kosmidis, I. and Firth, D. (2010).
%A generic algorithm for reducing bias in parametric estimation.
%\emph{Electronic Journal of Statistics},
%\bold{4}, 1097--1112.



Miranda-Soberanis, V. F. and Yee, T. W. (2019).
New link functions for distribution--specific
quantile regression based on vector generalized linear and
additive models.
\emph{Journal of Probability and Statistics},
\bold{5}, 1--11.



}

\seealso{
  \code{\link{Links}},
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{vglmff-class}},
  \code{\link{UtilitiesVGAM}},
  \code{\link{multilogitlink}},
  \code{\link{multinomial}},
  \pkg{VGAMextra}.


% \code{\link{normal.vcm}},


  
}
\author{T. W. Yee}

%\note{
%  See \code{\link{Links}} regarding a major
%  change in link functions,
%  for version 0.9-0 and higher
%  (released during the 2nd half of 2012).



%}

\examples{
# Example 1
cumulative()
cumulative(link = "probitlink", reverse = TRUE, parallel = TRUE)

# Example 2
wdata <- data.frame(x2 = runif(nn <- 1000))
wdata <- transform(wdata,
         y = rweibull(nn, shape = 2 + exp(1 + x2), scale = exp(-0.5)))
fit <- vglm(y ~ x2, weibullR(lshape = logofflink(offset = -2), zero = 2),
            data = wdata)
coef(fit, mat = TRUE)

# Example 3; multivariate (multiple) response
\dontrun{
ndata <- data.frame(x = runif(nn <- 500))
ndata <- transform(ndata,
           y1 = rnbinom(nn, exp(1), mu = exp(3+x)),  # k is size
           y2 = rnbinom(nn, exp(0), mu = exp(2-x)))
fit <- vglm(cbind(y1, y2) ~ x, negbinomial(zero = -2), ndata)
coef(fit, matrix = TRUE)
}
# Example 4
\dontrun{
# fit1 and fit2 are equivalent
fit1 <- vglm(ymatrix ~ x2 + x3 + x4 + x5,
             cumulative(parallel = FALSE ~ 1 + x3 + x5), cdata)
fit2 <- vglm(ymatrix ~ x2 + x3 + x4 + x5,
             cumulative(parallel = TRUE ~ x2 + x4), cdata)
}

# Example 5
udata <- data.frame(x2 = rnorm(nn <- 200))
udata <- transform(udata,
           x1copy = 1,  # Copy of the intercept
           x3 = runif(nn),
           y1 = rnorm(nn, 1 - 3*x2, sd = exp(1 + 0.2*x2)),
           y2 = rnorm(nn, 1 - 3*x2, sd = exp(1)))
args(uninormal)
fit1 <- vglm(y1 ~ x2, uninormal, udata)            # This is okay
fit2 <- vglm(y2 ~ x2, uninormal(zero = 2), udata)  # This is okay
fit4 <- vglm(y2 ~ x2 + x1copy + x3,
             uninormal(zero = NULL, vfl = TRUE,
                       Form2 = ~ x1copy + x3 - 1), udata)
coef(fit4,  matrix = TRUE)  # VFL model

# This creates potential conflict
clist <- list("(Intercept)" = diag(2), "x2" = diag(2))
fit3 <- vglm(y2 ~ x2, uninormal(zero = 2), data = udata,
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

# Example 7 ('zero' argument with character input)
set.seed(123); n <- 1000
ldata <- data.frame(x2 = runif(n))
ldata <- transform(ldata, y1 = rlogis(n, loc = 5*x2, scale = exp(2)))
ldata <- transform(ldata, y2 = rlogis(n, loc = 5*x2, scale = exp(1*x2)))
ldata <- transform(ldata, w1 = runif(n))
ldata <- transform(ldata, w2 = runif(n))
fit7 <- vglm(cbind(y1, y2) ~ x2,
#        logistic(zero = "location1"),  # location1 is intercept-only
#        logistic(zero = "location2"),
#        logistic(zero = "location*"),  # Not okay... all is unmatched
#        logistic(zero = "scale1"),
#        logistic(zero = "scale2"),
#        logistic(zero = "scale"),  # Both scale parameters are matched
         logistic(zero = c("location", "scale2")),  # All but scale1
#        logistic(zero = c("LOCAT", "scale2")),  # Only scale2 is matched
#        logistic(zero = c("LOCAT")),  # Nothing is matched
#        trace = TRUE,
#        weights = cbind(w1, w2),
         weights = w1,
         data = ldata)
coef(fit7, matrix = TRUE)
}

\keyword{models}

