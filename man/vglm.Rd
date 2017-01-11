\name{vglm}
\alias{vglm}
%\alias{vglm.fit}
\title{Fitting Vector Generalized Linear Models }
\description{
  \code{vglm} is used to fit vector generalized linear models (VGLMs).
  This is a very large class of models that includes
  generalized linear models (GLMs) as a special case.


}
\usage{
vglm(formula, family, data = list(), weights = NULL, subset = NULL,
     na.action = na.fail, etastart = NULL, mustart = NULL,
     coefstart = NULL, control = vglm.control(...), offset = NULL,
     method = "vglm.fit", model = FALSE, x.arg = TRUE, y.arg = TRUE,
     contrasts = NULL, constraints = NULL, extra = list(),
     form2 = NULL, qr.arg = TRUE, smart = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{

  \item{formula}{
  a symbolic description of the model to be fit.
  The RHS of the formula is applied to each linear
  predictor. Different variables in each linear predictor
  can be chosen by specifying constraint matrices.


  }
  \item{family}{
  a function of class \code{"vglmff"} (see \code{\link{vglmff-class}})
  describing what statistical model is to be fitted. This is called a
  ``\pkg{VGAM} family function''.  See \code{\link{CommonVGAMffArguments}}
  for general information about many types of arguments found in this
  type of function.
  The argument name \code{"family"} is used loosely and for
  the ease of existing \code{\link[stats]{glm}} users;
  there is no concept of a formal ``error distribution'' for VGLMs.
  Possibly the argument name should be better \code{"model"}
  but unfortunately
  that name has already been taken.



  }
  \item{data}{
  an optional data frame containing the variables in the model.
  By default the variables are taken from
  \code{environment(formula)}, typically the environment
  from which \code{vglm} is called.




  }
  \item{weights}{
  an optional vector or matrix of (prior fixed and known) weights
  to be used in the fitting process.
  If the \pkg{VGAM} family function handles multiple responses
  (\eqn{Q > 1} of them, say) then
  \code{weights} can be a matrix with \eqn{Q} columns.
  Each column matches the respective response.
  If it is a vector (the usually case) then it is recycled into a
  matrix with \eqn{Q} columns.
  The values of \code{weights} must be positive; try setting
  a very small value such as \code{1.0e-8} to effectively
  delete an observation.




% 20140507:

  Currently the \code{weights} argument does not support sampling
  weights from complex sampling designs.
  And currently sandwich estimators are not computed in any shape or form.
  The present weights are multiplied by the corresponding
  log-likelihood contributions.






% If \code{weights} is a matrix,
% then it should be must be in \emph{matrix-band} form, whereby the
% first \eqn{M} columns of the matrix are the diagonals,
% followed by the upper-diagonal band, followed by the
% band above that, etc. In this case, there can be up to
% \eqn{M(M+1)} columns, with the last column corresponding
% to the (1,\eqn{M}) elements of the weight matrices.


  }
  \item{subset}{
  an optional logical vector specifying a subset of
  observations to
  be used in the fitting process.



  }
  \item{na.action}{
  a function which indicates what should happen when
  the data contain \code{NA}s.
  The default is set by the \code{na.action} setting
  of \code{\link[base]{options}}, and is \code{na.fail} if that is unset.
  The ``factory-fresh'' default is \code{na.omit}.


  }
  \item{etastart}{
  starting values for the linear predictors.
  It is a \eqn{M}-column matrix with the same number of rows as the response.
  If \eqn{M = 1} then it may be a vector.
  Note that \code{etastart} and the output of \code{predict(fit)}
  should be comparable.
  Here, \code{fit} is the fitted object.


  }
  \item{mustart}{
  starting values for the fitted values.
  It can be a vector or a matrix;
  if a matrix, then it has the same number of rows as the response.
  Usually \code{mustart} and the output of \code{fitted(fit)}
  should be comparable.
  Some family functions do not make use of this argument.


  }
  \item{coefstart}{
  starting values for the coefficient vector.
  The length and order must match that of \code{coef(fit)}.


  }
  \item{control}{
  a list of parameters for controlling the fitting process.
  See \code{\link{vglm.control}} for details.


  }
  \item{offset}{
   a vector or \eqn{M}-column matrix of offset values.
   These are \emph{a priori} known and are added to the
   linear/additive predictors during fitting.


  }
  \item{method}{
  the method to be used in fitting the model.  The default (and
  presently only) method \code{vglm.fit()} uses iteratively
  reweighted least squares (IRLS).


  }
  \item{model}{
  a logical value indicating whether the
  \emph{model frame}
  should be assigned in the \code{model} slot.


  }
  \item{x.arg, y.arg}{
  logical values indicating whether
  the model matrix and response vector/matrix used in the fitting
  process should be assigned in the \code{x} and \code{y} slots.
  Note the model matrix is the LM model matrix; to get the VGLM
  model matrix type \code{model.matrix(vglmfit)} where
  \code{vglmfit} is a \code{vglm} object.


  }
  \item{contrasts}{
  an optional list. See the \code{contrasts.arg}
  of \code{\link{model.matrix.default}}.


  }
  \item{constraints}{
  an optional list  of constraint matrices.
  The components of the list must be named with the term it corresponds to
  (and it must match in character format exactly).
  There are two types of input: \code{"lm"}-type and \code{"vlm"}-type.
  The former is a subset of the latter.
  The former has a matrix for each term of the LM matrix.
  The latter has a matrix for each column of the VLM matrix.
  After fitting, the \code{\link{constraints}}
  extractor function may be applied; it returns
  the \code{"vlm"}-type list of constraint matrices
  by default.  If \code{"lm"}-type are returned by
  \code{\link{constraints}} then these can be fed into this
  argument and it should give the same model as before.


  Each constraint matrix must have \eqn{M} rows, and be of
  full-column rank.  By default, constraint matrices are
  the \eqn{M} by \eqn{M} identity matrix unless arguments
  in the family function itself override these values, e.g.,
  \code{parallel} (see  \code{\link{CommonVGAMffArguments}}).
  If \code{constraints} is used it must contain \emph{all}
  the terms; an incomplete list is not accepted.


  }
  \item{extra}{
  an optional list with any extra information that might be needed by
  the \pkg{VGAM} family function.


  }
  \item{form2}{
  The second (optional) formula.
  If argument \code{xij} is used (see \code{\link{vglm.control}}) then
  \code{form2} needs to have \emph{all} terms in the model.
  Also, some \pkg{VGAM} family functions such as \code{\link{micmen}}
  use this argument to input the regressor variable.
  If given, the slots \code{@Xm2} and \code{@Ym2} may be assigned.
  Note that smart prediction applies to terms in \code{form2} too.


  }
  \item{qr.arg}{
  logical value indicating whether the slot \code{qr}, which
  returns the QR decomposition of the VLM model matrix,
  is returned on the object.


  }
  \item{smart}{
  logical value indicating whether smart prediction
  (\code{\link{smartpred}}) will be used.


  }
  \item{\dots}{
  further arguments passed into \code{\link{vglm.control}}.


  }

}
\details{
  A vector generalized linear model (VGLM) is loosely defined
  as a statistical model that is a function of \eqn{M} linear
  predictors.
  The central formula is given by
  \deqn{\eta_j = \beta_j^T x}{%
         eta_j = beta_j^T x}
  where \eqn{x}{x} is a vector of explanatory variables
  (sometimes just a 1 for an intercept),
  and
  \eqn{\beta_j}{beta_j} is a vector of regression coefficients
  to be estimated.
  Here, \eqn{j=1,\ldots,M}, where \eqn{M} is finite.
  Then one can write
  \eqn{\eta=(\eta_1,\ldots,\eta_M)^T}{eta=(eta_1,\ldots,\eta_M)^T}
  as a vector of linear predictors.


  Most users will find \code{vglm} similar in flavour to
  \code{\link[stats]{glm}}.
  The function \code{vglm.fit} actually does the work.

% If more than one of \code{etastart}, \code{start} and \code{mustart}
% is specified, the first in the list will be used.

}
\value{
  An object of class \code{"vglm"}, which has the
  following slots. Some of these may not be assigned to save
  space, and will be recreated if necessary later.
  \item{extra}{the list \code{extra} at the end of fitting.}
  \item{family}{the family function (of class \code{"vglmff"}).}
  \item{iter}{the number of IRLS iterations used.}
  \item{predictors}{a \eqn{M}-column matrix of linear predictors.}
  \item{assign}{a named list which matches the columns and the
       (LM) model matrix terms.}
  \item{call}{the matched call.}
  \item{coefficients}{a named vector of coefficients.}
  \item{constraints}{
    a named list of constraint matrices used in the fitting.
  }
  \item{contrasts}{the contrasts used (if any).}
  \item{control}{list of control parameter used in the fitting.}
  \item{criterion}{list of convergence criterion evaluated at the
                   final IRLS iteration.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{df.total}{the total degrees of freedom.}
  \item{dispersion}{the scaling parameter.}
  \item{effects}{the effects.}
  \item{fitted.values}{
  the fitted values, as a matrix.
  This is often the mean but may be quantiles, or the location
  parameter, e.g., in the Cauchy model.

  }
  \item{misc}{a list to hold miscellaneous parameters.}
  \item{model}{the model frame.}
  \item{na.action}{a list holding information about missing values.}
  \item{offset}{if non-zero, a \eqn{M}-column matrix of offsets.}
  \item{post}{a list where post-analysis results may be put.}
  \item{preplot}{used by \code{\link{plotvgam}}, the plotting parameters
        may be put here.}
  \item{prior.weights}{
  initially supplied weights
  (the \code{weights} argument).
  Also see \code{\link{weightsvglm}}.

  }
  \item{qr}{the QR decomposition used in the fitting.}
  \item{R}{the \bold{R} matrix in the QR decomposition used in the fitting.}
  \item{rank}{numerical rank of the fitted model.}
  \item{residuals}{the \emph{working} residuals at the final IRLS iteration.}
  \item{ResSS}{residual sum of squares at the final IRLS iteration with
  the adjusted dependent vectors and weight matrices.}
  \item{smart.prediction}{
  a list of data-dependent parameters (if any)
  that are used by smart prediction.

  }
  \item{terms}{the \code{\link[stats]{terms}} object used.}
  \item{weights}{the working weight matrices at the final IRLS iteration.
    This is in matrix-band form.}
  \item{x}{the model matrix (linear model LM, not VGLM).}
  \item{xlevels}{the levels of the factors, if any, used in fitting.}
  \item{y}{the response, in matrix form.}


  This slot information is repeated at \code{\link{vglm-class}}.


}
\references{


Yee, T. W. (2015)
Vector Generalized Linear and Additive Models:
With an Implementation in R.
New York, USA: \emph{Springer}.



Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.



  Yee, T. W. (2014)
  Reduced-rank vector generalized linear models with two linear predictors.
  \emph{Computational Statistics and Data Analysis},
  \bold{71}, 889--902.



Yee, T. W. (2008)
The \code{VGAM} Package.
\emph{R News}, \bold{8}, 28--39.


%  Documentation accompanying the \pkg{VGAM} package at
%  \url{http://www.stat.auckland.ac.nz/~yee}
%  contains further information and examples.


}


\author{ Thomas W. Yee }
\note{
  This function can fit a wide variety of statistical models. Some of
  these are harder to fit than others because of inherent numerical
  difficulties associated with some of them. Successful model fitting
  benefits from cumulative experience. Varying the values of arguments
  in the \pkg{VGAM} family function itself is a good first step if
  difficulties arise, especially if initial values can be inputted.
  A second, more general step, is to vary the values of arguments in
  \code{\link{vglm.control}}.
  A third step is to make use of arguments such as \code{etastart},
  \code{coefstart} and \code{mustart}.


  Some \pkg{VGAM} family functions end in \code{"ff"} to avoid
  interference with other functions, e.g., \code{\link{binomialff}},
  \code{\link{poissonff}}, \code{\link{gaussianff}},
  \code{gammaff}. This is because \pkg{VGAM} family
  functions are incompatible with \code{\link[stats]{glm}}
  (and also \code{\link[gam]{gam}} in the \pkg{gam} library and
  \code{\link[mgcv]{gam}} in the \pkg{mgcv} library).


  The smart prediction (\code{\link{smartpred}}) library is incorporated
  within the \pkg{VGAM} library.


  The theory behind the scaling parameter is currently being made more
  rigorous, but it it should give the same value as the scale parameter
  for GLMs.


  In Example 5 below, the \code{xij} argument to illustrate covariates
  that are specific to a linear predictor. Here, \code{lop}/\code{rop} are
  the ocular pressures of the left/right eye (artificial data).
  Variables \code{leye} and \code{reye} might be the presence/absence of
  a particular disease on the LHS/RHS eye respectively.
  See
  \code{\link{vglm.control}}
  and
  \code{\link{fill}}
  for more details and examples.


}

%~Make other sections like WARNING with \section{WARNING }{....} ~
\section{WARNING}{
  See warnings in \code{\link{vglm.control}}.
  Also, see warnings under \code{weights} above regarding
  sampling weights from complex sampling designs.


}



\seealso{
  \code{\link{vglm.control}},
  \code{\link{vglm-class}},
  \code{\link{vglmff-class}},
  \code{\link{smartpred}},
  \code{vglm.fit},
  \code{\link{fill}},
  \code{\link{rrvglm}},
  \code{\link{vgam}}.
  Methods functions include
  \code{\link{coefvlm}},
  \code{\link{constraints.vlm}},
  \code{\link{hatvaluesvlm}},
  \code{\link{linkfun.vglm}},
  \code{\link{plotvglm}},
  \code{\link{predictvglm}},
  \code{summary.vglm},
  \code{AIC.vglm},
  \code{\link{lrtest_vglm}},
  etc.


}

\examples{
# Example 1. See help(glm)
print(d.AD <- data.frame(treatment = gl(3, 3),
                         outcome = gl(3, 1, 9),
                         counts = c(18,17,15,20,10,20,25,13,12)))
vglm.D93 <- vglm(counts ~ outcome + treatment, family = poissonff,
                 data = d.AD, trace = TRUE)
summary(vglm.D93)


# Example 2. Multinomial logit model
pneumo <- transform(pneumo, let = log(exposure.time))
vglm(cbind(normal, mild, severe) ~ let, multinomial, data = pneumo)


# Example 3. Proportional odds model
fit3 <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo)
coef(fit3, matrix = TRUE)
constraints(fit3)
model.matrix(fit3, type = "lm")  # LM model matrix
model.matrix(fit3)               # Larger VGLM (or VLM) model matrix


# Example 4. Bivariate logistic model
fit4 <- vglm(cbind(nBnW, nBW, BnW, BW) ~ age, binom2.or, coalminers)
coef(fit4, matrix = TRUE)
depvar(fit4)  # Response are proportions
weights(fit4, type = "prior")


# Example 5. The use of the xij argument (simple case).
# The constraint matrix for 'op' has one column.
nn <- 1000
eyesdat <- round(data.frame(lop = runif(nn),
                            rop = runif(nn),
                             op = runif(nn)), digits = 2)
eyesdat <- transform(eyesdat, eta1 = -1 + 2 * lop,
                              eta2 = -1 + 2 * lop)
eyesdat <- transform(eyesdat,
           leye = rbinom(nn, size = 1, prob = logit(eta1, inverse = TRUE)),
           reye = rbinom(nn, size = 1, prob = logit(eta2, inverse = TRUE)))
head(eyesdat)
fit5 <- vglm(cbind(leye, reye) ~ op,
             binom2.or(exchangeable = TRUE, zero = 3),
             data = eyesdat, trace = TRUE,
             xij = list(op ~ lop + rop + fill(lop)),
             form2 = ~  op + lop + rop + fill(lop))
coef(fit5)
coef(fit5, matrix = TRUE)
constraints(fit5)
}
\keyword{models}
\keyword{regression}

%eyesdat$leye <- ifelse(runif(n) < 1/(1+exp(-1+2*eyesdat$lop)), 1, 0)
%eyesdat$reye <- ifelse(runif(n) < 1/(1+exp(-1+2*eyesdat$rop)), 1, 0)
%coef(fit, matrix = TRUE, compress = FALSE)



% 20090506 zz Put these examples elsewhere:
%
%# Example 6. The use of the xij argument (complex case).
%# Here is one method to handle the xij argument with a term that
%# produces more than one column in the model matrix.
%# The constraint matrix for 'op' has essentially one column.
%POLY3 <- function(x, ...) {
%    # A cubic; ensures that the basis functions are the same.
%    poly(c(x,...), 3)[1:length(x),]
%    head(poly(c(x,...), 3), length(x), drop = FALSE)
%}
%
%fit6 <- vglm(cbind(leye, reye) ~ POLY3(op), trace = TRUE,
%            fam = binom2.or(exchangeable = TRUE, zero=3),  data=eyesdat,
%            xij = list(POLY3(op) ~ POLY3(lop,rop) + POLY3(rop,lop) +
%                                   fill(POLY3(lop,rop))),
%            form2 = ~  POLY3(op) + POLY3(lop,rop) + POLY3(rop,lop) +
%                       fill(POLY3(lop,rop)))
%coef(fit6)
%coef(fit6, matrix = TRUE)
%head(predict(fit6))
%\dontrun{
%plotvgam(fit6, se = TRUE)  # Wrong since it plots against op, not lop.
%}
%
%
%# Example 7. The use of the xij argument (simple case).
%# Each constraint matrix has 4 columns.
%ymat <- rdiric(n <- 1000, shape=c(4,7,3,1))
%mydat <- data.frame(x1=runif(n), x2=runif(n), x3=runif(n), x4=runif(n),
%                   z1=runif(n), z2=runif(n), z3=runif(n), z4=runif(n),
%                   X2=runif(n), Z2=runif(n))
%mydat <- round(mydat, dig=2)
%fit7 <- vglm(ymat ~ X2 + Z2, data=mydat, crit="c",
%           fam = dirichlet(parallel = TRUE),  # Intercept is also parallel.
%           xij = list(Z2 ~ z1 + z2 + z3 + z4,
%                      X2 ~ x1 + x2 + x3 + x4),
%           form2 =  ~ Z2 + z1 + z2 + z3 + z4 +
%                      X2 + x1 + x2 + x3 + x4)
%head(model.matrix(fit7, type="lm"))   # LM model matrix
%head(model.matrix(fit7, type="vlm"))  # Big VLM model matrix
%coef(fit7)
%coef(fit7, matrix = TRUE)
%max(abs(predict(fit7)-predict(fit7, new=mydat)))  # Predicts correctly
%summary(fit7)




