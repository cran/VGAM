\name{vglm}
\alias{vglm}
%\alias{vglm.fit}
\title{Fitting Vector Generalized Linear Models }
\description{
  \code{vglm} is used to fit vector generalized linear models (VGLMs).
  This is a large class of models that includes
  generalized linear models (GLMs) as special cases.

}
\usage{
vglm(formula, family, data = list(), weights = NULL, subset = NULL, 
     na.action = na.fail, etastart = NULL, mustart = NULL, 
     coefstart = NULL, control = vglm.control(...), offset = NULL, 
     method = "vglm.fit", model = FALSE, x.arg = TRUE, y.arg = TRUE, 
     contrasts = NULL, constraints = NULL, extra = list(), 
     qr.arg = FALSE, smart = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  In the following, \eqn{M} is the number of linear predictors.

  \item{formula}{
  a symbolic description of the model to be fit.
  The RHS of the formula is applied to each linear predictor. Different
  variables in each linear predictor can be chosen by specifying
  constraint matrices. 

  }
  \item{family}{
  a function of class \code{"vglmff"} describing
  what statistical model is to be fitted.
  These are called ``\pkg{VGAM} family functions''.

  }
  \item{data}{
  an optional data frame containing the variables in the model.
  By default the variables are taken from
  \code{environment(formula)}, typically the environment from which
  \code{vglm} is called.

  }
  \item{weights}{
  an optional vector or matrix of (prior) weights 
  to be used in the fitting process.
  If \code{weights} is a matrix, then it must be in
  \emph{matrix-band} form, whereby the first \eqn{M} 
  columns of the matrix are the
  diagonals, followed by the upper-diagonal band, followed by the
  band above that, etc. In this case, there can be up to \eqn{M(M+1)}
  columns, with the last column corresponding to the (1,\eqn{M}) elements
  of the weight matrices.

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
  It is a \eqn{M}-column matrix. If \eqn{M=1} then it may be a vector.

  }
  \item{mustart}{
  starting values for the 
  fitted values. It can be a vector or a matrix. 
  Some family functions do not make use of this argument.

  }
  \item{coefstart}{
  starting values for the coefficient vector.

  }
  \item{control}{
  a list of parameters for controlling the fitting process. 
  See \code{\link{vglm.control}} for details.

  }
  \item{offset}{
   a vector or \eqn{M}-column matrix of offset values.  These are \emph{a
   priori} known and are added to the linear predictors during fitting.

  }
  \item{method}{
  the method to be used in fitting the model.  The default (and
  presently only) method \code{vglm.fit} uses iteratively reweighted
  least squares (IRLS).

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
  The components of the list must be named with the term it corresponds
  to (and it must match in character format exactly). 
  Each constraint matrix must have \eqn{M} rows, and be of full-column
  rank. By default, constraint matrices are the \eqn{M} by \eqn{M}
  identity
  matrix unless arguments in the family function itself override
  these values. 
  If \code{constraints} is used it must contain \emph{all} the
  terms; an incomplete list is not accepted.

  }
  \item{extra}{
  an optional list with any extra information that might be needed by
  the \pkg{VGAM} family function.

  }
  \item{qr.arg}{
  logical value indicating whether
  the slot \code{qr}, which returns the QR decomposition of the
  VLM model matrix, is returned on the object.

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
  Here, \eqn{j=1,\ldots,M} where \eqn{M} is finite.
  Then one can write
  \eqn{\eta=(\eta_1,\ldots,\eta_M)^T}{eta=(eta_1,\ldots,\eta_M)^T}
  as a vector of linear predictors.

  Most users will find \code{vglm} similar in flavour to
  \code{\link[stats]{glm}}.  The function \code{vglm.fit} actually does
  the work.

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
  \item{fitted.values}{the fitted values, as a matrix.
  This may be missing or consist entirely
  of \code{NA}s, e.g., the Cauchy model.
  }
  \item{misc}{a list to hold miscellaneous parameters.}
  \item{model}{the model frame.}
  \item{na.action}{a list holding information about missing values.}
  \item{offset}{if non-zero, a \eqn{M}-column matrix of offsets.}
  \item{post}{a list where post-analysis results may be put.}
  \item{preplot}{used by \code{\link{plotvgam}}, the plotting parameters
        may be put here.}
  \item{prior.weights}{initially supplied weights.}
  \item{qr}{the QR decomposition used in the fitting.}
  \item{R}{the \bold{R} matrix in the QR decomposition used in the fitting.}
  \item{rank}{numerical rank of the fitted model.}
  \item{residuals}{the \emph{working} residuals at the final IRLS iteration.}
  \item{rss}{residual sum of squares at the final IRLS iteration with
  the adjusted dependent vectors and weight matrices.}
  \item{smart.prediction}{
  a list of data-dependent parameters (if any)
  that are used by smart prediction.

  }
  \item{terms}{the \code{\link[stats]{terms}} object used.}
  \item{weights}{the weight matrices at the final IRLS iteration.
    This is in matrix-band form.}
  \item{x}{the model matrix (linear model LM, not VGLM).}
  \item{xlevels}{the levels of the factors, if any, used in fitting.}
  \item{y}{the response, in matrix form.}

  This slot information is repeated at \code{\link{vglm-class}}.
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

The \pkg{VGAM} library can be downloaded
starting from \url{http://www.stat.auckland.ac.nz/~yee}.
Other \pkg{VGAM} resources and documentation can be found there. 
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

  The smart prediction (\code{\link{smartpred}}) library is packed with
  the \pkg{VGAM} library.

  The theory behind the scaling parameter is currently being made more
  rigorous, but it it should give the same value as the scale parameter
  for GLMs.

  In Example 5 below, the \code{xij} argument to illustrate covariates
  that are specific to a linear predictor. Here, \code{lop}/\code{rop} are
  the ocular pressures of the left/right eye (artificial data). Variables
  \code{leye} and \code{reye} might be the presence/absence of a particular
  disease on the LHS/RHS eye respectively.  See \code{\link{fill}}
  for more details and examples.

}

%~Make other sections like WARNING with \section{WARNING }{....} ~

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
  \code{coef.vlm},
  \code{\link{predict.vglm}},
  \code{summary.vglm},
  etc.
}

\examples{
# Example 1. Dobson (1990) Page 93: Randomized Controlled Trial :
counts = c(18,17,15,20,10,20,25,13,12)
outcome = gl(3,1,9)
treatment = gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
vglm.D93 = vglm(counts ~ outcome + treatment, family=poissonff)
summary(vglm.D93)


# Example 2. Multinomial logit model
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
vglm(cbind(normal, mild, severe) ~ let, multinomial, pneumo)


# Example 3. Proportional odds model
fit = vglm(cbind(normal,mild,severe) ~ let, cumulative(par=TRUE), pneumo)
coef(fit, matrix=TRUE) 
constraints(fit) 
fit@x # LM model matrix
model.matrix(fit) # Larger VGLM model matrix


# Example 4. Bivariate logistic model 
data(coalminers)
fit = vglm(cbind(nBnW, nBW, BnW, BW) ~ age, binom2.or, coalminers, trace=TRUE)
coef(fit, matrix=TRUE)
fit@y


# Example 5. The use of the xij argument
n = 1000
eyes = data.frame(lop = runif(n), rop = runif(n)) 
eyes = transform(eyes, 
                 leye = ifelse(runif(n) < logit(-1+2*lop, inverse=TRUE), 1, 0),
                 reye = ifelse(runif(n) < logit(-1+2*rop, inverse=TRUE), 1, 0))
fit = vglm(cbind(leye,reye) ~ lop + rop + fill(lop),
           binom2.or(exchangeable=TRUE, zero=3),
           xij = op ~ lop + rop + fill(lop), data=eyes)
coef(fit)
coef(fit, matrix=TRUE)
coef(fit, matrix=TRUE, compress=FALSE)

# Here's one method to handle the xij argument with a term that
# produces more than one column in the model matrix. 
POLY3 = function(x, ...) {
    # A cubic 
    poly(c(x,...), 3)[1:length(x),]
}

fit = vglm(cbind(leye,reye) ~ POLY3(lop,rop) + POLY3(rop,lop) + fill(POLY3(lop,rop)),
           binom2.or(exchangeable=TRUE, zero=3),  data=eyes,
           xij = POLY3(op) ~ POLY3(lop,rop) + POLY3(rop,lop) + 
                             fill(POLY3(lop,rop)))
coef(fit)
coef(fit, matrix=TRUE)
coef(fit, matrix=TRUE, compress=FALSE)
predict(fit)[1:4,]
} 
\keyword{models}
\keyword{regression}

%eyes$leye = ifelse(runif(n) < 1/(1+exp(-1+2*eyes$lop)), 1, 0)
%eyes$reye = ifelse(runif(n) < 1/(1+exp(-1+2*eyes$rop)), 1, 0)
