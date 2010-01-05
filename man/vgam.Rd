\name{vgam}
\alias{vgam}
%\alias{vgam.fit}
\title{ Fitting Vector Generalized Additive Models }
% 15/2/03; based a lot from vglm.Rd 
\description{
  Fit a vector generalized additive model (VGAM).  This is a large class
  of models that includes generalized additive models (GAMs) and vector
  generalized linear models (VGLMs) as special cases.

}
\usage{
vgam(formula, family, data = list(), weights = NULL, subset = NULL, 
     na.action = na.fail, etastart = NULL, mustart = NULL, 
     coefstart = NULL, control = vgam.control(...), offset = NULL, 
     method = "vgam.fit", model = FALSE, x.arg = TRUE, y.arg = TRUE, 
     contrasts = NULL, constraints = NULL, 
     extra = list(), qr.arg = FALSE, smart = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  % The following comes from vglm.Rd but with minor tweaks

  \item{formula}{
  a symbolic description of the model to be fit.
  The RHS of the formula is applied to each linear/additive predictor.
  Different
  variables in each linear/additive predictor can be chosen by specifying
  constraint matrices.

  }
  \item{family}{
  a function of class \code{"vglmff"} (see \code{\link{vglmff-class}})
  describing what statistical model is to be fitted. This is called a
  ``\pkg{VGAM} family function''.  See \code{\link{CommonVGAMffArguments}}
  for general information about many types of arguments found in this
  type of function.

  }
  \item{data}{
  an optional data frame containing the variables in the model.
  By default the variables are taken from
  \code{environment(formula)}, typically the environment from which
  \code{vgam} is called.

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
  starting values for the linear/additive predictors.
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
  See \code{\link{vgam.control}} for details.

  }
  \item{offset}{
  a vector or \eqn{M}-column matrix of offset values.
  These are \emph{a priori} known and are added to the linear/additive
  predictors during fitting.

  }
  \item{method}{
  the method to be used in fitting the model.
  The default (and presently only) method \code{vgam.fit}
  uses iteratively reweighted least squares (IRLS).

  }
  \item{model}{
  a logical value indicating whether the \emph{model frame} should be
  assigned in the \code{model} slot.

  }
  \item{x.arg, y.arg}{
  logical values indicating whether the model matrix and response
  vector/matrix used in the fitting process should be assigned in the
  \code{x} and \code{y} slots.  Note the model matrix is the LM model
  matrix; to get the VGAM model matrix type \code{model.matrix(vgamfit)}
  where \code{vgamfit} is a \code{vgam} object.

  }
  \item{contrasts}{
  an optional list. See the \code{contrasts.arg} of
  \code{\link{model.matrix.default}}.

  }
  \item{constraints}{
  an optional list  of constraint matrices.  The components of the list
  must be named with the term it corresponds to (and it must match in
  character format exactly).  Each constraint matrix must have \eqn{M} rows, and
  be of full-column rank. By default, constraint matrices are the \eqn{M}
  by \eqn{M} identity matrix unless arguments in the family function
  itself override these values.  If \code{constraints} is used it must
  contain \emph{all} the terms; an incomplete list is not accepted.

  }
  \item{extra}{
  an optional list with any extra information that might be needed by
  the \pkg{VGAM} family function.

  }
  \item{qr.arg}{
  logical value indicating whether the slot \code{qr}, which returns
  the QR decomposition of the VLM model matrix, is returned on the object.

  }
  \item{smart}{
  logical value indicating whether smart prediction
  (\code{\link{smartpred}}) will be used.

  }
  \item{\dots}{
  further arguments passed into \code{\link{vgam.control}}.

  }

}
\details{
  A vector generalized additive model (VGAM) is loosely defined
  as a statistical model that is a function of \eqn{M} additive predictors.
  The central formula is given by
  \deqn{\eta_j = \sum_{k=1}^p f_{(j)k}(x_k)}{%
         eta_j = sum_{k=1}^p f_{(j)k}(x_k)}
  where \eqn{x_k}{x_k} is the \eqn{k}th explanatory variable
  (almost always \eqn{x_1=1} for the intercept term),
  and
  \eqn{f_{(j)k}} are smooth functions of \eqn{x_k} that are estimated
  by smoothers. The first term in the summation is just the intercept.
  Currently only one type of smoother is
  implemented and this is called a \emph{vector (cubic smoothing spline)
  smoother}.
  Here, \eqn{j=1,\ldots,M} where \eqn{M} is finite.
  If all the functions are constrained to be linear then the resulting
  model is a vector generalized linear model (VGLM).
  VGLMs are best fitted with \code{\link{vglm}}.

  Vector (cubic smoothing spline) smoothers are represented
  by \code{s()} (see \code{\link[VGAM]{s}}).
  Local regression via \code{lo()} is \emph{not}
  supported. The results of \code{vgam} will differ from the S-PLUS and \R
  \code{gam} function (in the \pkg{gam} \R package) because \code{vgam}
  uses a different knot selection algorithm. In general, fewer knots
  are chosen because the computation becomes expensive when the number
  of additive predictors \eqn{M} is large.

  The underlying algorithm of VGAMs is iteratively
  reweighted least squares (IRLS) and modified vector backfitting
  using vector splines. B-splines are used as the basis functions
  for the vector (smoothing) splines. 
  \code{vgam.fit} is the function that actually does the work.
  The smoothing code is based on F. O'Sullivan's BART code. 

%  If more than one of \code{etastart}, \code{start} and \code{mustart}
%  is specified, the first in the list will be used.

  A closely related methodology based on VGAMs called
  \emph{constrained additive ordination} (CAO)
  first forms a linear combination of the explanatory variables 
  (called \emph{latent variables}) and then fits a GAM to these.
  This is implemented in the function \code{\link{cao}} for a very
  limited choice of family functions.

}
\value{
  An object of class \code{"vgam"}
  (see \code{\link{vgam-class}} for further information). 
}
\references{ 
Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.

Yee, T. W. (2008)
The \code{VGAM} Package.
\emph{R News}, \bold{8}, 28--39.

  Documentation accompanying the \pkg{VGAM} package at
  \url{http://www.stat.auckland.ac.nz/~yee}
  contains further information and examples.

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
  \code{\link{vgam.control}}.
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

}



%~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{
\code{\link{vgam.control}},
\code{\link{vgam-class}},
\code{\link{vglmff-class}},
\code{\link{plotvgam}},
\code{\link{vglm}},
\code{\link[VGAM]{s}},
\code{\link{vsmooth.spline}},
\code{\link{cao}}.
}

\examples{
# Nonparametric proportional odds model 
pneumo = transform(pneumo, let=log(exposure.time))
vgam(cbind(normal,mild,severe) ~ s(let), cumulative(par=TRUE), pneumo)

# Nonparametric logistic regression 
fit = vgam(agaaus ~ s(altitude, df=2), binomialff, hunua)
\dontrun{
plot(fit, se=TRUE)
}

# Fit two species simultaneously 
fit2 = vgam(cbind(agaaus, kniexc) ~ s(altitude, df=c(2,3)),
            binomialff(mv=TRUE), hunua)
coef(fit2, mat=TRUE)   # Not really interpretable 
\dontrun{
plot(fit2, se=TRUE, overlay=TRUE, lcol=1:2, scol=1:2)

ooo = with(hunua, order(altitude))
with(hunua, matplot(altitude[ooo], fitted(fit2)[ooo,], type="l", lwd=2,
     xlab="Altitude (m)", ylab="Probability of presence", las=1,
     main="Two plant species' response curves", ylim=c(0,.8)))
with(hunua, rug(altitude))
}
}
\keyword{models}
\keyword{regression}
\keyword{smooth}
