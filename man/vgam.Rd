\name{vgam}
\alias{vgam}
%\alias{vgam.fit}
\title{ Fitting Vector Generalized Additive Models }
% 20030215; This file is based a lot from vglm.Rd
\description{
  Fit a vector generalized additive model (VGAM).
  Both 1st-generation VGAMs (based on backfitting) and
  2nd-generation VGAMs (based on P-splines, with automatic
  smoothing parameter selection) are implemented.
  This is a large class
  of models that includes generalized additive models (GAMs) and vector
  generalized linear models (VGLMs) as special cases.


}
\usage{
vgam(formula, family, data = list(), weights = NULL, subset = NULL,
     na.action = na.fail, etastart = NULL, mustart = NULL,
     coefstart = NULL, control = vgam.control(...), offset = NULL,
     method = "vgam.fit", model = FALSE, x.arg = TRUE, y.arg = TRUE,
     contrasts = NULL, constraints = NULL,
     extra = list(), form2 = NULL, qr.arg = FALSE, smart = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  % The following comes from vglm.Rd but with minor tweaks

  \item{formula}{
  a symbolic description of the model to be fit.
  The RHS of the formula is applied to each linear/additive predictor,
  and should include at least one
     \code{\link[VGAM]{sm.os}} term
  or \code{\link[VGAM]{sm.ps}} term
  or \code{\link[VGAM]{s}} term.
  Mixing both together is not allowed.
  Different variables in each linear/additive predictor
  can be chosen by specifying constraint matrices.


  }
  \item{family}{
  Same as for \code{\link{vglm}}.


  }
  \item{data}{
  an optional data frame containing the variables in the model.
  By default the variables are taken from
  \code{environment(formula)}, typically the environment from which
  \code{vgam} is called.


  }
  \item{weights, subset, na.action}{
  Same as for \code{\link{vglm}}.
  Note that \code{subset} may be unreliable and to get around
  this problem it is best to use \code{\link[base]{subset}} to create
  a new smaller data frame and feed in the smaller data frame.
  See below for an example.
  This is a bug that needs fixing.


  }
  \item{etastart, mustart, coefstart}{
  Same as for \code{\link{vglm}}.


  }
  \item{control}{
  a list of parameters for controlling the fitting process.
  See \code{\link{vgam.control}} for details.


  }
  \item{method}{
  the method to be used in fitting the model.
  The default (and presently only) method \code{vgam.fit}
  uses iteratively reweighted least squares (IRLS).


  }
  \item{constraints, model, offset}{
  Same as for \code{\link{vglm}}.


  }
  \item{x.arg, y.arg}{
  logical values indicating whether the model matrix and response
  vector/matrix used in the fitting process should be assigned in the
  \code{x} and \code{y} slots.  Note the model matrix is the LM model
  matrix; to get the VGAM model matrix type \code{model.matrix(vgamfit)}
  where \code{vgamfit} is a \code{vgam} object.


  }
  \item{contrasts, extra, form2, qr.arg, smart}{
  Same as for \code{\link{vglm}}.


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
  by smoothers.
  The first term in the summation is just the intercept.
  Currently
  two types of smoothers are
  implemented:
  \code{\link[VGAM]{s}} represents
  the older and more traditional one, called a
  \emph{vector (cubic smoothing spline) smoother} and is
  based on Yee and Wild (1996);
  it is more similar to the \R{} package \pkg{gam}.
  The newer one is represented by
  \code{\link[VGAM]{sm.os}} and
  \code{\link[VGAM]{sm.ps}}, and these are
  based on O-splines and P-splines---they allow automatic
  smoothing parameter selection; it is more similar
  to the \R{} package \pkg{mgcv}.




  In the above, \eqn{j=1,\ldots,M} where \eqn{M} is finite.
  If all the functions are constrained to be linear then
  the resulting model is a vector generalized linear model
  (VGLM).  VGLMs are best fitted with \code{\link{vglm}}.




  Vector (cubic smoothing spline) smoothers are represented
  by \code{s()} (see \code{\link[VGAM]{s}}). Local
  regression via \code{lo()} is \emph{not} supported. The
  results of \code{vgam} will differ from the \code{gam()}
  (in the \pkg{gam}) because \code{vgam()} uses a different
  knot selection algorithm. In general, fewer knots are
  chosen because the computation becomes expensive when
  the number of additive predictors \eqn{M} is large.



  Second-generation VGAMs are based on the
  O-splines and P-splines.
  The latter is due to Eilers and Marx (1996).
  Backfitting is not required, and estimation is performed using IRLS.
  The function \code{\link{sm.os}} represents a \emph{smart}
  implementation of O-splines.
  The function \code{\link{sm.ps}} represents a \emph{smart}
  implementation of P-splines.
  Written G2-VGAMs or P-VGAMs, this methodology should not be used
  unless the sample size is reasonably large.
  Usually an UBRE predictive criterion is optimized
  (at each IRLS iteration)
  because the
  scale parameter for VGAMs is usually assumed to be known.
  This search for optimal smoothing parameters does not always converge,
  and neither is it totally reliable.
  G2-VGAMs implicitly set \code{criterion = "coefficients"} so that
  convergence occurs when the change in the regression coefficients
  between 2 IRLS iterations is sufficiently small.
  Otherwise the search for the optimal smoothing parameters might
  cause the log-likelihood to decrease between 2 IRLS iterations.
  Currently \emph{outer iteration} is implemented,
  by default,
  rather than \emph{performance iteration} because the latter
  is more easy to converge to a local solution; see
  Wood (2004) for details.
  One can use \emph{performance iteration}
  by setting \code{Maxit.outer = 1} in
  \code{\link{vgam.control}}.


% outeriter



%  A suggested rule-of-thumb is at least 500 observations.



  The underlying algorithm of VGAMs is IRLS.
  First-generation VGAMs (called G1-VGAMs)
  are estimated by modified vector backfitting
  using vector splines. O-splines are used as the basis functions
  for the vector (smoothing) splines, which are a lower dimensional
  version of natural B-splines.
  The function \code{vgam.fit()} actually does the
  work.  The smoothing code is based on F. O'Sullivan's
  BART code.






%  If more than one of \code{etastart}, \code{start} and \code{mustart}
%  is specified, the first in the list will be used.


  A closely related methodology based on VGAMs called
  \emph{constrained additive ordination} (CAO) first forms
  a linear combination of the explanatory variables (called
  \emph{latent variables}) and then fits a GAM to these.
  This is implemented in the function \code{\link{cao}}
  for a very limited choice of family functions.


}
\value{
  For G1-VGAMs and G2-VGAMs, an object of class
  \code{"vgam"} or
  \code{"pvgam"}
  respectively
  (see \code{\link{vgam-class}}
   and \code{\link{pvgam-class}}
  for further information).


}
\references{



Wood, S. N. (2004).
Stable and efficient multiple smoothing parameter estimation
for generalized additive models.
\emph{J. Amer. Statist. Assoc.}, \bold{99}(467): 673--686.



Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.



Yee, T. W. (2008)
The \code{VGAM} Package.
\emph{R News}, \bold{8}, 28--39.



Yee, T. W. (2015)
Vector Generalized Linear and Additive Models:
With an Implementation in R.
New York, USA: \emph{Springer}.



Yee, T. W. (2016).
Comments on ``Smoothing parameter and model selection for
general smooth models''
by Wood, S. N. and Pya, N. and Safken, N.,
\emph{J. Amer. Statist. Assoc.}, \bold{110}(516).



%Yee, T. W. and Somchit, C. and Wild, C. J. (2016)
%Generation-2
%vector generalized additive models.
%Manuscript in preparation.



%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.



%Wood, S. N. (2004).
%Stable and efficient multiple smoothing parameter estimation
%for generalized additive models.
%\emph{J. Amer. Statist. Assoc.}, \bold{99}(467): 673--686.



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



  Some \pkg{VGAM} family functions end in \code{"ff"}
  to avoid interference with other functions, e.g.,
  \code{\link{binomialff}}, \code{\link{poissonff}},
  \code{\link{gaussianff}}, \code{gammaff}. This is
  because \pkg{VGAM} family functions are incompatible with
  \code{\link[stats]{glm}} (and also \code{\link[gam]{gam}}
  in the \pkg{gam} library and \code{\link[mgcv]{gam}}
  in the \pkg{mgcv} library).


  The smart prediction (\code{\link{smartpred}}) library
  is packed with the \pkg{VGAM} library.


  The theory behind the scaling parameter is currently being
  made more rigorous, but it it should give the same value
  as the scale parameter for GLMs.


}



%~Make other sections like WARNING with \section{WARNING }{....} ~
\section{WARNING}{
  For G1-VGAMs,
  currently \code{vgam} can only handle constraint matrices \code{cmat},
  say, such that \code{crossprod(cmat)} is diagonal.
  It can be detected by \code{\link{is.buggy}}.
  VGAMs with constraint matrices that have non-orthogonal columns should
  be fitted with
  \code{\link{sm.os}} or
  \code{\link{sm.ps}} terms
  instead of \code{\link{s}}.


%  This is a bug that I will try to fix up soon;



  See warnings in \code{\link{vglm.control}}.


}

\seealso{
  \code{\link{is.buggy}},
  \code{\link{vgam.control}},
  \code{\link{vgam-class}},
  \code{\link{vglmff-class}},
  \code{\link{plotvgam}},
  \code{\link{summaryvgam}},
  \code{\link{summarypvgam}},
  \code{\link{sm.os}},
  \code{\link{sm.ps}},
  \code{\link[VGAM]{s}},
  \code{\link[mgcv]{magic}}.
  \code{\link{vglm}},
  \code{\link{vsmooth.spline}},
  \code{\link{cao}}.


}

\examples{# Nonparametric proportional odds model
pneumo <- transform(pneumo, let = log(exposure.time))
vgam(cbind(normal, mild, severe) ~ s(let),
     cumulative(parallel = TRUE), data = pneumo, trace = TRUE)

# Nonparametric logistic regression
hfit <- vgam(agaaus ~ s(altitude, df = 2), binomialff, data = hunua)
\dontrun{ plot(hfit, se = TRUE) }
phfit <- predict(hfit, type = "terms", raw = TRUE, se = TRUE)
names(phfit)
head(phfit$fitted)
head(phfit$se.fit)
phfit$df
phfit$sigma

# Fit two species simultaneously
hfit2 <- vgam(cbind(agaaus, kniexc) ~ s(altitude, df = c(2, 3)),
              binomialff(multiple.responses = TRUE), data = hunua)
coef(hfit2, matrix = TRUE)  # Not really interpretable
\dontrun{
plot(hfit2, se = TRUE, overlay = TRUE, lcol = 3:4, scol = 3:4)
ooo <- with(hunua, order(altitude))
with(hunua, matplot(altitude[ooo], fitted(hfit2)[ooo,], ylim = c(0, 0.8),
     xlab = "Altitude (m)", ylab = "Probability of presence", las = 1,
     main = "Two plant species' response curves", type = "l", lwd = 2))
with(hunua, rug(altitude))
}

# The 'subset' argument does not work here. Use subset() instead.
set.seed(1)
zdata <- data.frame(x2 = runif(nn <- 500))
zdata <- transform(zdata, y = rbinom(nn, 1, 0.5))
zdata <- transform(zdata, subS = runif(nn) < 0.7)
sub.zdata <- subset(zdata, subS)  # Use this instead
if (FALSE)
  fit4a <- vgam(cbind(y, y) ~ s(x2, df = 2),
                binomialff(multiple.responses = TRUE),
                data = zdata, subset = subS)  # This fails!!!
fit4b <- vgam(cbind(y, y) ~ s(x2, df = 2),
              binomialff(multiple.responses = TRUE),
              data = sub.zdata)  # This succeeds!!!
fit4c <- vgam(cbind(y, y) ~ sm.os(x2),
              binomialff(multiple.responses = TRUE),
              data = sub.zdata)  # This succeeds!!!
\dontrun{par(mfrow = c(2, 2))
plot(fit4b, se = TRUE, shade = TRUE, shcol = "pink")
plot(fit4c, se = TRUE, shade = TRUE, shcol = "pink")
}
}
\keyword{models}
\keyword{regression}
\keyword{smooth}



