\name{cao}
\alias{cao}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting Constrained Additive Ordination (CAO) }
\description{
  A constrained additive ordination (CAO) model is fitted using the
  \emph{reduced-rank vector generalized additive model} (RR-VGAM)
  framework.

}
\usage{
cao(formula, family, data = list(),
    weights = NULL, subset = NULL, na.action = na.fail,
    etastart = NULL, mustart = NULL, coefstart = NULL,
    control = cao.control(...), offset = NULL,
    method = "cao.fit", model = FALSE, x.arg = TRUE, y.arg = TRUE,
    contrasts = NULL, constraints = NULL,
    extra = NULL, qr.arg = FALSE, smart = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{formula}{
    a symbolic description of the model to be fit.  The RHS of the
    formula is used to construct the latent variables, upon which the
    smooths are applied.  All the variables in the formula are used
    for the construction of latent variables except for those specified
    by the argument \code{noRRR}, which is itself a formula.  The LHS
    of the formula contains the response variables, which should be a
    matrix with each column being a response (species).


  }
  \item{family}{
  a function of class \code{"vglmff"} (see \code{\link{vglmff-class}})
  describing what statistical model is to be fitted. This is called a
  ``\pkg{VGAM} family function''.  See \code{\link{CommonVGAMffArguments}}
  for general information about many types of arguments found in this
  type of function.
    See \code{\link{cqo}} for a list of those presently implemented.


  }

  \item{data}{
    an optional data frame containing the variables in the model.
    By default the variables are taken from \code{environment(formula)},
    typically the environment from which \code{cao} is called.


  }
  \item{weights}{
    an optional vector or matrix of (prior) weights to be used in the
    fitting process.  For \code{cao}, this argument currently should
    not be used.


  }
  \item{subset}{
    an optional logical vector specifying a subset of observations to
    be used in the fitting process.


  }
  \item{na.action}{
    a function which indicates what should happen when the data contain
    \code{NA}s.  The default is set by the \code{na.action} setting of
    \code{\link[base]{options}}, and is \code{na.fail} if that is unset.
    The ``factory-fresh'' default is \code{na.omit}.


  }
  \item{etastart}{
    starting values for the linear predictors.  It is a \eqn{M}-column
    matrix. If \eqn{M=1} then it may be a vector.  For \code{cao},
    this argument currently should not be used.


    }
  \item{mustart}{
    starting values for the fitted values. It can be a vector or a
    matrix.  Some family functions do not make use of this argument.
    For \code{cao}, this argument currently should not be used.


  }
  \item{coefstart}{
    starting values for the coefficient vector.  For \code{cao}, this
    argument currently should not be used.


  }
  \item{control}{
    a list of parameters for controlling the fitting process.
    See \code{\link{cao.control}} for details.


  }
  \item{offset}{
    a vector or \eqn{M}-column matrix of offset values.  These are
    \emph{a priori} known and are added to the linear predictors during
    fitting.  For \code{cao}, this argument currently should not be used.


  }
  \item{method}{
    the method to be used in fitting the model.  The default
    (and presently only) method \code{cao.fit} uses iteratively
    reweighted least squares (IRLS) within FORTRAN code called from
    \code{\link[stats]{optim}}.


  }
  \item{model}{
    a logical value indicating whether the \emph{model frame} should be
    assigned in the \code{model} slot.


  }
  \item{x.arg, y.arg}{
    logical values indicating whether the model matrix and response
    vector/matrix used in the fitting process should be assigned in the
    \code{x} and \code{y} slots.  Note the model matrix is the linear
    model (LM) matrix.


  }
  \item{contrasts}{
    an optional list. See the \code{contrasts.arg} of
    \code{\link{model.matrix.default}}.


  }
  \item{constraints}{
    an optional list  of constraint matrices.  For \code{cao}, this
    argument currently should not be used.  The components of the list
    must be named with the term it corresponds to (and it must match
    in character format).  Each constraint matrix must have \eqn{M}
    rows, and be of full-column rank. By default, constraint matrices
    are the \eqn{M} by \eqn{M} identity matrix unless arguments in the
    family function itself override these values.  If \code{constraints}
    is used it must contain \emph{all} the terms; an incomplete list is
    not accepted.


  }
  \item{extra}{
    an optional list with any extra information that might be needed by
    the family function.  For \code{cao}, this argument currently should
    not be used.


  }
  \item{qr.arg}{
    For \code{cao}, this argument currently should not be used.


  }
  \item{smart}{
    logical value indicating whether smart prediction
    (\code{\link{smartpred}}) will be used.


  }
  \item{\dots}{
    further arguments passed into \code{\link{cao.control}}.


  }
}
\details{
  The arguments of \code{cao} are a mixture of those from
  \code{\link{vgam}} and \code{\link{cqo}}, but with some extras
  in \code{\link{cao.control}}. Currently, not all of the
  arguments work properly.


  CAO can be loosely be thought of as the result of fitting generalized
  additive models (GAMs) to several responses (e.g., species) against
  a very small number of latent variables.  Each latent variable is a
  linear combination of the explanatory variables; the coefficients
  \bold{C} (called \eqn{C} below) are called \emph{constrained
  coefficients} or \emph{canonical coefficients}, and are interpreted as
  weights or loadings. The \bold{C} are estimated by maximum likelihood
  estimation.  It is often a good idea to apply \code{\link[base]{scale}}
  to each explanatory variable first.


  For each response (e.g., species), each latent variable is smoothed
  by a cubic smoothing spline, thus CAO is data-driven. If each smooth
  were a quadratic then CAO would simplify to \emph{constrained quadratic
  ordination} (CQO; formerly called \emph{canonical Gaussian ordination}
  or CGO).
  If each smooth were linear then CAO would simplify to \emph{constrained
  linear ordination} (CLO). CLO can theoretically be fitted with
  \code{cao} by specifying \code{df1.nl=0}, however it is more efficient
  to use \code{\link{rrvglm}}.


  Currently, only \code{Rank=1} is implemented, and only
  \code{noRRR = ~1} models are handled.


% Poisson and binary responses are implemented (viz.,
% \code{\link{poissonff}}, \code{\link{binomialff}}), and
% dispersion parameters for these must be assumed known.  Hence using
% \code{\link{quasipoissonff}} and \code{\link{quasibinomialff}} will
% currently fail.  Also, currently, only \code{noRRR = ~ 1} models are
% handled.


  With binomial data, the default formula is
  \deqn{logit(P[Y_s=1]) = \eta_s = f_s(\nu), \ \ \ s=1,2,\ldots,S}{%
        logit(P[Y_s=1]) =  eta_s = f_s(\nu), \ \ \ s=1,2,\ldots,S}
  where \eqn{x_2}{x_2} is a vector of environmental variables, and
  \eqn{\nu=C^T x_2}{nu=C^T x_2} is a \eqn{R}-vector of latent variables.
  The \eqn{\eta_s}{eta_s} is an additive predictor for species \eqn{s},
  and it models the probabilities of presence as an additive model on
  the logit scale.  The matrix \eqn{C} is estimated from the data, as
  well as the smooth functions \eqn{f_s}.  The argument \code{noRRR = ~
  1} specifies that the vector \eqn{x_1}{x_1}, defined for RR-VGLMs
  and QRR-VGLMs, is simply a 1 for an intercept.
  Here, the intercept in the model is absorbed into the functions.
  A \code{\link{cloglog}} link may be preferable over a
  \code{\link{logit}} link.


  With Poisson count data, the formula is
  \deqn{\log(E[Y_s]) = \eta_s = f_s(\nu)}{%
         log(E[Y_s]) =  eta_s = f_s(\nu)}
  which models the mean response as an additive models on the log scale.


  The fitted latent variables (site scores) are scaled to have
  unit variance.  The concept of a tolerance is undefined for
  CAO models, but the optimums and maximums are defined. The generic
  functions \code{\link{Max}} and \code{\link{Opt}} should work for
  CAO objects, but note that if the maximum occurs at the boundary then
  \code{\link{Max}} will return a \code{NA}.  Inference for CAO models
  is currently undeveloped.


}
\value{
  An object of class \code{"cao"}
  (this may change to \code{"rrvgam"} in the future).
  Several generic functions can be applied to the object, e.g.,
  \code{\link{Coef}}, \code{\link{concoef}}, \code{\link{lvplot}},
  \code{\link{summary}}.


}

\references{
Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.



% Documentation accompanying the \pkg{VGAM} package at
% \url{http://www.stat.auckland.ac.nz/~yee}
% contains further information and examples.



}
\author{T. W. Yee}
\note{
  CAO models are computationally expensive, therefore
  setting \code{trace = TRUE} is a good idea, as well as running it
  on a simple random sample of the data set instead.


  Sometimes the IRLS algorithm does not converge within the FORTRAN
  code. This results in warnings being issued.  In particular, if an
  error code of 3 is issued, then this indicates the IRLS algorithm has
  not converged. One possible remedy is to increase or decrease
  the nonlinear degrees of freedom so that the curves become more or
  less flexible, respectively.


}
\section{Warning }{
  CAO is very costly to compute. With version 0.7-8 it took 28 minutes on
  a fast machine. I hope to look at ways of speeding things up in the
  future.


  Use \code{\link[base:Random]{set.seed}} just prior to calling
  \code{cao()} to make your results reproducible.
  The reason for this is finding the optimal
  CAO model presents a difficult optimization problem, partly because the
  log-likelihood function contains many local solutions. To obtain the
  (global) solution the user is advised to try \emph{many} initial values.
  This can be done by setting \code{Bestof} some appropriate value
  (see \code{\link{cao.control}}). Trying many initial values becomes
  progressively more important as the nonlinear degrees of freedom of
  the smooths increase.


% The code is a little fragile at this stage, so the function might
% hang/lock up in the microsoft Windows version.


  Currently the dispersion parameter for a
  \code{\link{gaussianff}} CAO model is estimated slightly differently
  and may be slightly biassed downwards (usually a little too small).


}

\seealso{
  \code{\link{cao.control}},
  \code{Coef.cao},
  \code{\link{cqo}},
  \code{\link{latvar}},
  \code{\link{Opt}},
  \code{\link{Max}},
  \code{persp.cao},
  \code{\link{poissonff}},
  \code{\link{binomialff}},
  \code{\link{negbinomial}},
  \code{\link{gamma2}},
  \code{\link{gaussianff}},
  \code{\link[base:Random]{set.seed}},
  \code{\link[gam]{gam}},
  \code{\link[VGAMdata]{trapO}}.

}

\examples{
\dontrun{
hspider[, 1:6] <- scale(hspider[, 1:6])  # Standardized environmental vars
set.seed(149)  # For reproducible results
ap1 <- cao(cbind(Pardlugu, Pardmont, Pardnigr, Pardpull) ~
           WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
           family = poissonff, data = hspider, Rank = 1,
           df1.nl = c(Pardpull= 2.7, 2.5),
           Bestof = 7, Crow1positive = FALSE)
sort(deviance(ap1, history = TRUE))  # A history of all the iterations

Coef(ap1)
concoef(ap1)

par(mfrow = c(2, 2))
plot(ap1)  # All the curves are unimodal; some quite symmetric

par(mfrow = c(1, 1), las = 1)
index <- 1:ncol(depvar(ap1))
lvplot(ap1, lcol = index, pcol = index, y = TRUE)

trplot(ap1, label = TRUE, col = index)
abline(a = 0, b = 1, lty = 2)

trplot(ap1, label = TRUE, col = "blue", log = "xy", which.sp = c(1, 3))
abline(a = 0, b = 1, lty = 2)

persp(ap1, col = index, lwd = 2, label = TRUE)
abline(v = Opt(ap1), lty = 2, col = index)
abline(h = Max(ap1), lty = 2, col = index)
}
}
\keyword{models}
\keyword{regression}


