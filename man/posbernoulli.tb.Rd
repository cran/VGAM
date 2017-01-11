\name{posbernoulli.tb}
%\alias{posbernoulli}
\alias{posbernoulli.tb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Bernoulli Family Function with Time and Behavioural Effects }
\description{
  Fits a GLM/GAM-like model to multiple
  Bernoulli responses where
  each row in the capture history matrix response has at least one success
  (capture).
  Sampling occasion effects and behavioural effects are accommodated.


}
\usage{
posbernoulli.tb(link = "logit", parallel.t = FALSE ~ 1,
                parallel.b = FALSE ~ 0, drop.b = FALSE ~ 1,
                type.fitted = c("likelihood.cond", "mean.uncond"),
                imethod = 1, iprob = NULL,
                p.small = 1e-4, no.warning = FALSE,
                ridge.constant = 0.01, ridge.power = -4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, imethod, iprob}{
  See \code{\link{CommonVGAMffArguments}} for information.

  }



  \item{parallel.t, parallel.b, drop.b}{
  A logical, or formula with a logical as the response.
  See \code{\link{CommonVGAMffArguments}} for information.
  The \code{parallel.}-type arguments
  specify whether the constraint matrices have a parallelism
  assumption for the temporal and behavioural effects.
  Argument \code{parallel.t} means parallel with respect to time, and
  matches the same argument name in \code{\link{posbernoulli.t}}.


  Suppose the model is intercept-only.
  Setting \code{parallel.t = FALSE ~ 0} results in the \eqn{M_b} model.
  Setting \code{drop.b = FALSE ~ 0} results in the \eqn{M_t} model
  because it drops columns off the constraint matrices corresponding
  to any behavioural effect.
  Setting \code{parallel.t = FALSE ~ 0} and
  setting \code{parallel.b = FALSE ~ 0} results in the \eqn{M_b} model.
  Setting \code{parallel.t = FALSE ~ 0},
          \code{parallel.b = FALSE ~ 0} and
          \code{drop.b = FALSE ~ 0} results in the \eqn{M_0} model.
  Note the default for \code{parallel.t} and \code{parallel.b}
  may be unsuitable for data sets which have a large \eqn{\tau}
  because of the large number of parameters; it might be too flexible.
  If it is desired to have the behaviour affect some of
  the other covariates then set \code{drop.b = TRUE ~ 0}.


  The default model has a different intercept for each
  sampling occasion, a time-parallelism assumption for all other covariates,
  and a dummy variable representing a
  single behavioural effect (also in the intercept).


  The most flexible model is to set
  \code{parallel.b = TRUE  ~ 0},
  \code{parallel.t = TRUE  ~ 0} and
  \code{drop.b = TRUE ~ 0}.
  This means that all possible temporal and behavioural effects are
  estimated, for the intercepts and other covariates.
  Such a model is \emph{not} recommended; it will contain a lot of paramters.


  }
  \item{type.fitted}{
  Character, one of the choices for the type of fitted value returned.
  The default is the first one.
  Partial matching is okay.
  For \code{"likelihood.cond"}:
  the probability defined by the conditional likelihood.
  For \code{"mean.uncond"}: the unconditional mean, which should
  agree with \code{\link[base]{colMeans}} applied to the response
  matrix for intercept-only models.


  }
  \item{ridge.constant, ridge.power}{
  Determines the ridge parameters at each IRLS iteration.
  They are the constant and power (exponent) for the ridge adjustment
  for the working weight matrices (the capture probability block
  matrix, hence the first \eqn{\tau} diagonal values).
  At iteration \eqn{a} of the IRLS algorithm
  a positive value is added to the first \eqn{\tau}{tau}
  diagonal elements of the working weight matrices to make
  them positive-definite. This adjustment is the
  mean of the diagonal elements of \code{wz} multipled by
  \eqn{K \times a^p}{K * a^p}
  where \eqn{K} is \code{ridge.constant} and \eqn{p} is \code{ridge.power}.
  This is always positive but decays to zero as iterations proceed
  (provided \eqn{p} is negative etc.).


  }
  \item{p.small, no.warning}{
  See \code{\link{posbernoulli.t}}.


  }


}
\details{
  This model
  (commonly known as \eqn{M_{tb}}/\eqn{M_{tbh}} in the capture--recapture literature)
  operates on a response matrix of 0s and 1s (\eqn{n \times \tau}{n x tau}).
  See \code{\link{posbernoulli.t}}
  for information that is in common.
  It allows time and behavioural effects to be modelled.


  Evidently,
  the expected information matrix (EIM) seems \emph{not}
  of full rank (especially in early iterations), so
  \code{ridge.constant} and \code{ridge.power} are used to
  \emph{try} fix up the problem.
  The default link functions
  are \eqn{(logit \,p_{c1},\ldots,logit \,p_{c\tau},logit \,p_{r2},\ldots,logit \,p_{r\tau})^T}{
           (logit p_{c1},\ldots,logit p_{c,tau},logit p_{r2},\ldots,logit p_{r,tau})^T}
  where the subscript \eqn{c} denotes capture,
        the subscript \eqn{r} denotes recapture,
  and it is not possible to recapture the animal at sampling occasion 1.
  Thus \eqn{M = 2\tau - 1}{M=2*tau-1}.
  The parameters are currently prefixed by \code{pcapture} and \code{precapture}
  for the capture and recapture probabilities.
  This \pkg{VGAM} family function may be further modified in the future.


% Not surprisingly,
% the fitted values are similar to \code{\link{posbernoulli.t}} and
% \code{\link{posbernoulli.b}}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

  See \code{\link{posbernoulli.t}}.


}
\author{ Thomas W. Yee. }

\note{
  It is a good idea to apply the parallelism assumption to each sampling
  occasion except possibly with respect to the intercepts.
  Also, a simple behavioural effect such as being modelled using the intercept
  is recommended; if the behavioural effect is not parallel and/or
  allowed to apply to other covariates
  then there will probably be too many parameters, and hence,
  numerical problems. See \code{M_tbh.1} below.


%Data-wise, at each sampling occasion, the \eqn{M_{tb}} model requires at least
%one first capture and at least one noncapture.


% If not all of the \eqn{2^{\tau}-1}{2^(tau) - 1} combinations of
% the response matrix are not present then it pays to add
% such rows to the response matrix and assign a small but
% positive prior weight.
% For example, if \eqn{\tau=2}{tau=2} then there should be
% (0,1) rows,
% (1,0) rows and
% (1,1) rows present in the response matrix.


  It is a good idea to monitor convergence.
  Simpler models such as the \eqn{M_0}/\eqn{M_h} models
  are best fitted with \code{\link{posbernoulli.t}} or
  \code{\link{posbernoulli.b}} or
  \code{\link{posbinomial}}.



% yettodo:
% Some time in the future it might be possible to allow for a
% different tau value for each row.
% Then the response would be a matrix padded with NAs on the RHS.



}


\seealso{
  \code{\link{posbernoulli.b}} (including \code{N.hat}),
  \code{\link{posbernoulli.t}},
  \code{\link{posbinomial}},
  \code{\link{Select}},
  \code{\link{fill1}},
  \code{\link{Huggins89table1}},
  \code{\link{Huggins89.t1}},
  \code{\link{deermice}},
  \code{\link{prinia}}.


}

\examples{
\dontrun{
# Example 1: simulated data
nTimePts <- 5  # (aka tau == # of sampling occasions)
nnn <- 1000   # Number of animals
pdata <- rposbern(n = nnn, nTimePts = nTimePts, pvars = 2)
dim(pdata); head(pdata)

M_tbh.1 <- vglm(cbind(y1, y2, y3, y4, y5) ~ x2,
                posbernoulli.tb, data = pdata, trace = TRUE)
coef(M_tbh.1)  # First element is the behavioural effect
coef(M_tbh.1, matrix = TRUE)
constraints(M_tbh.1, matrix = TRUE)
summary(M_tbh.1, presid = FALSE)  # Standard errors are approximate
head(fitted(M_tbh.1))
head(model.matrix(M_tbh.1, type = "vlm"), 21)
dim(depvar(M_tbh.1))

M_tbh.2 <- vglm(cbind(y1, y2, y3, y4, y5) ~ x2,
                posbernoulli.tb(parallel.t = FALSE ~  0),
                data = pdata, trace = TRUE)
coef(M_tbh.2)  # First element is the behavioural effect
coef(M_tbh.2, matrix = TRUE)
constraints(M_tbh.2, matrix = TRUE)
summary(M_tbh.2, presid = FALSE)  # Standard errors are approximate
head(fitted(M_tbh.2))
head(model.matrix(M_tbh.2, type = "vlm"), 21)
dim(depvar(M_tbh.2))

# Example 2: deermice subset data
fit1 <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight,
             posbernoulli.t, data = deermice, trace = TRUE)
coef(fit1)
coef(fit1, matrix = TRUE)
constraints(fit1, matrix = TRUE)
summary(fit1, presid = FALSE)  # Standard errors are approximate

# fit1 is the same as Fit1 (a M_{th} model):
Fit1 <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight,
             posbernoulli.tb(drop.b = TRUE ~ sex + weight,
                             parallel.t = TRUE),  # No parallelism for the intercept
             data = deermice, trace = TRUE)
constraints(Fit1)
}
}
\keyword{models}
\keyword{regression}




%\section{Warning }{
%  As this model is likely to be overparameterized, probably this
%  function should not be used (for now?).
%
%% From Jakub:
%  Estimation for the population size (and its SE) for the
%  \eqn{M_{tb}} and \eqn{M_{tbh}} model may be wrong.
%  But models
%  \eqn{M_{0}},
%  \eqn{M_{h}},
%  \eqn{M_{b}},
%  \eqn{M_{bh}},
%  \eqn{M_{t}},
%  \eqn{M_{th}}
%  seem fine.
%
%  Inference, especially using standard errors, may be fraught here
%  because the EIM is, strictly speaking, not of full rank.
%  A similar adjustment is made by \code{\link{zipebcom}}.
%  It is a good idea to monitor convergence.
%  The \eqn{M_0}/\eqn{M_h} models are best fitted with
%  \code{\link{posbernoulli.t}} or \code{\link{posbinomial}} because
%  the standard errors are more accurate.
%
%
%}

%yyy <- depvar(fit1)
%if (length(table(4 * yyy[, 1] + 2 * yyy[, 2] + 1 * yyy[, 3])) != 2^(ncol(yyy)) - 1)
%  warning("not every combination is represented by a row in the response matrix")


