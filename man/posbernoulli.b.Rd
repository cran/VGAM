\name{posbernoulli.b}
%\alias{posbernoulli}
\alias{posbernoulli.b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Bernoulli Family Function with Behavioural Effects }
\description{
  Fits a GLM-/GAM-like model to multiple Bernoulli responses where
  each row in the capture history matrix response has at least one success
  (capture).
  Capture history behavioural effects are accommodated.


}
\usage{
posbernoulli.b(link = "logit", drop.b = FALSE ~ 1,
               type.fitted = c("likelihood.cond", "mean.uncond"), I2 = FALSE,
               ipcapture = NULL, iprecapture = NULL,
               p.small = 1e-4, no.warning = FALSE)
}
%- maybe also 'usage' for other objects documented here.
% apply.parint = FALSE,
\arguments{
  \item{link, drop.b, ipcapture, iprecapture}{

  See \code{\link{CommonVGAMffArguments}} for information about
  these arguments.
  By default the parallelism assumption does not apply to the intercept.
  With an intercept-only model
  setting \code{drop.b = TRUE ~ 1} results in the \eqn{M_0}/\eqn{M_h} model.


% it just deletes the 2nd column of the constraint matrix corresponding
% to the intercept.


% The default value of \code{zero} means that the behavioural
% effect is modelled as the difference between the
% two intercepts.


% That is, it is modelled through the intercept, and a
% negative value of the second linear/additive predictor means trap shy, etc.



  }
  \item{I2}{
  Logical.
  This argument is used for terms that are not parallel.
  If \code{TRUE} then the constraint matrix \code{diag(2)}
  (the general default constraint matrix in \pkg{VGAM}) is used,
  else \code{cbind(0:1, 1)}. The latter means the first element/column
  corresponds to the behavioural effect.
  Consequently it and its standard error etc. can be accessed
  directly without subtracting two quantities.


  }
  \item{type.fitted}{
  Details at \code{\link{posbernoulli.tb}}.


  }
  \item{p.small, no.warning}{
  See \code{\link{posbernoulli.t}}.


  }


}
\details{
  This model
  (commonly known as \eqn{M_b}/\eqn{M_{bh}} in the
  capture--recapture literature)
  operates on a capture history matrix response of 0s and 1s
  (\eqn{n \times \tau}{n x tau}).
  See \code{\link{posbernoulli.t}} for details,
  e.g., common assumptions with other models.
  Once an animal is captured for the first time,
  it is marked/tagged so that its future
  capture history can be recorded. The effect of the recapture
  probability is modelled through a second linear/additive predictor.
  It is well-known that some species of animals are affected by capture,
  e.g., trap-shy or trap-happy. This \pkg{VGAM} family function
  \emph{does} allow the capture history to be modelled via such
  behavioural effects.
  So does \code{\link{posbernoulli.tb}} but \code{\link{posbernoulli.t}} cannot.



% If \code{drop.b = TRUE} the parallelism does not apply to the intercept.


  The number of linear/additive predictors is \eqn{M = 2},
  and the default links
  are \eqn{(logit \,p_c, logit \,p_r)^T}{(logit p_c, logit p_r)^T}
  where \eqn{p_c} is the probability of capture and
        \eqn{p_r} is the probability of recapture.
  The fitted value returned is of the same dimension as
  the response matrix, and depends on the capture history:
  prior to being first captured, it is \code{pcapture}.
  Afterwards, it is \code{precapture}.


  By default, the constraint matrices for the intercept term
  and the other covariates are set up so that \eqn{p_r}
  differs from \eqn{p_c} by a simple binary effect,
  on a logit scale.
  However, this difference (the behavioural effect) is more
  directly estimated by having \code{I2 = FALSE}.
  Then it allows an estimate of the trap-happy/trap-shy effect;
  these are positive/negative values respectively.
  If \code{I2 = FALSE} then
  the (nonstandard) constraint matrix used is \code{cbind(0:1, 1)},
  meaning the first element can be interpreted as the behavioural effect.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}

%\section{Warning }{
%
%  See \code{\link{posbernoulli.t}}.
%
%
%}

\references{


  See \code{\link{posbernoulli.t}}.


}
\author{ Thomas W. Yee. }

\note{
  The dependent variable is \emph{not} scaled to row proportions.
  This is the same as \code{\link{posbernoulli.t}}
  and \code{\link{posbernoulli.tb}}
  but different from \code{\link{posbinomial}}
  and \code{\link{binomialff}}.


% Monitor convergence by setting \code{trace = TRUE}.


% To  fit \eqn{M_{tb}}{M_tb} and \eqn{M_{tbh}}{M_tbh}
% use \code{\link{posbernoulli.t}} with the \code{xij}
% argument of \code{\link{vglm.control}}.



}

\seealso{
  \code{\link{posbernoulli.t}} and
  \code{\link{posbernoulli.tb}} (including estimating \eqn{N}),
  \code{\link{deermice}},
  \code{\link{dposbern}},
  \code{\link{rposbern}},
  \code{\link{posbinomial}},
  \code{\link{aux.posbernoulli.t}},
  \code{\link{prinia}}.


% \code{\link{huggins91}}.
% \code{\link{vglm.control}} for \code{xij},


}

\examples{
# deermice data ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Fit a M_b model
M.b <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1,
            posbernoulli.b, data = deermice, trace = TRUE)
coef(M.b)["(Intercept):1"]  # Behavioural effect on the logit scale
coef(M.b, matrix = TRUE)
constraints(M.b, matrix = TRUE)
summary(M.b, presid = FALSE)

# Fit a M_bh model
M.bh <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight,
             posbernoulli.b, data = deermice, trace = TRUE)
coef(M.bh, matrix = TRUE)
coef(M.bh)["(Intercept):1"]  # Behavioural effect on the logit scale
constraints(M.bh)  # (2,1) element of "(Intercept)" is for the behavioural effect
summary(M.bh, presid = FALSE)  # Significant positive (trap-happy) behavioural effect
# Approx. 95 percent confidence for the behavioural effect:
SE.M.bh <- coef(summary(M.bh))["(Intercept):1", "Std. Error"]
coef(M.bh)["(Intercept):1"] + c(-1, 1) * 1.96 * SE.M.bh

# Fit a M_h model
M.h <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight,
            posbernoulli.b(drop.b = TRUE ~ sex + weight),
            data = deermice, trace = TRUE)
coef(M.h, matrix = TRUE)
constraints(M.h, matrix = TRUE)
summary(M.h, presid = FALSE)

# Fit a M_0 model
M.0 <- vglm(cbind(    y1 + y2 + y3 + y4 + y5 + y6,
                  6 - y1 - y2 - y3 - y4 - y5 - y6) ~ 1,
            posbinomial, data = deermice, trace = TRUE)
coef(M.0, matrix = TRUE)
summary(M.0, presid = FALSE)


# Simulated data set ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
set.seed(123); nTimePts <- 5; N <- 1000  # N is the popn size
pdata <- rposbern(n = N, nTimePts = nTimePts, pvars = 2, is.popn = TRUE)
nrow(pdata)  # Less than N (because some animals were never captured)
# The truth: xcoeffs are c(-2, 1, 2) and cap.effect = +1

M.bh.2 <- vglm(cbind(y1, y2, y3, y4, y5) ~ x2,
               posbernoulli.b, data = pdata, trace = TRUE)
coef(M.bh.2)
coef(M.bh.2, matrix = TRUE)
constraints(M.bh.2, matrix = TRUE)
summary(M.bh.2, presid = FALSE)
head(depvar(M.bh.2))    # Capture history response matrix
head(M.bh.2@extra$cap.hist1)  # Info on its capture history
head(M.bh.2@extra$cap1)  # When it was first captured
head(fitted(M.bh.2))     # Depends on capture history
(trap.effect <- coef(M.bh.2)["(Intercept):1"])  # Should be +1
head(model.matrix(M.bh.2, type = "vlm"), 21)
head(pdata)
summary(pdata)
dim(depvar(M.bh.2))
vcov(M.bh.2)

M.bh.2@extra$N.hat     # Estimate of the population size; should be about N
M.bh.2@extra$SE.N.hat  # SE of the estimate of the population size
# An approximate 95 percent confidence interval:
round(M.bh.2@extra$N.hat + c(-1, 1) * 1.96 *  M.bh.2@extra$SE.N.hat, 1)
}
\keyword{models}
\keyword{regression}

%# Compare the models using a LRT
%lrtest(M.bh, M.h)
%(wald.pvalue <- 2 * pnorm(abs(summary(M.bh)@coef3["(Intercept):2", "z value"]),
%                          lower.tail = FALSE))  # Two-sided pvalue



