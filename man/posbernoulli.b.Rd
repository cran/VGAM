\name{posbernoulli.b}
%\alias{posbernoulli}
\alias{posbernoulli.b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Bernoulli Family Function with Behavioural Effects }
\description{
  Fits a GLM-like model to multiple Bernoulli responses where
  each row in the capture history matrix response has at least one success
  (capture).
  Capture history behavioural effects are accommodated.


}
\usage{
posbernoulli.b(link = "logit", parallel.b = FALSE, apply.parint = TRUE,
               icap.prob = NULL, irecap.prob = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, parallel.b, apply.parint, icap.prob, irecap.prob}{

  See \code{\link{CommonVGAMffArguments}} for information about
  these arguments.
  With an intercept-only model
  setting \code{parallel.b = TRUE} results in the \eqn{M_0} model;
  it just deletes the 2nd column of the constraint matrix corresponding
  to the intercept.


% The default value of \code{zero} means that the behavioural
% effect is modelled as the difference between the
% two intercepts.


% That is, it is modelled through the intercept, and a
% negative value of the second linear/additive predictor means trap shy, etc.



  }
}
\details{
  This model
  (commonly known as \eqn{M_b} in the capture--recapture literature)
  operates on a capture history matrix response of 0s and 1s.
  See \code{\link{posbernoulli.t}} for details.


  Each sampling occasion has the same probability and this is modelled here.
  But once an animal is captured, it is marked so that its future
  capture history can be recorded. The effect of the recapture
  probability is modelled through a second linear/additive predictor,
  and this usually differs from the first linear/additive predictor
  by just a different intercept (because \code{parallel.b = TRUE}
  but the parallelism does not apply to the intercept).


  It is well-known that some species of animals are affected by capture,
  e.g., trap-shy or trap-happy. This \pkg{VGAM} family function
  \emph{does} allow the capture history to be modelled via such
  behavioural effects.


  See \code{\link{posbernoulli.t}} for other information,
  e.g., common assumptions.


  The number of linear/additive predictors is \eqn{M = 2},
  and the default links
  are \eqn{(logit \,p_c, logit \,p_r)^T}{(logit p_c, logit p_r)^T}
  where \eqn{p_c} is the probability of capture and
        \eqn{p_r} is the probability of recapture.
  The fitted value returned is of the same dimension as
  the response matrix, and depends on the capture history:
  prior to being first captured, it is \code{cap.prob}.
  Afterwards, it is \code{recap.prob}.


  By default, the constraint matrix for the intercept term
  is set up so that \eqn{p_r} differs from \eqn{p_c} by a
  simple binary effect. This allows an estimate of the
  trap-happy/trap-shy effect.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}

\section{Warning }{

  See \code{\link{posbernoulli.tb}}.


}

\references{


  See \code{\link{posbernoulli.t}}.


}
\author{ Thomas W. Yee. }

\note{
  When the number of sampling occasions is large
  the computation becomes increasingly slower.
  Monitor convergence by setting \code{trace = TRUE}.


  The dependent variable is \emph{not} scaled to row proportions.
  This is the same as \code{\link{posbernoulli.t}}
  but different from \code{\link{posbinomial}}
  and \code{\link{binomialff}}.


}

\seealso{ 
  \code{\link{posbernoulli.t}} (including estimating \eqn{N}),
  \code{\link{posbernoulli.tb}},
  \code{\link{Perom}},
  \code{\link{dposbern}},
  \code{\link{rposbern}},
  \code{\link{posbinomial}}.
% \code{\link{huggins91}}.
% \code{\link{vglm.control}} for \code{xij},


}

\examples{
# Perom data ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Fit a M_b model
M_b <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1,
            data = Perom, posbernoulli.b, trace = TRUE)
coef(M_b, matrix = TRUE)
constraints(M_b, matrix = TRUE)
summary(M_b)

# Fit a M_bh model
M_bh <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight,
             posbernoulli.b, trace = TRUE, data = Perom)
coef(M_bh, matrix = TRUE)
constraints(M_bh)  # (2,2) element of "(Intercept)" is the behavioural effect
summary(M_bh)  # Estimate of behavioural effect is positive (trap-happy)

# Fit a M_h model
M_h <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight,
            data = Perom,
            posbernoulli.t(parallel.t = TRUE), trace = TRUE)
coef(M_h, matrix = TRUE)
constraints(M_h, matrix = TRUE)
summary(M_h)

# Fit a M_0 model
M_0 <- vglm(cbind(    y1 + y2 + y3 + y4 + y5 + y6,
                  6 - y1 - y2 - y3 - y4 - y5 - y6) ~ 1,
            data = Perom, posbinomial, trace = TRUE)
coef(M_0, matrix = TRUE)
constraints(M_0, matrix = TRUE)
summary(M_0)


# Simulated data set ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
set.seed(123); nTimePts <- 5; N <- 1000
hdata <- rposbern(n = N, nTimePts = nTimePts, pvars = 2,
                  is.popn = TRUE)  # N is the popn size
nrow(hdata)  # Less than N
# The truth: xcoeffs are c(-2, 1, 2) and cap.effect = -1

model1 <- vglm(cbind(y1, y2, y3, y4, y5) ~ x2,
               posbernoulli.b, data = hdata, trace = TRUE)
coef(model1)
coef(model1, matrix = TRUE)
constraints(model1, matrix = TRUE)
summary(model1)
head(depvar(model1))    # Capture history response matrix
head(model1@extra$cap.hist1)  # Info on its capture history
head(model1@extra$cap1)  # When it was first captured
head(fitted(model1))     # Depends on capture history
(trap.effect <- coef(model1)["(Intercept):2"])  # Should be -1
head(model.matrix(model1, type = "vlm"), 21)
head(hdata)
summary(hdata)
dim(depvar(model1))
vcov(model1)

model1@extra$N.hat     # Estimate of the population size; should be about N
model1@extra$SE.N.hat  # SE of the estimate of the population size
# An approximate 95 percent confidence interval:
round(model1@extra$N.hat + c(-1, 1) * 1.96 *  model1@extra$SE.N.hat, 1)
}
\keyword{models}
\keyword{regression}

%# Compare the models using a LRT
%lrtest(M_bh, M_h)
%(wald.pvalue <- 2 * pnorm(abs(summary(M_bh)@coef3["(Intercept):2", "z value"]),
%                          lower.tail = FALSE))  # Two-sided pvalue



