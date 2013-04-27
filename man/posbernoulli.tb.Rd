\name{posbernoulli.tb}
%\alias{posbernoulli}
\alias{posbernoulli.tb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Bernoulli Family Function with Time and Behavioural Effects (experimental) }
\description{
  Fits a GLM-like model to multiple
  (currently only two or three)
  Bernoulli responses where
  each row in the capture history matrix response has at least one success
  (capture).
  Sampling occasion effects and behavioural effects are accommodated.
  However, this function only handles two and three sampling occasions.


}
\usage{
posbernoulli.tb(link = "logit", parallel.t = FALSE, parallel.b = FALSE,
                apply.parint = FALSE, imethod = 1, iprob = NULL,
                dconst = 0.1, dpower = -2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, imethod, iprob, parallel.t, parallel.b, apply.parint}{
  See \code{\link{CommonVGAMffArguments}} for information.
  But \code{parallel.t} and \code{parallel.b} must each be
  logicals only.
  Argument \code{parallel.t} means parallel with respect to time, and
  matches the same argument name in \code{\link{posbernoulli.t}}.



  Suppose the model is intercept-only.
  Setting \code{parallel.t = TRUE} results in the \eqn{M_b} model.
  Setting \code{parallel.b = TRUE} results in the \eqn{M_t} model.
  Setting \code{parallel.t = TRUE} and
  setting \code{parallel.b = TRUE} results in the \eqn{M_0} model.
  Note the default for \code{parallel.t} and \code{parallel.b}
  (both \code{FALSE})
  may be unsuitable for most data sets which have a large \eqn{\tau}
  because of the large number of parameters; it can be too flexible.
  Note that adding covariates will result in a \eqn{M_{tbh}} model.


  }
  \item{dconst, dpower}{
  Decay constants and power (exponent) for the ridge adjustment
  for the working weight matrices.
  At iteration \eqn{t} of the IRLS algorithm
  a positive value is added to the first \eqn{\tau}{tau}
  diagonal elements of the working weight matrices to make
  them positive-definite. This adjustment is \eqn{K \times t^p}{K * t^p}
  where \eqn{K} is \code{dconst} and \eqn{p} is \code{dpower}.
  This is always positive but decays to zero as iterations proceed
  (provided \eqn{p} is negative etc.).


  }
}
\details{
  This model
  (commonly known as \eqn{M_{tb}} in the capture--recapture literature)
  operates on a response matrix of 0s and 1s.
  See \code{\link{posbernoulli.t}}
  for information that is in common.


  This \pkg{VGAM} family function is \emph{experimental only}.
  When finished, it should allow time and behavioural effects to be modelled.
  Evidently, the expected information matrix (EIM) is \emph{not} of
  full rank, so \code{dconst} and \code{dpower} are used to
  \emph{try} fix up the problem.
  The default link functions
  are \eqn{(logit \,p_{c1},\ldots,logit \,p_{c\tau},logit \,p_{r2},\ldots,logit \,p_{r\tau})^T}{
           (logit p_{c1},\ldots,logit p_{c,tau},logit p_{r2},\ldots,logit p_{r,tau})^T}
  where the subscript \eqn{c} denotes capture,
        the subscript \eqn{r} denotes recapture,
  and it is not possible to recapture the animal at sampling occasion 1.
  Thus \eqn{M = 2\tau - 1}{M=2*tau-1}.
  The parameters are currently prefixed by \code{cap.prob} and \code{recap.prob}
  for the capture and recapture probabilities.


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
  It is a good idea to constrain the probabilities for each sampling
  occasion to be equal, and also allow the behavioural effect to
  be modelled using the intercept.
  See \code{M_tbh.1} below.


  The current restriction of handling only \eqn{\tau=2}{tau=2} and
  \eqn{\tau=3}{tau=3} sampling occasions
  is unfortunate and more work is needed to extend this to four or more.


Data-wise, at each sampling occasion, the \eqn{M_{tb}} model requires at least
one first capture and at least one noncapture.

% If not all of the \eqn{2^{\tau}-1}{2^(tau) - 1} combinations of
% the response matrix are not present then it pays to add
% such rows to the response matrix and assign a small but
% positive prior weight.
% For example, if \eqn{\tau=2}{tau=2} then there should be
% (0,1) rows,
% (1,0) rows and
% (1,1) rows present in the response matrix.


}

\section{Warning }{

  As this model is likely to be overparameterized, probably this
  function should not be used (for now?).


% From Jakub:
  Estimation for the population size (and its SE) for the
  \eqn{M_{tb}} model may be wrong.
  Models \eqn{M_{tbh}} and \eqn{M_{th}} may be wrong.
  But models \eqn{M_{bh}}, \eqn{M_{h}}, \eqn{M_{b}},
  \eqn{M_{t}}, \eqn{M_{0}} seem fine.



  Inference, especially using standard errors, may be fraught here
  because the EIM is, strictly speaking, not of full rank.
  A similar adjustment is made by \code{\link{zipebcom}}.
  It is a good idea to monitor convergence.
  The \eqn{M_0} model is best fitted with \code{\link{posbernoulli.b}}
  or \code{\link{posbernoulli.t}} or \code{\link{posbinomial}} because
  the standard errors are more accurate.


}

\seealso{ 
  \code{\link{posbernoulli.b}} (including \eqn{\widehat{N}}),
  \code{\link{posbernoulli.t}},
  \code{\link{posbinomial}}.


}

\examples{
\dontrun{
# Example 1: simulated data
set.seed(123)
nTimePts <- 2  # Must be 2 or 3 currently (aka tau == # of sampling occasions)
nnn <- 10000   # Number of animals
pdata <- rposbern(n = nnn, nTimePts = nTimePts, pvars = 2)
dim(pdata)
head(pdata)

clist <- list("(Intercept)" = cbind(1, c(0, 0, 1)),  # Capture effect is last coln
              x2            = rbind(1, 1, 1))
M_tbh.1 <- vglm(cbind(y1, y2) ~ x2,
                constraints = clist, trace = TRUE,
                posbernoulli.tb, data = pdata)
summary(M_tbh.1)

coef(M_tbh.1)
coef(M_tbh.1, matrix = TRUE)
constraints(M_tbh.1, matrix = TRUE)
summary(M_tbh.1)  # Standard errors are very approximate
head(fitted(M_tbh.1))
head(model.matrix(M_tbh.1, type = "vlm"), 21)
dim(depvar(M_tbh.1))


# Example 2: Perom subset data
Hlist <- list("(Intercept)" = cbind(1, c(0, 0, 0, 1, 1)),
              sex           = rbind(1, 1, 1, 1, 1),
              weight        = rbind(1, 1, 1, 1, 1))
Psubset <- subset(Perom, y1 + y2 + y3 > 0)
head(Psubset)

fit1 <- vglm(cbind(y1, y2, y3) ~ sex + weight, constraints = Hlist,
             posbernoulli.tb, data = Psubset, trace = TRUE)
coef(fit1)
coef(fit1, matrix = TRUE)
summary(fit1)  # Standard errors are very approximate

# fit1 is the same as Fit1:
Fit1 <- vglm(cbind(y1, y2, y3) ~ sex + weight, data = Psubset,
             posbernoulli.tb(parallel.t = TRUE), trace = TRUE)
constraints(Fit1)  # Same as Hlist

yyy <- depvar(fit1)
if (length(table(4 * yyy[, 1] + 2 * yyy[, 2] + 1 * yyy[, 3])) != 2^(ncol(yyy))-1)
  warning("not every combination is represented by a row in the response matrix")
}
}
\keyword{models}
\keyword{regression}

