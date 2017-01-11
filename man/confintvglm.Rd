\name{confintvglm}
%\name{confint}
% \alias{confint}
\alias{confintvglm}
\alias{confintrrvglm}
\alias{confintvgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence Intervals for Parameters of VGLMs }
\description{
  Computes confidence intervals (CIs)
  for one or more parameters in a fitted model.
  Currently the object must be a
  \code{"\link{vglm}"} object.


}
% confint(object, parm, level = 0.95, \dots)
\usage{
confintvglm(object, parm, level = 0.95, method = c("wald", "profile"),
            trace = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A fitted model object.


  }
  \item{parm, level, \dots}{Same as \code{\link[stats]{confint}}.
  }
  \item{method}{Character.
  The default is the first method.
  Abbreviations are allowed.
  Currently \code{"profile"} is basically working;
  and it is likely to be more accurate especially for
  small samples, as it is based on a profile log likelihood,
  however it is computationally intensive.


  }
  \item{trace}{
    Logical. If \code{TRUE} then one can monitor the
    computation as it progresses (because it is expensive).
    The default is the orginal model's \code{trace} value
    (see \code{\link{vglm.control}}).
    Setting \code{FALSE} suppresses all intermediate output.


  }
}
\details{
  The default for
  this methods function is based on \code{\link[stats]{confint.default}}
  and assumes
  asymptotic normality. In particular,
  the \code{\link[VGAM:coefvlm]{coef}} and
  \code{vcov} methods functions are used for
  \code{\link[VGAM]{vglm-class}} objects.



  When \code{method = "profile"} the function \code{profilevglm()}
  is called to do the profiling. The code is very heavily
  based on \code{\link[MASS]{profile.glm}}
  which was originally written by
  D. M. Bates and W. N. Venables (For S in 1996)
  and subsequently corrected by B. D. Ripley.
  Sometimes the profiling method can give problems, for
  example, \code{\link{cumulative}} requires the \eqn{M}
  linear predictors not to intersect in the data cloud.
  Such numerical problems are less common when
  \code{method = "wald"}, however, it is well-known
  that inference based on profile likelihoods is generally
  more accurate than Wald, especially when the sample size
  is small.
  The deviance (\code{deviance(object)}) is used if possible,
  else the difference
  \code{2 * (logLik(object) - ell)} is computed,
  where \code{ell} are the values of the loglikelihood on a grid.



  For
  Wald CIs and
  \code{\link[VGAM]{rrvglm-class}}
  objects, currently an error message is produced because
  I haven't gotten around to write the methods function;
  it's not too hard, but am too busy!
  An interim measure is to
  coerce the object into a \code{"\link{vglm}"} object,
  but then the confidence intervals will tend to be too narrow because
  the estimated constraint matrices are treated as known.



  For
  Wald CIs and
  \code{\link[VGAM]{vgam-class}}
  objects, currently an error message is produced because
  the theory is undeveloped.



}
\value{
  Same as \code{\link[stats]{confint}}.



}
%\references{
%}
\author{
  Thomas Yee adapted \code{\link[stats]{confint.lm}}
  to handle \code{"vglm"} objects, for Wald-type
  confidence intervals.
  Also, \code{\link[MASS]{profile.glm}}
  was originally written by
  D. M. Bates and W. N. Venables (For S in 1996)
  and subsequently corrected by B. D. Ripley.
  This function effectively calls \code{confint.profile.glm()}
  in \pkg{MASS}.



  }

\note{
  The order of the values of argument \code{method} may change
  in the future without notice.
  The functions
  \code{plot.profile.glm}
  and
  \code{pairs.profile.glm}
  from \pkg{MASS}
  appear to work with output from this function.


}
%\section{Warning }{
%}


\seealso{
  \code{\link{vcovvlm}},
  \code{\link{summaryvglm}},
  \code{\link[stats]{confint}},
  \code{\link[MASS]{profile.glm}},
  \code{plot.profile.glm},
  \code{pairs.profile.glm}.



}

\examples{
# Example 1: this is based on a glm example
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3, 1, 9); treatment <- gl(3, 3)
 glm.D93 <-  glm(counts ~ outcome + treatment, family = poisson())
vglm.D93 <- vglm(counts ~ outcome + treatment, family = poissonff)
confint(glm.D93) # needs MASS to be present on the system
confint.default(glm.D93)  # based on asymptotic normality
confint(vglm.D93)
confint(vglm.D93) - confint(glm.D93)    # Should be all 0s
confint(vglm.D93) - confint.default(glm.D93)  # based on asympt. normality

# Example 2: simulated negative binomial data with multiple responses
ndata <- data.frame(x2 = runif(nn <- 100))
ndata <- transform(ndata, y1 = rnbinom(nn, mu = exp(3+x2), size = exp(1)),
                          y2 = rnbinom(nn, mu = exp(2-x2), size = exp(0)))
fit1 <- vglm(cbind(y1, y2) ~ x2, negbinomial, data = ndata, trace = TRUE)
coef(fit1)
coef(fit1, matrix = TRUE)
confint(fit1)
confint(fit1, "x2:1")  #  This might be improved to "x2" some day...
\dontrun{
confint(fit1, method = "profile")  # Computationally expensive
confint(fit1, "x2:1", method = "profile", trace = FALSE)
}

fit2 <- rrvglm(y1 ~ x2, negbinomial(zero = NULL), data = ndata)
confint(as(fit2, "vglm"))  # Too narrow (SEs are biased downwards)
}
\keyword{models}
\keyword{regression}

