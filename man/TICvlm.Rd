\name{TIC}
\alias{TIC}
\alias{TICvlm}
%\alias{TICvglm}
%\alias{TICvgam}
%\alias{TICrrvglm}
%\alias{TICqrrvglm}
%\alias{TICrrvgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Takeuchi's Information Criterion }
\description{
  Calculates the Takeuchi information criterion
  for a fitted model object
  for which a log-likelihood value has been obtained.

}
\usage{
    TIC(object, \dots)
    TICvlm(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \pkg{VGAM} object having
    class \code{\link{vglm-class}}.


%   , for example,

  
  }
  \item{\dots}{
  Other possible arguments fed into
  \code{logLik} in order to compute the log-likelihood.


  }
%  \item{corrected}{
%  Logical, perform the finite sample correction?


%  }
%  \item{k}{
%  Numeric, the penalty per parameter to be used;
%  the default is the classical TIC.


%  }
}
\details{
  The following formula is used for VGLMs:
  \eqn{-2 \mbox{log-likelihood} + 2 trace(V K)}{-2*log-likelihood
    + 2 * trace(V K)},
  where \eqn{V} is the inverse of the EIM from the fitted model,
  and \eqn{K} is the outer product of the score vectors.
  Both \eqn{V} and \eqn{K} are order-\eqn{p.VLM} matrices.
  One has \eqn{V} equal to \code{vcov(object)},
  and \eqn{K} is computed by taking the outer product of
  the output from the \code{deriv} slot multiplied by the
  large VLM matrix and then taking their sum.
  Hence for the huge majority of models,
  the penalty is computed at the MLE and is empirical in nature.
  Theoretically, if the fitted model is the true model then
  AIC equals TIC.

  

  When there are prior weights the score vectors are divided
  by the square root of these,
  because \eqn{ (a_i U_i / \sqrt{a_i})^2 = a_i U_i^2}.



%  This is the function \code{TICvlm()}.


  This code relies on the log-likelihood being defined, and computed,
  for the object.
  When comparing fitted objects, the smaller the TIC, the better the fit.
  The log-likelihood and hence the TIC is only defined up to an additive
  constant.


  Currently
  any estimated scale parameter (in GLM parlance) is ignored by
  treating its value as unity.
  Also,
  currently
  this function is written only for \code{\link{vglm}} objects and
  not \code{\link{vgam}} or \code{\link{rrvglm}}, etc., objects.


}
\value{
  Returns a numeric TIC value.


}
\author{T. W. Yee. }
\note{
  TIC has not been defined for RR-VGLMs, QRR-VGLMs, etc., yet.


  See \code{\link{AICvlm}} about models
  such as \code{\link{posbernoulli.tb}}
  that require \code{posbinomial(omit.constant = TRUE)}.



}

\references{

Takeuchi, K. (1976)
Distribution of informational statistics and a criterion of model
fitting.  (In Japanese).
\emph{Suri-Kagaku} (Mathematic Sciences),
\bold{153}, 12--18.

%Distribution of informational statistics and a criterion of model
%fitting.
%Suri-Kagaku (Mathematic Sciences) 153, 12--18. (In Japanese).


Burnham, K. P. and Anderson, D. R. (2002)
\emph{Model Selection and Multi-Model Inference: A Practical
      Information-Theoretic Approach},
2nd ed. New York, USA: Springer.


}

\section{Warning }{
  This code has not been double-checked.
  The general applicability of \code{TIC} for the VGLM/VGAM classes
  has not been developed fully.
  In particular, \code{TIC} should not be run on some \pkg{VGAM} family
  functions because of violation of certain regularity conditions, etc.



  Some authors note that quite large sample sizes are needed
  for this IC to work reasonably well.

% Sociological Methods and Research article, p.270.


% Some authors note that numerical instability may occur for this IC.
  

}

\seealso{
  VGLMs are described in \code{\link{vglm-class}};
  \code{\link[stats]{AIC}},
  \code{\link{AICvlm}}.
  \code{\link{BICvlm}}.


  
%  VGAMs are described in \code{\link{vgam-class}};
%  RR-VGLMs are described in \code{\link{rrvglm-class}};

  
}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit1 <- vglm(cbind(normal, mild, severe) ~ let,
              cumulative(parallel = TRUE, reverse = TRUE), data = pneumo))
coef(fit1, matrix = TRUE)
TIC(fit1)
(fit2 <- vglm(cbind(normal, mild, severe) ~ let,
              cumulative(parallel = FALSE, reverse = TRUE), data = pneumo))
coef(fit2, matrix = TRUE)
TIC(fit2)
}
\keyword{models}
\keyword{regression}




%uiowa.edu 2011 thesis by Cristina Laura Acion:

%Shibata (1989) noted that the error incurred by this additional
%estimation can cause instability of the model selection results yielded
%by TIC. Therefore, TIC is not universally recommended (Burnham and
%Anderson, 2002).

%However, a data-dependent estimator might also be highly variable. This
%issue discourages some authors to recommend the use of TIC (Burnham and
%Anderson, 2002).
