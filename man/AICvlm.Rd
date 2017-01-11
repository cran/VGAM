\name{AICvlm}
\alias{AICvlm}
%\alias{AICvglm}
\alias{AICvgam}
\alias{AICrrvglm}
\alias{AICqrrvglm}
\alias{AICrrvgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Akaike's Information Criterion }
\description{
  Calculates the Akaike information criterion for a fitted model object
  for which a log-likelihood value has been obtained.

}
\usage{
    AICvlm(object, \dots, corrected = FALSE, k = 2)
   AICvgam(object, \dots, k = 2)
 AICrrvglm(object, \dots, k = 2)
AICqrrvglm(object, \dots, k = 2)
 AICrrvgam(object, \dots, k = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Some \pkg{VGAM} object, for example, having
  class \code{\link{vglmff-class}}.


  }
  \item{\dots}{
  Other possible arguments fed into
  \code{logLik} in order to compute the log-likelihood.


  }
  \item{corrected}{
  Logical, perform the finite sample correction?


  }
  \item{k}{
  Numeric, the penalty per parameter to be used;
  the default is the classical AIC.


  }
}
\details{
  The following formula is used for VGLMs:
  \eqn{-2 \mbox{log-likelihood} + k n_{par}}{-2*log-likelihood
  + k*npar}, where \eqn{n_{par}}{npar} represents the number of parameters
  in the fitted model, and \eqn{k = 2} for the usual AIC.
  One could assign \eqn{k = \log(n)} (\eqn{n} the number of observations)
  for the so-called BIC or SBC (Schwarz's Bayesian criterion).
  This is the function \code{AICvlm()}.


  This code relies on the log-likelihood being defined, and computed,
  for the object.
  When comparing fitted objects, the smaller the AIC, the better the fit.
  The log-likelihood and hence the AIC is only defined up to an additive
  constant.


  Any estimated scale parameter (in GLM parlance) is used as one
  parameter.


  For VGAMs and CAO the nonlinear effective degrees of freedom for each
  smoothed component is used. This formula is heuristic.
  These are the functions \code{AICvgam()} and \code{AICcao()}.


  The finite sample correction is usually recommended when the
  sample size is small or when the number of parameters is large.
  When the sample size is large their difference tends to be negligible.
  The correction is described in Hurvich and Tsai (1989), and is based
  on a (univariate) linear model with normally distributed errors.



}
\value{
  Returns a numeric value with the corresponding AIC (or BIC, or \dots,
  depending on \code{k}).


}
\author{T. W. Yee. }
\note{
  AIC has not been defined for QRR-VGLMs yet.


  Using AIC to compare \code{\link{posbinomial}} models
  with, e.g., \code{\link{posbernoulli.tb}} models,
  requires \code{posbinomial(omit.constant = TRUE)}.
  See \code{\link{posbinomial}} for an example.
  A warning is given if it suspects a wrong \code{omit.constant} value
  was used.



  Where defined,
  \code{AICc(...)} is the same as \code{AIC(..., corrected = TRUE)}.



}

\references{

Hurvich, C. M. and Tsai, C.-L. (1989)
Regression and time series model selection in small samples,
\emph{Biometrika},
\bold{76}, 297--307.



%  Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986).
%  \emph{Akaike Information Criterion Statistics}.
%  D. Reidel Publishing Company.
}

\section{Warning }{
  This code has not been double-checked.
  The general applicability of \code{AIC} for the VGLM/VGAM classes
  has not been developed fully.
  In particular, \code{AIC} should not be run on some \pkg{VGAM} family
  functions because of violation of certain regularity conditions, etc.


}

\seealso{
  VGLMs are described in \code{\link{vglm-class}};
  VGAMs are described in \code{\link{vgam-class}};
  RR-VGLMs are described in \code{\link{rrvglm-class}};
  \code{\link[stats]{AIC}},
  \code{\link{BICvlm}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit1 <- vglm(cbind(normal, mild, severe) ~ let,
              cumulative(parallel = TRUE, reverse = TRUE), data = pneumo))
coef(fit1, matrix = TRUE)
AIC(fit1)
AICc(fit1)  # Quick way
AIC(fit1, corrected = TRUE)  # Slow way
(fit2 <- vglm(cbind(normal, mild, severe) ~ let,
              cumulative(parallel = FALSE, reverse = TRUE), data = pneumo))
coef(fit2, matrix = TRUE)
AIC(fit2)
AICc(fit2)
AIC(fit2, corrected = TRUE)
}
\keyword{models}
\keyword{regression}

