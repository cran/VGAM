\name{AICvlm}
\alias{AICvlm}
%\alias{AICvglm}
\alias{AICvgam}
\alias{AICrrvglm}
%\alias{AICqrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Akaike's An Information Criterion }
\description{
  Calculates the Akaike information criterion for a fitted model object
  for which a log-likelihood value has been obtained.

}
\usage{
AICvlm(object, \dots, k = 2)
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
  One can assign \eqn{k = \log(n)} (\eqn{n} the number of observations)
  for the so-called BIC or SBC (Schwarz's Bayesian criterion).

  This code relies on the log-likelihood being defined, and computed,
  for the object.
  When comparing fitted objects, the smaller the AIC, the better the fit.
  The log-likelihood and hence the AIC is only defined up to an additive
  constant.

  Any estimated scale parameter (in GLM parlance) is used as one
  parameter.

  For VGAMs the nonlinear effective degrees of freedom for each
  smoothed component is used. This is heuristic.

}
\value{
  Returns a numeric value with the corresponding AIC (or BIC, or \dots,
  depending on \code{k}).

}
\author{T. W. Yee }
\note{
  AIC has not been defined for QRR-VGLMs yet.

}

%\references{
%  Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986).
%  \emph{Akaike Information Criterion Statistics}.
%  D. Reidel Publishing Company.
%}

\section{Warning }{
  This code has not been double-checked.
  The general applicability of \code{AIC} for the VGLM/VGAM classes
  has not been developed fully. 
  In particular, \code{AIC} should not be run on some \pkg{VGAM} family
  functions because of violation of certain regularity conditions, etc.

}

\seealso{
  VGLMs are described in \code{\link{vglm-class}};
  RR-VGLMs are described in \code{\link{rrvglm-class}};
  \code{\link[stats]{AIC}}.

}
\examples{
# Fit a proportional odds model
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
(fit1 = vglm(cbind(normal, mild, severe) ~ let,
           cumulative(parallel=TRUE, reverse=TRUE), pneumo))
coef(fit1, matrix=TRUE)
AIC(fit1)
(fit2 = vglm(cbind(normal, mild, severe) ~ let,
           cumulative(parallel=FALSE, reverse=TRUE), pneumo))
coef(fit2, matrix=TRUE)
AIC(fit2)
}
\keyword{models}
\keyword{regression}

