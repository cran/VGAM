\name{BICvlm}
\alias{BICvlm}
%\alias{BICvglm}
\alias{BICvgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian Information Criterion }
\description{
  Calculates the Bayesian information criterion (BIC) for
  a fitted model object for which a log-likelihood value
  has been obtained.


}
\usage{
BICvlm(object, \dots, k = log(nobs(object)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, \dots}{
  Same as \code{\link{AICvlm}}.


  }
  \item{k}{
  Numeric, the penalty per parameter to be used;
  the default is \code{log(n)} where
  \code{n} is the number of observations).


  }
}
\details{
  The so-called BIC or SBC (Schwarz's Bayesian criterion)
  can be computed by calling \code{\link{AICvlm}} with a
  different \code{k} argument.
  See \code{\link{AICvlm}} for information and caveats.


}
\value{
  Returns a numeric value with the corresponding BIC, or \dots,
  depending on \code{k}.


}
\author{T. W. Yee. }
\note{
  BIC, AIC and other ICs can have have many additive
  constants added to them. The important thing are the
  differences since the minimum value corresponds to the best model.
  Preliminary testing shows absolute differences
  with some \pkg{VGAM} family functions such as
  \code{\link{gaussianff}},
  however, they should agree with non-normal families.


  BIC has not been defined for QRR-VGLMs yet.


}

%\references{
%  Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986).
%  \emph{Akaike Information Criterion Statistics}.
%  D. Reidel Publishing Company.
%}


\section{Warning }{
  Like \code{\link{AICvlm}}, this code has not been double-checked.
  The general applicability of \code{BIC} for the VGLM/VGAM classes
  has not been developed fully.
  In particular, \code{BIC} should not be run on some \pkg{VGAM} family
  functions because of violation of certain regularity conditions, etc.


  Many \pkg{VGAM} family functions such as
  \code{\link{cumulative}} can have the number of
  observations absorbed into the prior weights argument
  (e.g., \code{weights} in \code{\link{vglm}}), either
  before or after fitting.  Almost all \pkg{VGAM} family
  functions can have the number of observations defined by
  the \code{weights} argument, e.g., as an observed frequency.
  \code{BIC} simply uses the number of rows of the model matrix, say,
  as defining \code{n}, hence the user must be very careful
  of this possible error.
  Use at your own risk!!


}

\seealso{
  \code{\link{AICvlm}},
  VGLMs are described in \code{\link{vglm-class}};
  VGAMs are described in \code{\link{vgam-class}};
  RR-VGLMs are described in \code{\link{rrvglm-class}};
  \code{\link[stats]{BIC}},
  \code{\link[stats]{AIC}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit1 <- vglm(cbind(normal, mild, severe) ~ let,
              cumulative(parallel = TRUE, reverse = TRUE), data = pneumo))
coef(fit1, matrix = TRUE)
BIC(fit1)
(fit2 <- vglm(cbind(normal, mild, severe) ~ let,
              cumulative(parallel = FALSE, reverse = TRUE), data = pneumo))
coef(fit2, matrix = TRUE)
BIC(fit2)

# These do not agree in absolute terms:
gdata <- data.frame(x2 = sort(runif(n <- 40)))
gdata <- transform(gdata, y1 = 1 + 2*x2 + rnorm(n, sd = 0.1))
fit.v <- vglm(y1 ~ x2, gaussianff, data = gdata)
fit.g <-  glm(y1 ~ x2, gaussian  , data = gdata)
fit.l <-   lm(y1 ~ x2, data = gdata)
c(BIC(fit.l), BIC(fit.g), BIC(fit.v))
c(AIC(fit.l), AIC(fit.g), AIC(fit.v))
c(AIC(fit.l) - AIC(fit.v),
  AIC(fit.g) - AIC(fit.v))
c(logLik(fit.l), logLik(fit.g), logLik(fit.v))
}
\keyword{models}
\keyword{regression}

