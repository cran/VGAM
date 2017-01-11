\name{weibull.mean}
\alias{weibull.mean}
%\alias{weibullff}
%\alias{weibull.lsh}
%\alias{weibull3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Weibull Distribution Family Function, Parameterized by the Mean }
\description{
  Maximum likelihood estimation of the 2-parameter Weibull distribution.
  The mean is one of the  parameters.
  No observations should be censored.

}
\usage{
weibull.mean(lmean = "loge", lshape = "loge", imean = NULL,
             ishape = NULL, probs.y = c(0.2, 0.5, 0.8),
             imethod = 1, zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmean, lshape}{
  Parameter link functions applied to the
  (positive) mean parameter (called \eqn{mu} below) and
  (positive) shape parameter (called \eqn{a} below).
  See \code{\link{Links}} for more choices.


  }
  \item{imean, ishape}{
  Optional initial values for the mean and shape parameters.


  }
  \item{imethod, zero, probs.y}{
  Details at \code{\link{CommonVGAMffArguments}}.

  }
}
\details{
  See \code{\link{weibullR}} for most of the details
  for this family function too.
  The mean of \eqn{Y} is \eqn{b \, \Gamma(1+ 1/a)}{b * gamma(1+ 1/a)}
  (returned as the fitted values),
  and this is the first parameter (a \code{\link{loge}}
  link is the default because it is positive).
  The other parameter is the positive shape paramter \eqn{a},
  also having a default \code{\link{loge}} link.


  This \pkg{VGAM} family function currently does not handle
  censored data.
  Fisher scoring is used to estimate the two parameters.
  Although the expected information matrices used here are valid
  in all regions of the parameter space,
  the regularity conditions for maximum
  likelihood estimation are satisfied only if \eqn{a>2}
  (according to Kleiber and Kotz (2003)).
  If this is violated then a warning message is issued.
  One can enforce \eqn{a>2} by
  choosing \code{lshape = logoff(offset = -2)}.
  Common values of the shape parameter lie between 0.5 and 3.5.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.



}
\author{ T. W. Yee }
\note{
  See \code{\link{weibullR}} for more details.
  This \pkg{VGAM} family function handles multiple responses.


}
%\section{Warning}{
% This function is under development to handle other censoring situations.
% The version of this function which will handle censored data will be


%}

\seealso{
    \code{\link{weibullR}},
    \code{\link[stats:Weibull]{dweibull}},
    \code{\link{truncweibull}},
    \code{\link{gev}},
    \code{\link{lognormal}},
    \code{\link{expexpff}},
    \code{\link{maxwell}},
    \code{\link{rayleigh}},
    \code{\link{gumbelII}}.


}
\examples{
wdata <- data.frame(x2 = runif(nn <- 1000))  # Complete data
wdata <- transform(wdata, mu     = exp(-1 + 1 * x2),
                          x3     = rnorm(nn),
                          shape1 = exp(1),
                          shape2 = exp(2))
wdata <- transform(wdata,
  y1 = rweibull(nn, shape = shape1, scale = mu / gamma(1 + 1/shape1)),
  y2 = rweibull(nn, shape = shape2, scale = mu / gamma(1 + 1/shape2)))
fit <- vglm(cbind(y1, y2) ~ x2 + x3, weibull.mean, data = wdata, trace = TRUE)
coef(fit, matrix = TRUE)
sqrt(diag(vcov(fit)))  # SEs
summary(fit, presid = FALSE)
}
\keyword{models}
\keyword{regression}
