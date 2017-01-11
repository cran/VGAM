\name{betageometric}
\alias{betageometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Beta-geometric Distribution Family Function }
\description{
  Maximum likelihood estimation for the beta-geometric distribution.

}
\usage{
betageometric(lprob = "logit", lshape = "loge",
              iprob = NULL,    ishape = 0.1,
              moreSummation = c(2, 100), tolerance = 1.0e-10, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lprob, lshape}{
  Parameter link functions applied to the
  parameters \eqn{p}{prob} and \eqn{\phi}{phi}
  (called \code{prob} and \code{shape} below).
  The former lies in the unit interval and the latter is positive.
  See \code{\link{Links}} for more choices.


  }
  \item{iprob, ishape}{
  Numeric.
  Initial values for the two parameters.
  A \code{NULL} means a value is computed internally.


  }
  \item{moreSummation}{
  Integer, of length 2.
  When computing the expected information matrix a series summation from
  0 to \code{moreSummation[1]*max(y)+moreSummation[2]} is made, in which the
  upper limit is an approximation to infinity.
  Here, \code{y} is the response.


  }
  \item{tolerance}{
  Positive numeric.
  When all terms are less than this then the series is deemed to have
  converged.


  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  If used, the value must be from the set \{1,2\}.


  }
}
\details{
  A random variable \eqn{Y} has a 2-parameter beta-geometric distribution
  if \eqn{P(Y=y) = p (1-p)^y}{P(Y=y) = prob * (1-prob)^y}
  for \eqn{y=0,1,2,\ldots}{y=0,1,2,...} where
  \eqn{p}{prob} are generated from a standard beta distribution with
  shape parameters \code{shape1} and \code{shape2}.
  The parameterization here is to focus on the parameters
  \eqn{p}{prob} and
  \eqn{\phi = 1/(shape1+shape2)}{phi = 1/(shape1+shape2)},
  where \eqn{\phi}{phi} is \code{shape}.
  The default link functions for these ensure that the appropriate range
  of the parameters is maintained.
  The mean of \eqn{Y} is
  \eqn{E(Y) = shape2 / (shape1-1) = (1-p) / (p-\phi)}{E(Y) =
       shape2 / (shape1-1) = (1-prob) / (prob-phi)}
  if \code{shape1 > 1}, and if so, then this is returned as
  the fitted values.


  The geometric distribution is a special case of the beta-geometric
  distribution with \eqn{\phi=0}{phi=0} (see \code{\link{geometric}}).
  However, fitting data from a geometric distribution may result in
  numerical problems because the estimate of \eqn{\log(\phi)}{log(phi)}
  will 'converge' to \code{-Inf}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Paul, S. R. (2005)
  Testing goodness of fit of the geometric distribution:
  an application to human fecundability data.
  \emph{Journal of Modern Applied Statistical Methods}, \bold{4}, 425--433.


}

\author{ T. W. Yee }
\note{
  The first iteration may be very slow;
  if practical, it is best for the \code{weights} argument of
  \code{\link{vglm}} etc. to be used rather than inputting a very
  long vector as the response, i.e., \code{vglm(y ~ 1, ..., weights = wts)}
  is to be preferred over \code{vglm(rep(y, wts) ~ 1, ...)}.
  If convergence problems occur try inputting some values of argument
  \code{ishape}.


  If an intercept-only model is fitted then the \code{misc} slot of the
  fitted object has list components \code{shape1} and \code{shape2}.


}

\seealso{
  \code{\link{geometric}},
  \code{\link{betaff}},
  \code{\link{rbetageom}}.


}
\examples{
bdata <- data.frame(y = 0:11, wts = c(227,123,72,42,21,31,11,14,6,4,7,28))
fitb <- vglm(y ~ 1, betageometric, data = bdata, weight = wts, trace = TRUE)
fitg <- vglm(y ~ 1,     geometric, data = bdata, weight = wts, trace = TRUE)
coef(fitb, matrix = TRUE)
Coef(fitb)
sqrt(diag(vcov(fitb, untransform = TRUE)))
fitb@misc$shape1
fitb@misc$shape2
# Very strong evidence of a beta-geometric:
pchisq(2 * (logLik(fitb) - logLik(fitg)), df = 1, lower.tail = FALSE)
}
\keyword{models}
\keyword{regression}

