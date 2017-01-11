\name{gamma2}
\alias{gamma2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 2-parameter Gamma Distribution }
\description{
  Estimates the 2-parameter gamma distribution
  by maximum likelihood estimation.

}
\usage{
gamma2(lmu = "loge", lshape = "loge",
       imethod = 1,  ishape = NULL,
       parallel = FALSE, deviance.arg = FALSE, zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
% apply.parint = FALSE,
\arguments{
  \item{lmu, lshape}{
  Link functions applied to the (positive) \emph{mu} and \emph{shape}
  parameters (called \eqn{\mu}{mu} and \eqn{a}{shape} respectively).
  See \code{\link{Links}} for more choices.


  }
  \item{ishape}{
  Optional initial value for \emph{shape}.
  A \code{NULL} means a value is computed internally.
  If a failure to converge occurs, try using this argument.
  This argument is ignored if used within \code{\link{cqo}}; see the
  \code{iShape} argument of \code{\link{qrrvglm.control}} instead.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method for the \eqn{\mu}{mu} parameter.
  If failure to converge occurs
  try another value (and/or specify a value for \code{ishape}).


  }
  \item{deviance.arg}{
  Logical. If \code{TRUE}, the deviance function
  is attached to the object. Under ordinary circumstances, it should
  be left alone because it really assumes the shape parameter is at
  the maximum likelihood estimate. Consequently, one cannot use that
  criterion to minimize within the IRLS algorithm.
  It should be set \code{TRUE} only when used with \code{\link{cqo}}
  under the fast algorithm.


  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


% An integer specifying which
% linear/additive predictor is to be modelled as an intercept only.
% If assigned, the single value should be either 1 or 2 or \code{NULL}.
% The default is to model \eqn{shape} as an intercept only.
% A value \code{NULL} means neither 1 or 2.

%   Integer valued vector, usually assigned \eqn{-2} or \eqn{2} if
%   used at all.  Specifies which of the two linear/additive predictors
%   are modelled as an intercept only. By default, the shape parameter
%   (after \code{lshape} is applied) is modelled as a single unknown
%   number that is estimated.  It can be modelled as a function of
%   the explanatory variables by setting \code{zero = NULL}.  A negative
%   value means that the value is recycled, so setting \eqn{-2} means
%   all shape parameters are intercept only.
%   See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{parallel}{
  Details at \code{\link{CommonVGAMffArguments}}.
  If \code{parallel = TRUE} then the constraint is not applied to the intercept.



  }
}
\details{
  This distribution can model continuous skewed responses.
  The density function is given by
  \deqn{f(y;\mu,a) = \frac{\exp(-a y / \mu) \times
               (a y / \mu)^{a-1}
               \times a}{
               \mu \times \Gamma(a)}}{%
   f(y;mu,shape) = exp(-shape * y / mu) y^(shape-1) shape^(shape) /
          [mu^(shape) * gamma(shape)]}
  for
  \eqn{\mu > 0}{mu > 0},
  \eqn{a > 0}{shape > 0}
  and \eqn{y > 0}.
  Here,
  \eqn{\Gamma(\cdot)}{gamma()} is the gamma
  function, as in \code{\link[base:Special]{gamma}}.
  The mean of \emph{Y} is \eqn{\mu=\mu}{mu=mu} (returned as the fitted
  values) with variance \eqn{\sigma^2 = \mu^2 / a}{sigma^2 =
  mu^2 / shape}.  If \eqn{0<a<1}{0<shape<1} then the density has a
  pole at the origin and decreases monotonically as \eqn{y} increases.
  If \eqn{a=1}{shape=1} then this corresponds to the exponential
  distribution.  If \eqn{a>1}{shape>1} then the density is zero at the
  origin and is unimodal with mode at \eqn{y = \mu - \mu / a}{y =
  mu - mu / shape}; this can be achieved with \code{lshape="loglog"}.


  By default, the two linear/additive predictors are
  \eqn{\eta_1=\log(\mu)}{eta1=log(mu)} and
  \eqn{\eta_2=\log(a)}{eta2=log(shape)}.
  This family function implements Fisher scoring and the working
  weight matrices are diagonal.


  This \pkg{VGAM} family function handles \emph{multivariate} responses,
  so that a matrix can be used as the response. The number of columns is
  the number of species, say, and \code{zero=-2} means that \emph{all}
  species have a shape parameter equalling a (different) intercept only.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
  The parameterization of this \pkg{VGAM} family function is the
  2-parameter gamma distribution described in the monograph


McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
\author{ T. W. Yee }
\note{
  The response must be strictly positive.
  A moment estimator for the shape parameter may be implemented in
  the future.


  If \code{mu} and \code{shape} are vectors, then \code{rgamma(n = n,
  shape = shape, scale = mu/shape)} will generate random gamma variates of this
  parameterization, etc.;
  see \code{\link[stats]{GammaDist}}.


  For \code{\link{cqo}} and \code{\link{cao}}, taking the logarithm
  of the response means (approximately) a \code{\link{gaussianff}} family
  may be used on the transformed data.


}

\seealso{
  \code{\link{gamma1}} for the 1-parameter gamma distribution,
  \code{\link{gammaR}} for another parameterization of
  the 2-parameter gamma distribution that is directly matched
  with \code{\link[stats]{rgamma}},
  \code{\link{bigamma.mckay}} for \emph{a} bivariate gamma distribution,
  \code{\link{expexpff}},
  \code{\link[stats]{GammaDist}},
  \code{\link{golf}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}},
  \code{\link{negloge}}.


}
\examples{
# Essentially a 1-parameter gamma
gdata <- data.frame(y = rgamma(n = 100, shape = exp(1)))
fit1 <- vglm(y ~ 1, gamma1, data = gdata)
fit2 <- vglm(y ~ 1, gamma2, data = gdata, trace = TRUE, crit = "coef")
coef(fit2, matrix = TRUE)
c(Coef(fit2), colMeans(gdata))

# Essentially a 2-parameter gamma
gdata <- data.frame(y = rgamma(n = 500, rate = exp(-1), shape = exp(2)))
fit2 <- vglm(y ~ 1, gamma2, data = gdata, trace = TRUE, crit = "coef")
coef(fit2, matrix = TRUE)
c(Coef(fit2), colMeans(gdata))
summary(fit2)
}
\keyword{models}
\keyword{regression}

