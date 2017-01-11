\name{gengamma.stacy}
\alias{gengamma.stacy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Gamma distribution family function }
\description{
  Estimation of the 3-parameter generalized gamma distribution proposed by
  Stacy (1962).

}
\usage{
gengamma.stacy(lscale = "loge", ld = "loge", lk = "loge",
               iscale = NULL, id = NULL, ik = NULL,
               gscale    = exp(-5:5), gshape1.d = exp(-5:5),
               gshape2.k = exp(-5:5), zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, ld, lk}{
  Parameter link function applied to each of the positive parameters
  \eqn{b}, \eqn{d} and \eqn{k}, respectively.
  See \code{\link{Links}} for more choices.


  }
  \item{iscale, id, ik}{
  Initial value for \eqn{b}, \eqn{d} and \eqn{k}, respectively.
  The defaults mean an initial value is determined internally for each.


  }
  \item{gscale, gshape1.d, gshape2.k}{
    See \code{\link{CommonVGAMffArguments}} for information.
    Replaced by \code{iscale}, \code{id} etc. if given.


  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


% An integer-valued vector specifying which
% linear/additive predictors are modelled as intercepts only.
% The values must be from the set \{1,2,3\}.
% The default value means none are modelled as intercept-only terms.


  }
}
\details{
  The probability density function can be written
  \deqn{f(y;b,d,k) = d   b^{-d k}   y^{d k-1}  \exp[-(y/b)^d] / \Gamma(k)}{%
        f(y;b,d,k) = d * b^(-d*k) * y^(d*k-1) * exp(-(y/b)^d) /  gamma(k)}
  for scale parameter \eqn{b > 0}, and \eqn{d > 0}, \eqn{k > 0},
  and \eqn{y > 0}.
  The mean of \eqn{Y}
  is \eqn{b \times \Gamma(k+1/d) / \Gamma(k)}{b*gamma(k+1/d)/gamma(k)}
  (returned as the fitted values),
  which equals \eqn{bk}{b*k} if \eqn{d=1}.


There are many special cases, as given in Table 1 of Stacey and Mihram (1965).
In the following, the parameters are in the order \eqn{b,d,k}.
The special cases are:
Exponential       \eqn{f(y;b,1,1)},
Gamma             \eqn{f(y;b,1,k)},
Weibull           \eqn{f(y;b,d,1)},
Chi Squared       \eqn{f(y;2,1,a/2)} with \eqn{a} degrees of freedom,
Chi               \eqn{f(y;\sqrt{2},2,a/2)}{f(y;sqrt(2),2,a/2)} with \eqn{a} degrees of freedom,
Half-normal       \eqn{f(y;\sqrt{2},2,1/2)}{f(y;sqrt(2),2,1/2)},
Circular normal   \eqn{f(y;\sqrt{2},2,1)}{f(y;sqrt(2),2,1)},
Spherical normal  \eqn{f(y;\sqrt{2},2,3/2)}{f(y;sqrt(2),2,3/2)},
Rayleigh          \eqn{f(y;c\sqrt{2},2,1)}{f(y;c sqrt(2),2,1)} where \eqn{c>0}.


}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Stacy, E. W. (1962)
  A generalization of the gamma distribution.
  \emph{Annals of Mathematical Statistics}, \bold{33}(3), 1187--1192.


  Stacy, E. W. and Mihram, G. A. (1965)
  Parameter estimation for a generalized gamma distribution.
  \emph{Technometrics}, \bold{7}, 349--358.


  Prentice, R. L. (1974)
  A log gamma model and its maximum likelihood estimation.
  \emph{Biometrika}, \bold{61}, 539--544.


}
\section{Warning }{
  Several authors have considered maximum likelihood estimation for the
  generalized gamma distribution and have found that the Newton-Raphson
  algorithm does not work very well and that the existence of solutions
  to the log-likelihood equations is sometimes in doubt.
  Although Fisher scoring is used here, it is likely that the same
  problems will be encountered. It appears that large samples
  are required, for example, the estimator of \eqn{k} became asymptotically
  normal only with 400 or more observations.
  It is not uncommon for maximum likelihood estimates to fail to converge
  even with two or three hundred observations.
  With covariates, even more observations are needed to increase the
  chances of convergence.
  Using covariates is not advised unless the sample size is at least
  a few thousand.


}
\author{ T. W. Yee }
\note{
  The notation used here differs from Stacy (1962) and Prentice (1974).
  Poor initial values may result in failure to converge so
  if there are covariates and there are convergence problems,
  try using the \code{zero} argument (e.g., \code{zero = 2:3})
  or the \code{ik} argument.


}
\seealso{
  \code{\link{rgengamma.stacy}},
  \code{\link{gamma1}},
  \code{\link{gamma2}},
  \code{\link{prentice74}},
  \code{\link{simulate.vlm}}.


}
\examples{
k <- exp(-1); Scale = exp(1)
gdata <- data.frame(y = rgamma(1000, shape = k, scale = Scale))
fit <- vglm(y ~ 1, gengamma.stacy, data = gdata, trace = TRUE)
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}

%# Another example
%gdata <- data.frame(x2 = runif(nn <- 5000))
%gdata <- transform(gdata, Scale = exp(1),
%                          d = exp( 0 + 1.2* x2),
%                          k = exp(-1 + 2  * x2))
%gdata <- transform(gdata, y = rgengamma.stacy(nn, scale = Scale, d = d, k = k))
%fit <- vglm(y ~ x2, gengamma.stacy(zero = 1, iscale = 6), data = gdata, trace = TRUE)
%fit <- vglm(y ~ x2, gengamma.stacy(zero = 1), data = gdata, trace = TRUE, maxit = 50)
%coef(fit, matrix = TRUE)
