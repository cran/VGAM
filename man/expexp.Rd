\name{expexp}
\alias{expexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exponentiated Exponential Distribution }
\description{
  Estimates the two parameters of the exponentiated exponential
  distribution by maximum likelihood estimation.

}
\usage{
expexp(lshape = "loge", lscale = "loge",
       ishape = 1.1, iscale = NULL, tolerance = 1.0e-6, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, lscale}{
  Parameter link functions for the
  \eqn{\alpha}{shape} and \eqn{\lambda}{scale} parameters.
  See \code{\link{Links}} for more choices.
  The defaults ensure both parameters are positive.

  }
  \item{ishape}{
  Initial value for the \eqn{\alpha}{shape} 
  parameter. If convergence fails try setting a different
  value for this argument.

  }
  \item{iscale}{
  Initial value for the \eqn{\lambda}{scale} parameter.
  By default, an initial value is chosen internally using
  \code{ishape}.

  }
  \item{tolerance}{
  Numeric. Small positive value for testing whether values
  are close enough to 1 and 2.

  }
  \item{zero}{ An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The default is none of them.
  If used, choose one value from the set \{1,2\}.

  }
}
\details{
  The exponentiated exponential distribution is an alternative
  to the Weibull and the gamma distributions.
  The formula for the density is
  \deqn{f(y;\alpha,\lambda) =
    \alpha \lambda (1-\exp(-\lambda y))^{\alpha-1}
    \exp(-\lambda y) }{%
    f(y;shape,scale) =
    shape scale (1-\exp(-scale y))^(shape-1)
    \exp(-scale y) }
  where \eqn{y>0},
  \eqn{\alpha>0}{shape>0} and
  \eqn{\lambda>0}{scale>0}.
  The mean of \eqn{Y} is
  \eqn{(\psi(\alpha+1)-\psi(1))/\lambda}{(psi(shape+1)-psi(1))/scale}
  (returned as the fitted values)
  where \eqn{\psi}{psi} is the digamma function.
  The variance of \eqn{Y} is
  \eqn{(\psi'(1)-\psi'(\alpha+1))/\lambda^2}{(psi'(1)-psi'(shape+1))/
    scale^2}
  where \eqn{\psi'}{psi'} is the trigamma function.


  This distribution has been called the two-parameter generalized
  exponential distribution by Gupta and Kundu (2006).
  A special case of the exponentiated exponential distribution:
  \eqn{\alpha=1}{shape=1} is the exponential distribution.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
  Gupta, R. D. and Kundu, D. (2001)
  Exponentiated exponential family: an alternative to
  gamma and Weibull distributions,
  \emph{Biometrical Journal},
  \bold{43},
  117--130.


  Gupta, R. D. and Kundu, D. (2006)
  On the comparison of Fisher information of the
  Weibull and GE distributions,
  \emph{Journal of Statistical Planning and Inference},
  \bold{136},
  3130--3144.


}
\author{ T. W. Yee }
\note{
  Fisher scoring is used, however, convergence is usually very slow.
  This is a good sign that there is a bug, but I have yet to check
  that the expected information is correct.
  Also, I have yet to implement Type-I right censored data using the
  results of Gupta and Kundu (2006).


  Another algorithm for fitting this model is implemented in
  \code{\link{expexp1}}.


}
\section{Warning }{
  Practical experience shows that reasonably good initial values really
  helps. In particular, try setting different values for the \code{ishape}
  argument if numerical problems are encountered or failure to convergence
  occurs. Even if convergence occurs try perturbing the initial value
  to make sure the global solution is obtained and not a local solution.
  The algorithm may fail if the estimate of the shape parameter is
  too close to unity.


}

\seealso{
  \code{\link{expexp1}},
  \code{\link{gamma2.ab}},
  \code{\link{weibull}},
  \code{\link{CommonVGAMffArguments}}.


}
\examples{
# A special case: exponential data
edata <- data.frame(y = rexp(n <- 1000))
fit <- vglm(y ~ 1, fam = expexp, edata, trace = TRUE, maxit = 99)
coef(fit, matrix = TRUE)
Coef(fit)


# Ball bearings data (number of million revolutions before failure)
bbearings <- c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60,
48.80, 51.84, 51.96, 54.12, 55.56, 67.80, 68.64, 68.64,
68.88, 84.12, 93.12, 98.64, 105.12, 105.84, 127.92,
128.04, 173.40)
fit <- vglm(bbearings ~ 1, fam = expexp(iscale = 0.05, ish = 5),
            trace = TRUE, maxit = 300)
coef(fit, matrix = TRUE)
Coef(fit)   # Authors get c(shape=5.2589, scale=0.0314)
logLik(fit) # Authors get -112.9763


# Failure times of the airconditioning system of an airplane
acplane <- c(23, 261, 87, 7, 120, 14, 62, 47,
225, 71, 246, 21, 42, 20, 5, 12, 120, 11, 3, 14,
71, 11, 14, 11, 16, 90, 1, 16, 52, 95)
fit <- vglm(acplane ~ 1, fam = expexp(ishape = 0.8, isc = 0.15),
            trace = TRUE, maxit = 99)
coef(fit, matrix = TRUE)
Coef(fit)   # Authors get c(shape=0.8130, scale=0.0145)
logLik(fit) # Authors get log-lik -152.264
}
\keyword{models}
\keyword{regression}

