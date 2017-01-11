\name{weibullR}
\alias{weibullR}
%\alias{weibullff}
%\alias{weibull.lsh}
%\alias{weibull3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Weibull Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter Weibull distribution.
  No observations should be censored.

}
\usage{
weibullR(lscale = "loge", lshape = "loge",
         iscale = NULL,   ishape = NULL, lss = TRUE, nrfs = 1,
         probs.y = c(0.2, 0.5, 0.8), imethod = 1, zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, lscale}{
  Parameter link functions applied to the
  (positive) shape parameter (called \eqn{a} below) and
  (positive) scale parameter (called \eqn{b} below).
  See \code{\link{Links}} for more choices.


  }
  \item{ishape, iscale}{
  Optional initial values for the shape and scale parameters.


  }
  \item{nrfs}{
  Currently this argument is ignored.
  Numeric, of length one, with value in \eqn{[0,1]}.
  Weighting factor between Newton-Raphson and Fisher scoring.
  The value 0 means pure Newton-Raphson, while 1 means pure Fisher scoring.
  The default value uses a mixture of the two algorithms, and retaining
  positive-definite working weights.


  }
  \item{imethod}{
  Initialization method used if there are censored observations.
  Currently only the values 1 and 2 are allowed.


  }
  \item{zero, probs.y, lss}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The Weibull density for a response \eqn{Y} is
  \deqn{f(y;a,b) = a y^{a-1} \exp[-(y/b)^a] / (b^a)}{%
        f(y;a,b) = a y^(a-1) * exp(-(y/b)^a) / [b^a]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{y > 0}.
  The cumulative distribution function is
  \deqn{F(y;a,b) = 1 - \exp[-(y/b)^a].}{%
        F(y;a,b) = 1 - exp(-(y/b)^a).}
  The mean of \eqn{Y} is \eqn{b \, \Gamma(1+ 1/a)}{b * gamma(1+ 1/a)}
  (returned as the fitted values),
  and the mode is at \eqn{b\,(1-1/a)^{1/a}}{b * (1- 1/a)^(1/a)} when
  \eqn{a>1}.
  The density is unbounded for \eqn{a<1}.
  The \eqn{k}th moment about the origin is
  \eqn{E(Y^k) = b^k \, \Gamma(1+ k/a)}{E(Y^k) = b^k * gamma(1+ k/a)}.
  The hazard function is \eqn{a t^{a-1} / b^a}{a * t^(a-1) / b^a}.


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


  Summarized in Harper et al. (2011),
  for inference, there are 4 cases to consider.
  If \eqn{a \le 1} then the MLEs are not consistent
  (and the smallest observation becomes a hyperefficient
  solution for the location parameter in the 3-parameter case).
  If \eqn{1 < a < 2} then MLEs exist but are not asymptotically normal.
  If \eqn{a = 2} then the MLEs exist and are normal and asymptotically
  efficient but with a slower convergence rate than when \eqn{a > 2}.
  If \eqn{a > 2} then MLEs have classical asymptotic properties.


  The 3-parameter (location is the third parameter) Weibull can
  be estimated by maximizing a profile log-likelihood (see,
  e.g., Harper et al. (2011) and Lawless (2003)), else try
  \code{\link{gev}} which is a better parameterization.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
\emph{Continuous Univariate Distributions},
2nd edition, Volume 1, New York: Wiley.


Lawless, J. F. (2003)
\emph{Statistical Models and Methods for Lifetime Data},
2nd ed.
{Hoboken, NJ, USA: John Wiley & Sons}.



Rinne, Horst. (2009)
\emph{The Weibull Distribution: A Handbook}.
Boca Raton, FL, USA: CRC Press.



Gupta, R. D. and Kundu, D. (2006)
On the comparison of Fisher information of the
Weibull and GE distributions,
\emph{Journal of Statistical Planning and Inference},
\bold{136},
3130--3144.


Harper, W. V. and Eschenbach, T. G. and James, T. R. (2011)
Concerns about Maximum Likelihood Estimation for
              the Three-Parameter {W}eibull Distribution:
              Case Study of Statistical Software,
\emph{The American Statistician},
\bold{65(1)},
{44--54}.


  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of nonregular cases.
  \emph{Biometrika}, \bold{72}, 67--90.


  Smith, R. L. and Naylor, J. C. (1987)
  A comparison of maximum likelihood and Bayesian estimators
  for the three-parameter Weibull distribution.
  \emph{Applied Statistics}, \bold{36}, 358--369.




}
\author{ T. W. Yee }
\note{
  Successful convergence depends on having reasonably good initial
  values. If the initial values chosen by this function are not good,
  make use the two initial value arguments.


  This \pkg{VGAM} family function handles multiple responses.


  The Weibull distribution is often an alternative to the lognormal
  distribution. The inverse Weibull distribution, which is that of
  \eqn{1/Y} where \eqn{Y} has a Weibull(\eqn{a,b}) distribution, is
  known as the log-Gompertz distribution.


  There are problems implementing the three-parameter Weibull
  distribution. These are because
  the classical regularity conditions for the
  asymptotic properties of the MLEs are not satisfied because the
  support of the distribution depends on one of the parameters.


  Other related distributions are the Maxwell and Rayleigh
  distributions.



}
\section{Warning}{
  This function is under development to handle other censoring situations.
  The version of this function which will handle censored data will be
  called \code{cenweibull()}. It is currently being written and will use
  \code{\link{SurvS4}} as input.
  It should be released in later versions of \pkg{VGAM}.


  If the shape parameter is less than two then misleading inference may
  result, e.g., in the \code{summary} and \code{vcov} of the object.


}

\seealso{
    \code{\link{weibull.mean}},
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
wdata <- transform(wdata,
            y1 = rweibull(nn, shape = exp(1), scale = exp(-2 + x2)),
            y2 = rweibull(nn, shape = exp(2), scale = exp( 1 - x2)))
fit <- vglm(cbind(y1, y2) ~ x2, weibullR, data = wdata, trace = TRUE)
coef(fit, matrix = TRUE)
vcov(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
