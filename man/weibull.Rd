\name{weibull}
\alias{weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Weibull Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter Weibull distribution.
  No observations should be censored.

}
\usage{
weibull(lshape = "loge", lscale = "loge", 
        eshape = list(), escale = list(),
        ishape = NULL,   iscale = NULL, nrfs = 1, imethod=1, zero = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, lscale}{
  Parameter link functions applied to the 
  (positive) shape parameter (called \eqn{a} below) and
  (positive) scale parameter (called \eqn{b} below).
  See \code{\link{Links}} for more choices.

  }
  \item{eshape, escale}{
  Extra argument for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  Parameter link functions applied to the 
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
  \item{zero}{
  An integer specifying which linear/additive predictor is to be modelled
  as an intercept only.  The value must be from the set \{1,2\},
  which correspond to the shape and scale parameters respectively.
  Setting \code{zero=NULL} means none of them.

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
  Although the Fisher information matrices used here are valid
  in all regions of the parameter space,
  the regularity conditions for maximum
  likelihood estimation are satisfied only if \eqn{a>2}
  (according to Kleiber and Kotz (2003)).
  If this is violated then a warning message is issued.
  One can enforce \eqn{a>2} by choosing \code{lshape = "logoff"}
  and \code{eshape=list(offset=-2)}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and Actuarial Sciences},
Hoboken, NJ: Wiley-Interscience.


Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
\emph{Continuous Univariate Distributions},
2nd edition, Volume 1, New York: Wiley.


Gupta, R. D. and Kundu, D. (2006)
On the comparison of Fisher information of the
Weibull and GE distributions,
\emph{Journal of Statistical Planning and Inference},
\bold{136},
3130--3144.


}
\author{ T. W. Yee }
\note{
  Successful convergence depends on having reasonably good initial
  values. If the initial values chosen by this function are not good,
  make use the two initial value arguments.


  The Weibull distribution is often an alternative to the lognormal
  distribution.  The inverse Weibull distribution, which is that of
  \eqn{1/Y} where \eqn{Y} has a Weibull(\eqn{a,b}) distribution, is
  known as the log-Gompertz distribution.

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
    \code{\link[stats:Weibull]{dweibull}},
    \code{\link{gev}},
    \code{\link{lognormal}},
    \code{\link{expexp}}.
}
\examples{
# Complete data
wdata = data.frame(x = runif(nn <- 1000))
wdata = transform(wdata, y = rweibull(nn, shape = exp(1+x), scale = exp(-2)))
fit = vglm(y ~ x, weibull, wdata, trace = TRUE)
coef(fit, mat = TRUE)
vcov(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
