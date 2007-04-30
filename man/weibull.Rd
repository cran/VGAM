\name{weibull}
\alias{weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Weibull Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter Weibull distribution.
  Allows for Type-I right censored data.

}
\usage{
weibull(lshape = "logoff", lscale = "loge", 
        eshape=if(lshape == "logoff") list(offset=-2) else list(),
        escale=list(),
        ishape = NULL, iscale = NULL, imethod=1, zero = 2)
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

  This \pkg{VGAM} family function handles Type-I right censored data as
  well as complete data.
  Fisher scoring is used to estimate the two parameters.
  The Fisher information matrices used here are only valid
  if \eqn{a>2}; these are where the regularity conditions for maximum
  likelihood estimation are satisfied.
  For this reason, the default link function for the shape parameter is
  a log-link with an offset value of \eqn{-2}.

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
  For censored data, numerical integration is used to compute the
  expected working weight matrices; this may fail if the data is
  censored `too much' and/or may be quite slow.
  See the example below on how to input censored data.

  The Weibull distribution is often an alternative to the lognormal
  distribution.  The inverse Weibull distribution, which is that of
  \eqn{1/Y} where \eqn{Y} has a Weibull(\eqn{a,b}) distribution, is
  known as the log-Gompertz distribution.

}
\section{Warning}{
  If the shape parameter is less than two then numerical problems may
  occur during the fitting and/or misleading inference may result in
  the \code{summary} of the object.

}

\seealso{
    \code{\link[stats:Weibull]{dweibull}},
    \code{\link{gev}},
    \code{\link{lognormal}},
    \code{\link{expexp}}.
}
\examples{
# Complete data
x = runif(n <- 1000)
y = rweibull(n, shape=2+exp(1+x), scale = exp(-0.5))
fit = vglm(y ~ x, weibull, tra=TRUE)
coef(fit, mat=TRUE)
Coef(fit)

# Type-I right censored data
cutpt = 0.6 # Making this too small results in numerical problems
rcensored = y > cutpt
cy = ifelse(rcensored, cutpt, y)
table(rcensored)
\dontrun{
par(mfrow=1:2)
hist(y, xlim=range(y))
hist(cy, xlim=range(y), main="Censored y")
}
cfit = vglm(cy ~ x, weibull, trace=TRUE, crit="l",
            extra=list(rightcensored=rcensored))
coef(cfit, mat=TRUE)
summary(cfit)
}
\keyword{models}
\keyword{regression}
