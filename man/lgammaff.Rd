\name{lgammaff}
\alias{lgammaff}
\alias{lgamma3ff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log-gamma Distribution Family Function }
\description{
  Estimation of the parameter of the standard and nonstandard log-gamma
  distribution.

}
\usage{
lgammaff(link = "loge", earg=list(), init.k = NULL)
lgamma3ff(llocation="identity", lscale="loge", lshape="loge",
          elocation=list(), escale=list(), eshape=list(),
          ilocation=NULL, iscale=NULL, ishape=1, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation}{
  Parameter link function applied to the
  location parameter \eqn{a}.
  See \code{\link{Links}} for more choices.

  }
  \item{lscale}{
  Parameter link function applied to the
  positive scale parameter \eqn{b}.
  See \code{\link{Links}} for more choices.

  }
  \item{link, lshape}{
  Parameter link function applied to 
  the positive shape parameter \eqn{k}. 
  See \code{\link{Links}} for more choices.

  }
  \item{earg, elocation, escale, eshape}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{init.k, ishape}{
  Initial value for \eqn{k}.
  If given, it must be positive. 
  If failure to converge occurs, try some other value.
  The default means an initial value is determined internally. }
  \item{ilocation, iscale}{ Initial value for \eqn{a} and \eqn{b}.
  The defaults mean an initial value is determined internally for each.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,3\}.
  The default value means none are modelled as intercept-only terms.

  }
}
\details{
  The probability density function of the standard log-gamma distribution
  is given by
  \deqn{f(y)=\exp[ky - \exp(y)] / \Gamma(k),}{%
        f(y) = exp[ky - exp(y)]/gamma(k),}
  for parameter \eqn{k>0}{k>0} and all real \eqn{y}.
  The mean of \eqn{Y} is \code{digamma(k)} (returned as
  the fitted values) and its variance is \code{trigamma(k)}.

  For the non-standard log-gamma distribution, one replaces \eqn{y}
  by \eqn{(y-a)/b}, where \eqn{a} is the location parameter
  and \eqn{b} is the positive scale parameter.
  Then the density function is 
    \deqn{f(y)=\exp[k(y-a)/b - \exp((y-a)/b)] / (b \Gamma(k)).}{%
          f(y) = exp[k(y-a)/b - exp((y-a)/b)]/(b*gamma(k)).}
  The mean and variance of \eqn{Y} are \code{a + b*digamma(k)} (returned as
  the fitted values) and \code{b^2 * trigamma(k)}, respectively.

}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{
Kotz, S. and Nadarajah, S. (2000)
\emph{Extreme Value Distributions: Theory and Applications},
pages 48--49,
London: Imperial College Press.

Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1995)
\emph{Continuous Univariate Distributions},
2nd edition, Volume 2, p.89,
New York: Wiley.

}

\author{ T. W. Yee }
\note{ 
  The standard log-gamma distribution can be viewed as a generalization
  of the standard type 1 extreme value density: when \eqn{k=1}
  the distribution of \eqn{-Y} is the standard type 1 extreme value
  distribution.

  The standard log-gamma distribution is fitted with \code{lgammaff}
  and the non-standard (3-parameter) log-gamma distribution is fitted
  with \code{lgamma3ff}.

}
\seealso{
\code{\link{rlgamma}},
\code{\link{ggamma}},
\code{\link{prentice74}},
\code{\link{gamma1}},
\code{\link[base:Special]{lgamma}}.
}
\examples{
y = rlgamma(n <- 100, k=exp(1))
fit = vglm(y ~ 1, lgammaff, trace=TRUE, crit="c")
summary(fit)
coef(fit, matrix=TRUE)
Coef(fit)

# Another example
x = runif(n <- 5000)
loc = -1 + 2*x
Scale = exp(1+x)
y = rlgamma(n, loc=loc, scale=Scale, k=exp(0))
fit = vglm(y ~ x, lgamma3ff(zero=3), trace=TRUE, crit="c")
coef(fit, matrix=TRUE)
}
\keyword{models}
\keyword{regression}
