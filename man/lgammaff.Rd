\name{lgamma1}
\alias{lgamma1}
\alias{lgamma3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log-gamma Distribution Family Function }
\description{
  Estimation of the parameter of the standard and nonstandard log-gamma
  distribution.

}
\usage{
lgamma1(lshape = "loge", ishape = NULL)
lgamma3(llocation = "identitylink", lscale = "loge", lshape = "loge",
        ilocation = NULL, iscale = NULL, ishape = 1,
        zero = c("scale", "shape"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{
  Parameter link function applied to the
  location parameter \eqn{a}
  and the positive scale parameter \eqn{b}.
  See \code{\link{Links}} for more choices.


  }
  \item{lshape}{
  Parameter link function applied to
  the positive shape parameter \eqn{k}.
  See \code{\link{Links}} for more choices.


  }
  \item{ishape}{
  Initial value for \eqn{k}.
  If given, it must be positive.
  If failure to converge occurs, try some other value.
  The default means an initial value is determined internally.


  }
  \item{ilocation, iscale}{ Initial value for \eqn{a} and \eqn{b}.
  The defaults mean an initial value is determined internally for each.


  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,3\}.
  The default value means none are modelled as intercept-only terms.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The probability density function of the standard log-gamma distribution
  is given by
  \deqn{f(y;k)=\exp[ky - \exp(y)] / \Gamma(k),}{%
        f(y;k) = exp[ky - exp(y)]/gamma(k),}
  for parameter \eqn{k>0}{k>0} and all real \eqn{y}.
  The mean of \eqn{Y} is \code{digamma(k)} (returned as
  the fitted values) and its variance is \code{trigamma(k)}.


  For the non-standard log-gamma distribution, one replaces \eqn{y}
  by \eqn{(y-a)/b}, where \eqn{a} is the location parameter
  and \eqn{b} is the positive scale parameter.
  Then the density function is
    \deqn{f(y)=\exp[k(y-a)/b - \exp((y-a)/b)] / (b \, \Gamma(k)).}{%
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
  The standard log-gamma distribution can be viewed as a
  generalization of the standard type 1 extreme value density:
  when \eqn{k = 1} the distribution of \eqn{-Y} is the standard
  type 1 extreme value distribution.


  The standard log-gamma distribution is fitted with \code{lgamma1}
  and the non-standard (3-parameter) log-gamma distribution is fitted
  with \code{lgamma3}.


}
\seealso{
  \code{\link{rlgamma}},
  \code{\link{gengamma.stacy}},
  \code{\link{prentice74}},
  \code{\link{gamma1}},
  \code{\link[base:Special]{lgamma}}.


}
\examples{
ldata <- data.frame(y = rlgamma(100, shape = exp(1)))
fit <- vglm(y ~ 1, lgamma1, data = ldata, trace = TRUE, crit = "coef")
summary(fit)
coef(fit, matrix = TRUE)
Coef(fit)

ldata <- data.frame(x2 = runif(nn <- 5000))  # Another example
ldata <- transform(ldata, loc = -1 + 2 * x2, Scale = exp(1))
ldata <- transform(ldata, y = rlgamma(nn, loc, scale = Scale, shape = exp(0)))
fit2 <- vglm(y ~ x2, lgamma3, data = ldata, trace = TRUE, crit = "c")
coef(fit2, matrix = TRUE)
}
\keyword{models}
\keyword{regression}
