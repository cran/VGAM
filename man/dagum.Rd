\name{dagum}
\alias{dagum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Dagum Distribution Family Function }
\description{
  Maximum likelihood estimation of the 3-parameter 
  Dagum distribution.
}
\usage{
dagum(link.a = "loge", link.scale = "loge", link.p = "loge",
      earg.a=list(), earg.scale=list(), earg.p=list(),
      init.a = NULL, init.scale = NULL, init.p = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.a, link.scale, link.p}{
  Parameter link functions applied to the
  (positive) parameters \code{a}, \code{scale}, and \code{p}.
  See \code{\link{Links}} for more choices.

  }
  \item{earg.a, earg.scale, earg.p}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{init.a, init.scale, init.p}{
  Optional initial values for \code{a}, \code{scale}, and \code{p}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2,3\} which correspond to
  \code{a}, \code{scale}, \code{p}, respectively.

  }
}
\details{
  The 3-parameter Dagum distribution is the 4-parameter
  generalized beta II distribution with shape parameter \eqn{q=1}.
It is known under various other names, such as the Burr III, inverse Burr,
beta-K, and 3-parameter kappa distribution.
It can be considered a generalized log-logistic distribution.
Some distributions which are special cases of the 3-parameter Dagum
are the inverse Lomax (\eqn{a=1}), Fisk (\eqn{p=1}),
and the inverse paralogistic (\eqn{a=p}).
  More details can be found in Kleiber and Kotz (2003).


The Dagum distribution has a cumulative distribution function
  \deqn{F(y) = [1 + (y/b)^{-a}]^{-p}}{%
        F(y) = [1 + (y/b)^(-a)]^(-p)}
which leads to a probability density function
  \deqn{f(y) = ap y^{ap-1} / [b^{ap} \{1 + (y/b)^a\}^{p+1}]}{%
        f(y) = ap y^(ap-1) / [b^(ap)  (1 + (y/b)^a)^(p+1)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{p > 0}, \eqn{y > 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and the others are shape parameters.
The mean is
  \deqn{E(Y) = b \, \Gamma(p + 1/a) \, \Gamma(1 - 1/a) / \Gamma(p)}{%
        E(Y) = b  gamma(p + 1/a)  gamma(1 - 1/a) /  gamma(p)}
provided \eqn{-ap < 1 < a}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ: Wiley-Interscience.

}

\author{ T. W. Yee }
\note{
If the self-starting initial values fail, try experimenting
with the initial value arguments, especially those whose
default value is not \code{NULL}.

From Kleiber and Kotz (2003), the MLE is rather sensitive to isolated
observations located sufficiently far from the majority of the data.
Reliable estimation of the scale parameter require \eqn{n>7000},
while estimates for \eqn{a} and \eqn{p} can be considered unbiased for
\eqn{n>2000} or 3000.

}

\seealso{
    \code{\link{Dagum}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{invlomax}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{invparalogistic}}.
}

\examples{
ddata = data.frame(y = rdagum(n=3000, 4, 6, 2))
fit = vglm(y ~ 1, dagum, ddata, trace=TRUE)
fit = vglm(y ~ 1, dagum(init.a=2.1), ddata, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
