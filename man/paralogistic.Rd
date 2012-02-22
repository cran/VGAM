\name{paralogistic}
\alias{paralogistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Paralogistic Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter 
  paralogistic distribution.
}
\usage{
paralogistic(lshape1.a = "loge", lscale = "loge", eshape1.a = list(),
             escale = list(), ishape1.a = 2, iscale = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1.a, lscale}{
  Parameter link functions applied to the
  (positive) shape parameter \code{a} and
  (positive) scale parameter \code{scale}.
  See \code{\link{Links}} for more choices.

  }
  \item{eshape1.a, escale}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ishape1.a, iscale}{
  Optional initial values for \code{a} and \code{scale}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2\} which correspond to
  \code{a}, \code{scale}, respectively.

  }
}
\details{
  The 2-parameter paralogistic distribution is the 4-parameter
  generalized beta II distribution with shape parameter \eqn{p=1} and
  \eqn{a=q}.
It is the 3-parameter Singh-Maddala distribution with \eqn{a=q}.
  More details can be found in Kleiber and Kotz (2003).


The 2-parameter paralogistic has density
  \deqn{f(y) = a^2 y^{a-1} / [b^a \{1 + (y/b)^a\}^{1+a}]}{%
        f(y) = a^2 y^(a-1) / [b^a (1 + (y/b)^a)^(1+a)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{y > 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and \eqn{a} is the shape parameter.
The mean is
  \deqn{E(Y) = b \, \Gamma(1 + 1/a) \, \Gamma(a - 1/a) / \Gamma(a)}{%
        E(Y) = b  gamma(1 + 1/a)  gamma(a - 1/a) /  gamma(a)}
provided \eqn{a > 1}; these are returned as the fitted values.


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


}

\author{ T. W. Yee }
\note{
If the self-starting initial values fail, try experimenting
with the initial value arguments, especially those whose
default value is not \code{NULL}.


}

\seealso{
    \code{\link{Paralogistic}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{fisk}},
    \code{\link{invlomax}},
    \code{\link{lomax}},
    \code{\link{invparalogistic}}.

}

\examples{
pdata = data.frame(y = rparalogistic(n = 3000, 4, 6))
fit = vglm(y ~ 1, paralogistic, pdata, trace = TRUE)
fit = vglm(y ~ 1, paralogistic(ishape1.a = 2.3, iscale = 5),
           pdata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
