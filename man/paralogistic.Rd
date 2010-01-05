\name{paralogistic}
\alias{paralogistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Paralogistic Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter 
  paralogistic distribution.
}
\usage{
paralogistic(link.a = "loge", link.scale = "loge", earg.a=list(),
             earg.scale=list(), init.a = 1, init.scale = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.a, link.scale}{
  Parameter link functions applied to the
  (positive) shape parameter \code{a} and
  (positive) scale parameter \code{scale}.
  See \code{\link{Links}} for more choices.

  }
  \item{earg.a, earg.scale}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{init.a, init.scale}{
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
provided \eqn{a > 1}.

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
pdat = data.frame(y = rparalogistic(n=3000, 4, 6))
fit = vglm(y ~ 1, paralogistic, pdat, trace=TRUE)
fit = vglm(y ~ 1, paralogistic(init.a=2.3, init.sc=5),
           pdat, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
