\name{sinmad}
\alias{sinmad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Singh-Maddala Distribution Family Function }
\description{
  Maximum likelihood estimation of the 3-parameter 
  Singh-Maddala distribution.
}
\usage{
sinmad(link.a = "loge", link.scale = "loge", link.q = "loge",
       init.a = NULL, init.scale = NULL, init.q = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.a, link.scale, link.q}{
  Parameter link functions applied to the
  (positive) parameters \code{a}, \code{scale}, and \code{q}.
  See \code{\link{Links}} for more choices.

  }
  \item{init.a, init.scale, init.q}{
  Optional initial values for \code{a}, \code{scale}, and \code{q}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2,3\} which correspond to
  \code{a}, \code{scale}, \code{q}, respectively.

  }
}
\details{
  The 3-parameter Singh-Maddala distribution is the 4-parameter
  generalized beta II distribution with shape parameter \eqn{p=1}.
It is known under various other names, such as the Burr XII (or
just the Burr distribution), Pareto IV,
beta-P, and generalized log-logistic distribution.
  More details can be found in Kleiber and Kotz (2003).

Some distributions which are special cases of the 3-parameter Singh-Maddala
are the Lomax (\eqn{a=1}), Fisk (\eqn{q=1}), and
paralogistic (\eqn{a=q}).

The Singh-Maddala distribution has density
  \deqn{f(y) = aq y^{a-1} / [b^a \{1 + (y/b)^a\}^{1+q}]}{%
        f(y) = aq y^(a-1) / [b^a (1 + (y/b)^a)^(1+q)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{q > 0}, \eqn{y > 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and the others are shape parameters.
The cumulative distribution function is
  \deqn{F(y) = 1 - [1 + (y/b)^a]^{-q}.}{%
        F(y) = 1 - [1 + (y/b)^a]^(-q).}
The mean is
  \deqn{E(Y) = b \, \Gamma(1 + 1/a) \, \Gamma(q - 1/a) / \Gamma(q)}{%
        E(Y) = b  gamma(1 + 1/a)  gamma(q - 1/a) /  gamma(q)}
provided \eqn{-a < 1 < aq}.

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
    \code{\link{Sinmad}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{fisk}},
    \code{\link{invlomax}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{invparalogistic}}.
}

\examples{
y = rsinmad(n=3000, 3, 5, 2)
fit = vglm(y ~ 1, sinmad, trace=TRUE)
fit = vglm(y ~ 1, sinmad, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
