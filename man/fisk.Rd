\name{fisk}
\alias{fisk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fisk Distribution family function }
\description{
  Maximum likelihood estimation of the 2-parameter 
  Fisk distribution.
}
\usage{
fisk(link.a = "loge", link.scale = "loge",
     init.a = NULL, init.scale = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.a, link.scale}{
  Parameter link functions applied to the
  (positive) parameters \code{a} and \code{scale}.
  See \code{\link{Links}} for more choices.

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
  The 2-parameter Fisk (aka log-logistic) distribution is the 4-parameter
  generalized beta II distribution with shape parameter \eqn{q=p=1}.
  It is also the 3-parameter Singh-Maddala distribution
  with shape parameter \eqn{q=1}, as well as the 
  Dagum distribution with \eqn{p=1}.
  More details can be found in Kleiber and Kotz (2003).

The Fisk distribution has density
  \deqn{f(y) = a y^{a-1} / [b^a \{1 + (y/b)^a\}^2]}{%
        f(y) = a y^(a-1) / [b^a (1 + (y/b)^a)^2]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{y > 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and \code{a} is a shape parameter.
The cumulative distribution function is
  \deqn{F(y) = 1 - [1 + (y/b)^a]^{-1} = [1 + (y/b)^{-a}]^{-1}.}{%
        F(y) = 1 - [1 + (y/b)^a]^(-1) = [1 + (y/b)^(-a)]^(-1).}
The mean is
  \deqn{E(Y) = b \, \Gamma(1 + 1/a) \, \Gamma(1 - 1/a)}{%
        E(Y) = b  gamma(1 + 1/a)  gamma(1 - 1/a)}
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
    \code{\link{Fisk}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{invlomax}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{invparalogistic}}.
}

\examples{
y = rfisk(n=200, 4, 6)
fit = vglm(y ~ 1, fisk, trace=TRUE)
fit = vglm(y ~ 1, fisk(init.a=3.3), trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
