\name{invlomax}
\alias{invlomax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inverse Lomax Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter 
  inverse Lomax distribution.
}
\usage{
invlomax(link.scale = "loge", link.p = "loge",
         init.scale = NULL, init.p = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.scale, link.p}{
  Parameter link functions applied to the
  (positive) scale parameter \code{scale} and
  (positive) shape parameter \code{p}.
  See \code{\link{Links}} for more choices.

  }
  \item{init.scale, init.p}{
  Optional initial values for \code{scale} and \code{p}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2\} which correspond to
  \code{scale}, \code{p}, respectively.

  }
}
\details{
  The 2-parameter inverse Lomax distribution is the 4-parameter
  generalized beta II distribution with shape parameters \eqn{a=q=1}.
  It is also the 3-parameter Dagum distribution
  with shape parameter \eqn{a=1}, as well as the
  beta distribution of the second kind with \eqn{q=1}.
  More details can be found in Kleiber and Kotz (2003).

The inverse Lomax distribution has density
  \deqn{f(y) = p y^{p-1} / [b^p \{1 + y/b\}^{p+1}]}{%
        f(y) = p y^(p-1) / [b^p (1 + y/b)^(p+1)]}
  for \eqn{b > 0}, \eqn{p > 0}, \eqn{y > 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and \code{p} is a shape parameter.
The mean does not seem to exist. 

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
    \code{\link{Invlomax}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{invparalogistic}}.
}

\examples{
y = rinvlomax(n=2000, 6, 2)
fit = vglm(y ~ 1, invlomax, trace=TRUE)
fit = vglm(y ~ 1, invlomax, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
