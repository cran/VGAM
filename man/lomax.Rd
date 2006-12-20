\name{lomax}
\alias{lomax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lomax Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter 
  Lomax distribution.
}
\usage{
lomax(link.scale = "loge", link.q = "loge",
      earg.scale=list(), earg.q=list(),
      init.scale = NULL, init.q = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.scale, link.q}{
  Parameter link function applied to the
  (positive) parameters \code{scale} and \code{q}.
  See \code{\link{Links}} for more choices.

  }
  \item{earg.scale, earg.q}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{init.scale, init.q}{
  Optional initial values for \code{scale} and \code{q}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2\} which correspond to
  \code{scale}, \code{q}, respectively.

  }
}
\details{
  The 2-parameter Lomax distribution is the 4-parameter
  generalized beta II distribution with shape parameters \eqn{a=p=1}.
It is probably more widely known as the Pareto (II) distribution.
  It is also the 3-parameter Singh-Maddala distribution
  with shape parameter \eqn{a=1}, as well as the
  beta distribution of the second kind with \eqn{p=1}.
  More details can be found in Kleiber and Kotz (2003).

The Lomax distribution has density
  \deqn{f(y) = q / [b \{1 + y/b\}^{1+q}]}{%
        f(y) = q / [b (1 + y/b)^(1+q)]}
  for \eqn{b > 0}, \eqn{q > 0}, \eqn{y > 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and \code{q} is a shape parameter.
The cumulative distribution function is
  \deqn{F(y) = 1 - [1 + (y/b)]^{-q}.}{%
        F(y) = 1 - [1 + (y/b)]^(-q).}
The mean is
  \deqn{E(Y) = b/(q-1)}{%
        E(Y) = b/(q-1)}
provided \eqn{q > 1}.

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
    \code{\link{Lomax}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{invlomax}},
    \code{\link{paralogistic}},
    \code{\link{invparalogistic}}.
}

\examples{
y = rlomax(n=2000, 6, 2)
fit = vglm(y ~ 1, lomax, trace=TRUE)
fit = vglm(y ~ 1, lomax, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
