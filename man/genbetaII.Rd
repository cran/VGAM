\name{genbetaII}
\alias{genbetaII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Beta Distribution of the Second Kind }
\description{
  Maximum likelihood estimation of the 4-parameter 
  generalized beta II distribution.
}
\usage{
genbetaII(link.a = "loge", link.scale = "loge",
          link.p = "loge", link.q = "loge",
          earg.a=list(), earg.scale=list(), earg.p=list(), earg.q=list(),
          init.a = NULL, init.scale = NULL, init.p = 1, init.q = 1,
          zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.a, link.scale, link.p, link.q}{
  Parameter link functions applied to the
  shape parameter \code{a},
  scale parameter \code{scale},
  shape parameter \code{p}, and
  shape parameter \code{q}.
  All four parameters are positive.
  See \code{\link{Links}} for more choices.

  }
  \item{earg.a, earg.scale, earg.p, earg.q}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{init.a, init.scale}{
  Optional initial values for \code{a} and \code{scale}.
  A \code{NULL} means a value is computed internally.

  }
  \item{init.p, init.q}{
  Optional initial values for \code{p} and \code{q}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2,3,4\} which correspond to
  \code{a}, \code{scale}, \code{p}, \code{q}, respectively.

  }
}
\details{
  This distribution is most useful for unifying a substantial
  number of size distributions. For example, the Singh-Maddala,
  Dagum, Fisk (log-logistic), Lomax (Pareto type II),
  inverse Lomax, beta distribution of the second kind
  distributions are all special cases.
  Full details can be found in Kleiber and Kotz (2003), and
  Brazauskas (2002).

The 4-parameter generalized beta II distribution has density
  \deqn{f(y) = a y^{ap-1} / [b^{ap} B(p,q) \{1 + (y/b)^a\}^{p+q}]}{%
        f(y) = a y^(ap-1) / [b^(ap) B(p,q) (1 + (y/b)^a)^(p+q)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{p > 0}, \eqn{q > 0}, \eqn{y > 0}.
Here \eqn{B} is the beta function, and 
\eqn{b} is the scale parameter \code{scale},
while the others are shape parameters.
The mean is 
  \deqn{E(Y) = b \, \Gamma(p + 1/a) \, \Gamma(q - 1/a) / (\Gamma(p) \, \Gamma(q))}{%
        E(Y) = b  gamma(p + 1/a)  gamma(q - 1/a) / ( gamma(p)  gamma(q))}
provided \eqn{-ap < 1 < aq}.

%The distribution is motivated by the incomplete beta function
%\eqn{B_y(p,q)} which is the integral from 0 to \eqn{y} of the integrand
%\eqn{u^{p-1} (1-u)^{q-1}}{u^(p-1) (1-u)^(q-1)} where \eqn{y>0}.

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

Brazauskas, V. (2002)
Fisher information matrix for the Feller-Pareto distribution.
\emph{Statistics & Probability Letters},
\bold{59}, 159--167.

}

\author{ T. W. Yee }
\note{
If the self-starting initial values fail, try experimenting
with the initial value arguments, especially those whose
default value is not \code{NULL}.

Successful convergence depends on having very
good initial values. This is rather difficult for this distribution!
More improvements could be made here.

}

\seealso{
    \code{\link{lino}},
    \code{\link{betaff}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{lomax}},
    \code{\link{invlomax}},
    \code{\link{paralogistic}},
    \code{\link{invparalogistic}}.
}

\examples{
gdata = data.frame(y = rsinmad(n=3000, 4, 6, 2)) # Not very good data!
fit = vglm(y ~ 1, genbetaII, gdata, trace=TRUE)
fit = vglm(y ~ 1, genbetaII(init.p=1.0, init.a=4, init.sc=7, init.q=2.3),
           gdata, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
