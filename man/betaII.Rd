\name{betaII}
\alias{betaII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Beta Distribution of the Second Kind }
\description{
  Maximum likelihood estimation of the 3-parameter 
  beta II distribution.
}
\usage{
betaII(link.scale = "loge", link.p = "loge", link.q = "loge",
       earg.scale=list(), earg.p=list(), earg.q=list(),
       init.scale = NULL, init.p = 1, init.q = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.scale, link.p, link.q}{
  Parameter link functions applied to the
  (positive) parameters \code{scale}, \code{p} and \code{q}.
  See \code{\link{Links}} for more choices.

  }
  \item{earg.scale, earg.p, earg.q}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{init.scale, init.p, init.q}{
  Optional initial values for \code{scale}, \code{p} and \code{q}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2,3\} which correspond to
  \code{scale}, \code{p}, \code{q}, respectively.

  }
}
\details{
  The 3-parameter beta II is the 4-parameter
  \emph{generalized} beta II distribution with shape parameter \eqn{a=1}.
  It is also known as the Pearson VI distribution.
  Other distributions which are special cases of the 3-parameter beta II include
  the Lomax (\eqn{p=1}) and inverse Lomax (\eqn{q=1}).
  More details can be found in Kleiber and Kotz (2003).

The beta II distribution has density
  \deqn{f(y) = y^{p-1} / [b^p B(p,q) \{1 + y/b\}^{p+q}]}{%
        f(y) = y^(p-1) / [b^p B(p,q) (1 + y/b)^(p+q)]}
  for \eqn{b > 0}, \eqn{p > 0}, \eqn{q > 0}, \eqn{y > 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and the others are shape parameters.
The mean is 
  \deqn{E(Y) = b \, \Gamma(p + 1) \, \Gamma(q - 1) / (\Gamma(p) \, \Gamma(q))}{%
        E(Y) = b  gamma(p + 1)  gamma(q - 1) / ( gamma(p)  gamma(q))}
provided \eqn{q > 1}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and Actuarial Sciences},
Hoboken, NJ: Wiley-Interscience.

}

\author{ T. W. Yee }
\note{
    If the self-starting initial values fail, try experimenting with
    the initial value arguments, especially those whose default value
    is not \code{NULL}.


}

\seealso{
    \code{\link{betaff}},
    \code{\link{genbetaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{invlomax}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{invparalogistic}}.

}

\examples{
bdata = data.frame(y = rsinmad(n=2000, a=1, 6, 2))  # Not genuine data!
fit = vglm(y ~ 1, betaII, bdata, trace=TRUE)
fit = vglm(y ~ 1, betaII(init.p=0.7, init.q=0.7), bdata, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
