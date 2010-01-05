\name{inv.gaussianff}
\alias{inv.gaussianff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inverse Gaussian Distribution Family Function }
\description{
  Estimates the two parameters of the inverse Gaussian distribution by
  maximum likelihood estimation.

}
\usage{
inv.gaussianff(lmu="loge", llambda="loge", emu=list(), elambda=list(),
               ilambda=1, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, llambda}{
  Parameter link functions for the \eqn{\mu}{mu} and
  \eqn{\lambda}{lambda} parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{emu, elambda}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ilambda}{ 
  Initial value for the \eqn{\lambda}{lambda} parameter. 

  }
  \item{zero}{ 
  An integer-valued vector specifying which
  linear/additive predictors \eqn{\eta_j}{eta_j}
  are modelled as intercepts only.
  The values must be from the set \{1,2\}.

  }
}
\details{
  The inverse Gaussian distribution has a density
  that can be written as
  \deqn{f(y;\mu,\lambda) = \sqrt{\lambda/(2\pi y^3)}
       \exp\left(-\lambda (y-\mu)^2/(2 \mu^2 y)\right)}{%
       f(y;mu,lambda) = sqrt(lambda/(2*pi*y^3)) *
       exp(-lambda*(y-mu)^2/(2*mu^2*y))
  }
  where \eqn{y>0},
  \eqn{\mu>0}{mu>0}, and
  \eqn{\lambda>0}{lambda>0}.
  The mean of \eqn{Y} is \eqn{\mu}{mu} and its variance is
  \eqn{\mu^3/\lambda}{mu^3/lambda}.
  By default, \eqn{\eta_1=\log(\mu)}{eta1=log(mu)} and
  \eqn{\eta_2=\log(\lambda)}{eta2=log(lambda)}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.

}
\references{ 
Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
\emph{Continuous Univariate Distributions},
2nd edition, Volume 1, New York: Wiley.

Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.
}
\author{ T. W. Yee }
\note{
  The inverse Gaussian distribution can be fitted (to a certain extent)
  using the usual GLM framework involving a scale parameter. This family
  function is different from that approach in that it estimates both
  parameters by full maximum likelihood estimation.

}

\seealso{ 
  \code{\link{Inv.gaussian}},
  \code{\link{wald}},
  \code{\link{bisa}}.

  The \R{} package \pkg{SuppDists} has several functions for evaluating
  the density, distribution function, quantile function and generating
  random numbers from the inverse Gaussian distribution.

}
\examples{
n = 1000
shape = exp(3)
y = rinv.gaussian(n=n, mu=exp(2), lambda=shape)
fit = vglm(y ~ 1, inv.gaussianff(ilam=shape), trace=TRUE)
coef(fit, matrix=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

