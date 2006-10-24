\name{gamma2.ab}
\alias{gamma2.ab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 2-parameter Gamma Distribution }
\description{ Estimates the 2-parameter gamma distribution
  by maximum likelihood estimation.
}
\usage{
gamma2.ab(lrate = "loge", lshape = "loge",
          irate=NULL, ishape=NULL, expected = TRUE, zero = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lrate, lshape}{
  Link functions applied to the (positive) \emph{rate} and \emph{shape}
  parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{expected}{
  Logical. Use Fisher scoring? The default is yes, otherwise
  Newton-Raphson is used.

  }
  \item{irate, ishape}{
  Optional initial values for \emph{rate} and \emph{shape}.
  A \code{NULL} means a value is computed internally.
  If a failure to converge occurs, try using these arguments.

  }
  \item{zero}{
  An integer specifying which
  linear/additive predictor is to be modelled as an intercept only.
  If assigned, the single value should be either 1 or 2 or \code{NULL}.
  The default is to model \eqn{shape} as an intercept only.
  A value \code{NULL} means neither 1 or 2.

  }
}
\details{
  The density function is given by
  \deqn{f(y) = \exp(-rate \times y) \times y^{shape-1} \times rate^{shape} /
               \Gamma(shape)}{%
   f(y) = exp(-rate * y) y^(shape-1) rate^(shape) / gamma(shape)}
  for \eqn{shape > 0}, \eqn{rate > 0} and \eqn{y > 0}.
  Here, \eqn{\Gamma(shape)}{gamma(shape)} is the gamma
  function, as in \code{\link[base:Special]{gamma}}.
  The mean of \emph{Y} is \eqn{\mu=shape/rate}{mu=shape/rate} 
  (returned as the fitted values) with variance
  \eqn{\sigma^2 = \mu^2 /shape = shape/rate^2}{sigma^2 =
       mu^2 /shape = shape/rate^2}.
  By default, the two linear/additive predictors are
  \eqn{\eta_1=\log(rate)}{eta1=log(rate)} and
  \eqn{\eta_2=\log(shape)}{eta2=log(shape)}.

  The argument \code{expected} refers to the type of information
  matrix. The expected information matrix corresponds to Fisher scoring
  and is numerically better here. The observed information matrix
  corresponds to the Newton-Raphson algorithm and may be withdrawn
  from the family function in the future.  If both algorithms work then
  the differences in the results are often not huge.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{
  Most standard texts on statistical distributions describe
  the 2-parameter gamma distribution, e.g.,
  
  Evans, M., Hastings, N. and Peacock, B. (2000)
  \emph{Statistical Distributions},
  New York: Wiley-Interscience, Third edition.

}
\author{ T. W. Yee }
\note{
  The parameters \eqn{rate} and \eqn{shape} match with the arguments
  \code{rate} and \code{shape} of \code{\link[stats]{rgamma}}.  Often,
  \eqn{scale=1/rate} is used.

  If \eqn{rate=1} use the family function \code{\link{gamma1}} to
  estimate \eqn{shape}.

}

\seealso{
  \code{\link{gamma1}} for the 1-parameter gamma distribution,
  \code{\link{gamma2}} for another parameterization of
  the 2-parameter gamma distribution,
  \code{\link{mckaygamma2}} for \emph{a} bivariate gamma distribution,
  \code{\link{expexp}}.

}
\examples{
# Essentially a 1-parameter gamma
y = rgamma(n <- 100, shape= exp(1))
fit1 = vglm(y ~ 1, gamma1, trace=TRUE, crit="c")
fit2 = vglm(y ~ 1, gamma2.ab, trace=TRUE, crit="c")
coef(fit1, matrix=TRUE)
Coef(fit1)
coef(fit2, matrix=TRUE)
Coef(fit2)


# Essentially a 2-parameter gamma
y = rgamma(n <- 500, rate=exp(1), shape=exp(2))
fit2 = vglm(y ~ 1, gamma2.ab, trace=TRUE, crit="c")
coef(fit2, matrix=TRUE)
Coef(fit2)
summary(fit2)
}
\keyword{models}
\keyword{regression}

