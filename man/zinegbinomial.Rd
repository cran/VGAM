\name{zinegbinomial}
\alias{zinegbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Negative Binomial Distribution Family Function }
\description{
  Fits a zero-inflated negative binomial distribution by
  full maximum likelihood estimation.

}
\usage{
zinegbinomial(lphi = "logit", lmunb = "loge", lsize = "loge", ephi = list(),
      emunb = list(), esize = list(), iphi = NULL, isize = NULL, zero = -3,
      imethod = 1, shrinkage.init = 0.95, nsimEIM = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi, lmunb, lsize}{ 
  Link functions for the parameters \eqn{\phi}{phi},
  the mean and \eqn{k}; see \code{\link{negbinomial}} for details,
  and \code{\link{Links}} for more choices.

  }
  \item{ephi, emunb, esize}{
  List. Extra arguments for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.

  }

  \item{iphi, isize}{
  Optional initial values for \eqn{\phi}{phi} and \eqn{k}{k}.
  The default is to compute an initial value internally for both.
  If a vector then recycling is used.

  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} or \code{3} which
  specifies the initialization method for the mean parameter.
  If failure to converge occurs try another value
  and/or else specify a value for \code{shrinkage.init}.

  }
  \item{zero}{ 
  Integers specifying which linear/additive predictor is modelled
  as intercepts only.  If given, their absolute values must be
  either 1 or 2 or 3, and the default is only the \eqn{k} parameters
  (one for each response).
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{shrinkage.init, nsimEIM}{ 
  See \code{\link{CommonVGAMffArguments}} for information.

  }
}
\details{
  This function uses simulation and Fisher scoring, and is based on
  \deqn{P(Y=0) =  \phi + (1-\phi) (k/(k+\mu))^k,}{%
        P(Y=0) =  \phi + (1-\phi) * (k/(k+\mu))^k,}
  and for \eqn{y=1,2,\ldots},
  \deqn{P(Y=y) =  (1-\phi) \, dnbinom(y, \mu, k).}{%
        P(Y=y) =  (1-\phi) * dnbinom(y, \mu, k).}
  The parameter \eqn{\phi}{phi} satisfies \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is \eqn{(1-\phi) \mu}{(1-phi)*munb}
  (returned as the fitted values).  By default, the three linear/additive
  predictors are \eqn{(logit(\phi), \log(\mu), \log(k))^T}{(logit(phi),
  log(munb), log(k))^T}.
  See \code{\link{negbinomial}}, another \pkg{VGAM} family function,
  for the formula of the probability density function and other details
  of the negative binomial distribution.


  Independent multivariate responses are handled.
  If so then arguments \code{iphi} and \code{isize} may be vectors
  with length equal to the number of responses.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
%\references{
%
%}
\author{ T. W. Yee }
\note{
  For intercept-models, the \code{misc} slot has a component called
  \code{p0} which is the estimate of \eqn{P(Y=0)}. Note that \eqn{P(Y=0)}
  is not the parameter \eqn{\phi}{phi}.


  If \eqn{k} is large then the use of \pkg{VGAM} family function
  \code{\link{zipoisson}} is probably preferable.
  This follows because the Poisson is the limiting distribution of a
  negative binomial as \eqn{k} tends to infinity.


  The zero-\emph{deflated} negative binomial distribution cannot
  be handled with this family function. It can be handled
  with the zero-altered negative binomial distribution; see
  \code{\link{zanegbinomial}}.


}

\section{Warning }{
  Numerical problems can occur, e.g., when the probability of
  zero is actually less than, not more than, the nominal
  probability of zero.
  Half-stepping is not uncommon.
  If failure to converge occurs, try using combinations of
  \code{imethod},
  \code{shrinkage.init},
  \code{iphi},
  \code{isize}, and/or
  \code{zero} if there are explanatory variables.

} 

\seealso{
  \code{\link{Zinegbin}},
  \code{\link{negbinomial}},
  \code{\link[stats:Poisson]{rpois}},
  \code{\link{CommonVGAMffArguments}}.

}
\examples{
# Example 1
nbdat <- data.frame(x2 = runif(nn <- 1000))
nbdat <- transform(nbdat, phi  = logit(-0.5 + 1 * x2, inverse = TRUE),
                          munb = exp(3 + x2),
                          size = exp(0 + 2*x2))
nbdat <- transform(nbdat, y1 = rzinegbin(nn, phi, mu = munb, size = size),
                          y2 = rzinegbin(nn, phi, mu = munb, size = size))
with(nbdat, table(y1)["0"] / sum(table(y1)))
fit <- vglm(cbind(y1, y2) ~ x2, zinegbinomial(zero = NULL), nbdat)
coef(fit, matrix = TRUE)
summary(fit)
head(cbind(fitted(fit), with(nbdat, (1-phi) * munb)))
round(vcov(fit), 3)


# Example 2: RR-ZINB could also be called a COZIVGLM-ZINB-2
nbdat <- data.frame(x2 = runif(nn <- 2000))
nbdat <- transform(nbdat, x3 = runif(nn))
nbdat <- transform(nbdat, eta1 = 3 + 1 * x2 + 2 * x3)
nbdat <- transform(nbdat, phi  = logit(-1.5 + 0.5 * eta1, inverse = TRUE),
                          munb = exp(eta1),
                          size = exp(4))
nbdat <- transform(nbdat, y1 = rzinegbin(nn, phi, mu = munb, size = size))
with(nbdat, table(y1)["0"] / sum(table(y1)))
rrzinb <- rrvglm(y1 ~ x2 + x3, zinegbinomial(zero = NULL), nbdat,
                 Index.corner = 2, szero = 3, trace = TRUE)
coef(rrzinb, matrix = TRUE)
Coef(rrzinb)
}
\keyword{models}
\keyword{regression}

