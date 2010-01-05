\name{zipoisson}
\alias{zipoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Poisson Distribution Family Function }
\description{
  Fits a zero-inflated Poisson distribution by full maximum likelihood
  estimation.

}
\usage{
zipoisson(lphi="logit", llambda = "loge", ephi=list(), elambda =list(),
          iphi = NULL, method.init=1, shrinkage.init=0.8, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi}{ 
  Link function for the parameter \eqn{\phi}{phi}.
  See \code{\link{Links}} for more choices.

  }
  \item{llambda}{ 
  Link function for the usual \eqn{\lambda}{lambda} parameter.
  See \code{\link{Links}} for more choices.

  }
  \item{ephi, elambda}{
  List. Extra argument for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.

  }

  \item{iphi}{
  Optional initial values for \eqn{\phi}{phi}, whose values must lie
  between 0 and 1. The default is to compute an initial value internally.
  If a vector then recycling is used.

  }
  \item{method.init}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method for \eqn{\lambda}{lambda}.
  If failure to converge occurs try another value
  and/or else specify a value for \code{shrinkage.init}
  and/or else specify a value for \code{iphi}.

  }
  \item{shrinkage.init}{
  How much shrinkage is used when initializing \eqn{\lambda}{lambda}.
  The value must be between 0 and 1 inclusive, and 
  a value of 0 means the individual response values are used,
  and a value of 1 means the median or mean is used.
  This argument is used in conjunction with \code{method.init}.

  }
  \item{zero}{ 
  An integer specifying which linear/additive predictor is modelled as
  intercepts only.  If given, the value must be either 1 or 2, and the
  default is none of them. Setting \code{zero=1} makes \eqn{\phi}{phi}
  a single parameter.

  }
}
\details{
  The model is a mixture of a Poisson distribution and the value 0;
  it has value 0 with probability \eqn{\phi}{phi} else is
  Poisson(\eqn{\lambda}{lambda}) distributed.
  The model can be written
  \deqn{P(Y=0) =  \phi + (1-\phi) \exp(-\lambda),}{%
        P(Y=0) =  phi + (1-phi) * exp(-lambda),}
  and for \eqn{y=1,2,\ldots},
  \deqn{P(Y=y) =  (1-\phi) \exp(-\lambda) \lambda^y / y!.}{%
        P(Y=y) =  (1-phi) * exp(-lambda) * lambda^y / y!.}
  Here, the parameter \eqn{\phi}{phi} satisfies
  \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is \eqn{(1-\phi) \lambda}{(1-phi)*lambda} and these
  are returned as the fitted values.  By default, the two linear/additive
  predictors are \eqn{(logit(\phi), \log(\lambda))^T}{(logit(phi),
  log(lambda))^T}.
  This function implements Fisher scoring.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{
  Thas, O. and Rayner, J. C. W. (2005)
  Smooth tests for the zero-inflated Poisson distribution.
  \emph{Biometrics},
  \bold{61}, 808--815.

  Data: Angers, J-F. and Biswas, A. (2003)
  A Bayesian analysis of zero-inflated generalized Poisson model.
  \emph{Computational Statistics & Data Analysis},
  \bold{42}, 37--46.

  Cameron, A. C. and Trivedi, P. K. (1998)
  \emph{Regression Analysis of Count Data}.
  Cambridge University Press: Cambridge.

}
\author{ T. W. Yee }
\note{
  For intercept-models, the \code{misc} slot has a component called
  \code{p0} which is the estimate of \eqn{P(Y=0)}. Note that \eqn{P(Y=0)}
  is not the parameter \eqn{\phi}{phi}.  This family function currently
  cannot handle a multivariate response.

  This family function is now recommended above \code{\link{yip88}}.

  The zero-\emph{deflated} Poisson distribution cannot be handled with
  this family function. It can be handled with the zero-altered Poisson
  distribution; see \code{\link{zapoisson}}.

}

\section{Warning }{
  Numerical problems can occur, e.g., when the probability of
  zero is actually less than, not more than, the nominal
  probability of zero.
  For example, in the Angers and Biswas (2003) data below,
  replacing 182 by 1 results in nonconvergence.
  Half-stepping is not uncommon.
  If failure to converge occurs, try using combinations of
  \code{method.init},
  \code{shrinkage.init},
  \code{iphi}, and/or
  \code{zero=1} if there are explanatory variables.

} 

\seealso{
  \code{\link{zapoisson}},
  \code{\link{Zipois}},
  \code{\link{yip88}},
  \code{\link{zipebcom}},
  \code{\link[stats:Poisson]{rpois}}.

}
\examples{
zipdat = data.frame(x = runif(nn <- 2000))
zipdat = transform(zipdat, phi = logit(-0.5 + 1*x, inverse=TRUE),
                           lambda = loge(0.5 + 2*x, inverse=TRUE))
zipdat = transform(zipdat, y = rzipois(nn, lambda, phi))
with(zipdat, table(y))
fit = vglm(y ~ x, zipoisson, zipdat, trace=TRUE)
coef(fit, matrix=TRUE)  # These should agree with the above values


# Another example: McKendrick (1926). Data from 223 Indian village households
cholera = data.frame(ncases = 0:4, # Number of cholera cases,
                     wfreq = c(168, 32, 16, 6, 1)) # Frequencies
fit = vglm(ncases ~ 1, zipoisson, wei=wfreq, data=cholera, trace=TRUE)
coef(fit, matrix=TRUE)
with(cholera, cbind(actual=wfreq, fitted= round(
      dzipois(ncases, lambda=Coef(fit)[2], phi=Coef(fit)[1]) * sum(wfreq), dig=2)))


# Another example: data from Angers and Biswas (2003)
abdat = data.frame(y = 0:7, w = c(182, 41, 12, 2, 2, 0, 0, 1))
abdat = subset(abdat, w>0)
fit = vglm(y ~ 1, zipoisson(lphi=probit, iphi=0.3), abdat, wei=w, trace=TRUE)
fit@misc$prob0  # Estimate of P(Y=0)
coef(fit, matrix=TRUE)
Coef(fit)  # Estimate of phi and lambda
fitted(fit)
with(abdat, weighted.mean(y, w)) # Compare this with fitted(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

