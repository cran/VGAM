\name{zipoisson}
\alias{zipoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Poisson Distribution Family Function }
\description{
  Fits a zero-inflated Poisson distribution using full maximum likelihood
  estimation.

}
\usage{
zipoisson(lphi="logit", llambda = "loge", 
          ephi=list(), elambda =list(),
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
  Optional initial value for \eqn{\phi}{phi}, whose value must lie
  between 0 and 1.  The default is to compute an initial value internally.

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
  This function uses Fisher scoring and is based on
  \deqn{P(Y=0) =  \phi + (1-\phi) \exp(-\lambda),}{%
        P(Y=0) =  phi + (1-phi) * exp(-lambda),}
  and for \eqn{y=1,2,\ldots},
  \deqn{P(Y=y) =  (1-\phi) \exp(-\lambda) \lambda^y / y!.}{%
        P(Y=y) =  (1-phi) * exp(-lambda) * lambda^y / y!.}
  The parameter \eqn{\phi}{phi} satisfies \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is \eqn{(1-\phi) \lambda}{(1-phi)*lambda} and these
  are returned as the fitted values.  By default, the two linear/additive
  predictors are \eqn{(logit(\phi), \log(\lambda))^T}{(logit(phi),
  log(lambda))^T}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
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
  Numerical problems can occur.
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
  \code{\link[stats:Poisson]{rpois}}.
}
\examples{
x = runif(n <- 2000)
phi = logit(-0.5 + 1*x, inverse=TRUE)
lambda = loge(0.5 + 2*x, inverse=TRUE)
y = rzipois(n, lambda, phi)
table(y)
fit = vglm(y ~ x, zipoisson, trace=TRUE)
coef(fit, matrix=TRUE)  # These should agree with the above values


# Another example: data from McKendrick (1926).
y = 0:4  # Number of cholera cases per household in an Indian village
w = c(168, 32, 16, 6, 1)  # Frequencies; there are 223=sum(w) households
fit = vglm(y ~ 1, zipoisson, wei=w, trace=TRUE)
coef(fit, matrix=TRUE)
cbind(actual=w, fitted=
      dzipois(y, lambda=Coef(fit)[2], phi=Coef(fit)[1]) * sum(w))


# Another example: data from Angers and Biswas (2003)
y = 0:7
w = c(182, 41, 12, 2, 2, 0, 0, 1)
y = y[w>0]
w = w[w>0]
fit = vglm(y ~ 1, zipoisson(lphi=probit, iphi=0.3), wei=w, tra=TRUE)
fit@misc$prob0  # Estimate of P(Y=0)
coef(fit, matrix=TRUE)
Coef(fit)  # Estimate of phi and lambda
fitted(fit)
weighted.mean(y,w) # Compare this with fitted(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

