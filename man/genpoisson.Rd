\name{genpoisson}
\alias{genpoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Poisson distribution }
\description{
  Estimation of the two parameters of a generalized Poisson distribution.
}
\usage{
genpoisson(llambda = "logit", ltheta = "loge",
           elambda=list(), etheta=list(),
           ilambda = 0.5, itheta = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llambda}{
  Parameter link function for \code{lambda}.
  See \code{\link{Links}} for more choices.
  The parameter lies in the unit interval, therefore the
  logit link is a reasonable default.

  }
  \item{ltheta}{
  Parameter link function for \code{theta}. 
  See \code{\link{Links}} for more choices.
  The parameter is positive, therefore the default is the log link. 

  }
  \item{elambda, etheta}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ilambda}{ Optional initial value for \code{lambda}. }
  \item{itheta}{ Optional initial value for \code{theta}. }
  \item{zero}{ An integer vector, containing the value 1 or 2. If so,
  \code{lambda} or \code{theta} respectively are modelled as an intercept only.
  The default value \code{NULL} means both linear/additive predictors
  are modelled as functions of the explanatory variables.
}
}
\details{
  See Consul (1989) for details.
}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{
Jorgensen, B. (1997)
\emph{The Theory of Dispersion Models}.
London: Chapman & Hall

Consul, P. C. (1989)
\emph{Generalized Poisson Distributions: Properties and Applications},
Marcel Dekker.

}
\author{ T. W. Yee }
\note{
  Convergence problems may occur when \code{lambda} is very close to 0. 
  This distribution is useful for dispersion modelling.
}
\seealso{
  \code{\link{poissonff}}.
}

\examples{
y = rpois(n=100, lam=24)
fit  = vglm(y ~ 1, genpoisson, trace=TRUE)
fitted(fit)[1:5]
mean(y)
summary(fit)
coef(fit, matrix=TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}
