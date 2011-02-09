\name{pospoisson}
\alias{pospoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Poisson Distribution Family Function }
\description{
  Fits a positive Poisson distribution.
}
\usage{
pospoisson(link = "loge", earg = list(),
           expected = TRUE, ilambda = NULL, method.init = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function for the usual mean (lambda) parameter of
  an ordinary Poisson distribution.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{expected}{
  Logical.
  Fisher scoring is used if \code{expected = TRUE}, else Newton-Raphson.

  }
  \item{ilambda}{
  Optional initial value for \eqn{\lambda}{lambda}.
  A \code{NULL} means a value is computed internally.

  }
  \item{method.init}{
  An integer with value \code{1} or \code{2} or \code{3} which 
  specifies the initialization method for \eqn{\lambda}{lambda}.
  If failure to converge occurs try another value
  and/or else specify a value for \code{ilambda}.

  }

}
\details{
  The positive Poisson  distribution is the ordinary Poisson
  distribution but with the probability of zero being zero.  Thus the
  other probabilities are scaled up (i.e., divided by \eqn{1-P[Y=0]}).
  The mean, \eqn{\lambda / (1 - \exp(-\lambda))}{lambda/(1-exp(-lambda))},
  can be obtained by the extractor function \code{fitted} applied to
  the object.

  A related distribution is the zero-inflated Poisson, in which the
  probability \eqn{P[Y=0]} involves another parameter \eqn{\phi}{phi}.
  See \code{\link{zipoisson}}.

}
\section{Warning }{
    Under- or over-flow may occur if the data is ill-conditioned.

}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.

}
\references{
Coleman, J. S. and James, J. (1961)
The equilibrium size distribution of freely-forming groups.
\emph{Sociometry}, \bold{24}, 36--45.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
  Yet to be done: a \code{quasi.pospoisson} which estimates a dispersion
  parameter.

  This family function can handle a multivariate response.

}
\seealso{ 
\code{\link{Pospois}},
\code{\link{posnegbinomial}},
\code{\link{poissonff}},
\code{\link{zipoisson}}.

}
\examples{
# Data from Coleman and James (1961)
cjdat = data.frame(y = 1:6, w = c(1486, 694, 195, 37, 10, 1))
fit = vglm(y ~ 1, pospoisson, cjdat, weights=w)
Coef(fit)
summary(fit)
fitted(fit)

# Artificial data
pdat = data.frame(x = runif(nn <- 1000))
pdat = transform(pdat, lambda = exp(1 - 2*x))
pdat = transform(pdat, y = rpospois(nn, lambda))
with(pdat, table(y))
fit = vglm(y ~ x, pospoisson, pdat, trace=TRUE, crit="c")
coef(fit, matrix=TRUE)
}
\keyword{models}
\keyword{regression}

