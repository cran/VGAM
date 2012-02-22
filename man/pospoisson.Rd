\name{pospoisson}
\alias{pospoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Poisson Distribution Family Function }
\description{
  Fits a positive Poisson distribution.
}
\usage{
pospoisson(link = "loge", earg = list(),
           expected = TRUE, ilambda = NULL, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, earg}{
  Link function and extra argument for the usual mean (lambda) parameter of
  an ordinary Poisson distribution.
  See \code{\link{Links}} for more choices.

  }
  \item{expected}{
  Logical.
  Fisher scoring is used if \code{expected = TRUE}, else Newton-Raphson.

  }
  \item{ilambda, imethod}{
  See \code{\link{CommonVGAMffArguments}} for more information.

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
  This family function can handle a multivariate response.

  Yet to be done: a \code{quasi.pospoisson} which estimates a dispersion
  parameter.


}
\seealso{ 
  \code{\link{Pospois}},
  \code{\link{posnegbinomial}},
  \code{\link{poissonff}},
  \code{\link{zipoisson}}.

}
\examples{
# Data from Coleman and James (1961)
cjdat = data.frame(y = 1:6, freq = c(1486, 694, 195, 37, 10, 1))
fit = vglm(y ~ 1, pospoisson, cjdat, weights = freq)
Coef(fit)
summary(fit)
fitted(fit)

pdat = data.frame(x2 = runif(nn <- 1000)) # Artificial data
pdat = transform(pdat, lambda = exp(1 - 2 * x2))
pdat = transform(pdat, y1 = rpospois(nn, lambda))
with(pdat, table(y1))
fit = vglm(y1 ~ x2, pospoisson, pdat, trace = TRUE, crit = "coef")
coef(fit, matrix=TRUE)
}
\keyword{models}
\keyword{regression}

