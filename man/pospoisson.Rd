\name{pospoisson}
\alias{pospoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Poisson Distribution Family Function }
\description{
  Fits a positive Poisson distribution.
}
\usage{
pospoisson(link = "loge")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function for the usual mean (lambda) parameter of
  an ordinary Poisson distribution.
  See \code{\link{Links}} for more choices.

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
y = 1:6
w = c(1486, 694, 195, 37, 10, 1) # Data from Coleman and James (1961)
fit = vglm(y ~ 1, pospoisson, weights=w)
Coef(fit)
summary(fit)
fitted(fit)

# Artificial data
x = runif(n <- 1000)
lambda = exp(2 + 3*x)
y = rpospois(n, lambda)
table(y)
fit = vglm(y ~ x, pospoisson, trace=TRUE, crit="c")
coef(fit, matrix=TRUE)
}
\keyword{models}
\keyword{regression}

