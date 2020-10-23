\name{pospoisson}
\alias{pospoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Poisson Distribution Family Function }
\description{
  Fits a positive Poisson distribution.
}
\usage{
pospoisson(link = "loglink", type.fitted = c("mean", "lambda", "prob0"),
    expected = TRUE, ilambda = NULL, imethod = 1, zero = NULL, gt.1 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function for the usual mean (lambda) parameter of
  an ordinary Poisson distribution.
  See \code{\link{Links}} for more choices.

  }
  \item{expected}{
  Logical.
  Fisher scoring is used if \code{expected = TRUE}, else Newton-Raphson.

  }
  \item{ilambda, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.

  }
  \item{type.fitted}{
    See \code{\link{CommonVGAMffArguments}} for details.


  }
  \item{gt.1}{
  Logical.
  Enforce \code{lambda > 1}? The default is to
  enforce \code{lambda > 0}.


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
Coleman, J. S. and James, J. (1961).
The equilibrium size distribution of freely-forming groups.
\emph{Sociometry}, \bold{24}, 36--45.


%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.


}
\author{ Thomas W. Yee }
\note{
  This family function can handle multiple responses.


  Yet to be done: a \code{quasi.pospoisson} which estimates a dispersion
  parameter.



}
\seealso{
  \code{\link[VGAM]{Gaitpois}},
  \code{\link{gaitpoisson}},
  \code{\link{posnegbinomial}},
  \code{\link{poissonff}},
  \code{\link{zapoisson}},
  \code{\link{zipoisson}},
  \code{\link{simulate.vlm}},
  \code{\link[VGAMdata]{otpospoisson}},
  \code{\link[VGAMdata]{Pospois}}.



}
\examples{
# Data from Coleman and James (1961)
cjdata <- data.frame(y = 1:6, freq = c(1486, 694, 195, 37, 10, 1))
fit <- vglm(y ~ 1, pospoisson, data = cjdata, weights = freq)
Coef(fit)
summary(fit)
fitted(fit)

pdata <- data.frame(x2 = runif(nn <- 1000))  # Artificial data
pdata <- transform(pdata, lambda = exp(1 - 2 * x2))
pdata <- transform(pdata, y1 = rgaitpois(nn, lambda, truncate = 0))
with(pdata, table(y1))
fit <- vglm(y1 ~ x2, pospoisson, data = pdata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}

% pdata <- transform(pdata, y1 = rpospois(nn, lambda))




