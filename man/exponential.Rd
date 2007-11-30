\name{exponential}
\alias{exponential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exponential Distribution }
\description{
  Maximum likelihood estimation for the exponential distribution.

}
\usage{
exponential(link = "loge", earg = list(), location = 0, expected = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function applied to the positive parameter \eqn{rate}.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{location}{
  Numeric of length 1, the known location parameter, \eqn{A}, say.

  }
  \item{expected}{
  Logical. If \code{TRUE} Fisher scoring is used,
  otherwise Newton-Raphson. The latter is usually faster.

  }

}
\details{
  The family function assumes the response \eqn{Y} has density
  \deqn{f(y) = \lambda \exp(-\lambda (y-A))}{%
        f(y) = rate * exp(-rate * (y-A)) }
  for \eqn{y > A}, where \eqn{A} is the known location parameter.
  By default, \eqn{A=0}.
  Then \eqn{E(Y) = A + 1/ \lambda}{E(Y) = A + 1/rate} and
  \eqn{Var(Y) = 1/ \lambda^2}{Var(Y) = 1/rate^2}.
}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.
}

\author{ T. W. Yee }
\note{ 
  Suppose \eqn{A=0}.
  For a fixed time interval, the number of events is 
  Poisson with mean \eqn{\lambda}{rate} if the time
  between events has a
  geometric distribution with mean \eqn{\lambda^{-1}}{1/rate}.
  The argument \code{rate} in \code{exponential} is the same as
  \code{\link[stats:Exponential]{rexp}} etc.
  The argument \code{lambda} in \code{\link{rpois}} is somewhat
  the same as \code{rate} here.

}
\seealso{
%   \code{\link{cexpon}},
    \code{\link{poissonff}},
    \code{\link{mix2exp}},
    \code{\link{freund61}}.
}

\examples{
nn = 100
x1 = runif(nn) - 0.5
x2 = runif(nn) - 0.5
eta = 0.2 - 0.7 * x1 + 1.9 * x2
rate = exp(eta)
y = rexp(nn, rate=rate)
stem(y)
fit = vglm(y ~ x1 + x2, exponential, trace=TRUE, crit="c") # slower
fit = vglm(y ~ x1 + x2, exponential(exp=FALSE), trace=TRUE, crit="c") # faster
coef(fit)
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

