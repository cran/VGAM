\name{geometric}
\alias{geometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Geometric Distribution }
\description{
  Maximum likelihood estimation for the geometric distribution.
}
\usage{
geometric(link = "logit", earg=list(), expected = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{ 
  Parameter link function applied to the
  parameter \eqn{p}{prob}, which lies in the unit interval.
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
}
\details{
  A random variable \eqn{Y} has a 1-parameter geometric distribution
  if \eqn{P(Y=y) = p (1-p)^y}{P(Y=y) = prob * (1-prob)^y}
  for \eqn{y=0,1,2,\ldots}{y=0,1,2,...}.
  Here, \eqn{p}{prob} is the probability of success,
  and \eqn{Y} is the number of (independent) trials that are fails
  until a success occurs.
  Thus the response \eqn{Y} should be a non-negative integer.
  The mean of \eqn{Y} is \eqn{E(Y) = (1-p)/p}{E(Y) = (1-prob)/prob}
  and its variance is \eqn{Var(Y) = (1-p)/p^2}{Var(Y) = (1-prob)/prob^2}.
  The geometric distribution is a special case of the 
  negative binomial distribution (see \code{\link{negbinomial}}).

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
%\note{
%
%}

\seealso{ 
  \code{\link{negbinomial}},
  \code{\link[stats]{Geometric}},
  \code{\link{betageometric}},
  \code{\link{rbetageom}}.
}
\examples{
x1 = runif(n <- 1000) - 0.5
x2 = runif(n) - 0.5
x3 = runif(n) - 0.5
eta = 0.2 - 0.7 * x1 + 1.9 * x2
prob = logit(eta, inverse=TRUE)
y = rgeom(n, prob)
table(y)
fit = vglm(y ~ x1 + x2 + x3, geometric, trace=TRUE, crit="coef")
coef(fit)
coef(fit, mat=TRUE)
summary(fit)
}
\keyword{models}
\keyword{regression}
