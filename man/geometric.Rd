\name{geometric}
\alias{geometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Geometric Distribution }
\description{
  Maximum likelihood estimation for the geometric distribution.
}
\usage{
geometric(link = "logit", earg = list(), expected = TRUE, imethod = 1,
          iprob = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, earg}{ 
  Parameter link function and extra argument applied to the
  parameter \eqn{p}{prob}, which lies in the unit interval.
  See \code{\link{Links}} for more choices,
  and \code{earg} in \code{\link{Links}} for general information.

  }
  \item{expected}{ 
  Logical. 
  Fisher scoring is used if \code{expected = TRUE}, else Newton-Raphson.

  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} or \code{3} which
  specifies the initialization method for the probability.
  If failure to converge occurs try another value.

  }
  \item{iprob}{
  Optional initial value.
  See \code{\link{CommonVGAMffArguments}} for more details.

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
  If \eqn{Y} has a geometric distribution with parameter \eqn{p}{prob} then
  \eqn{Y+1} has a positive-geometric distribution with the same parameter.


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
  \code{\link{expgeometric}},
  \code{\link{zageometric}},
  \code{\link{zigeometric}},
  \code{\link{rbetageom}}.


}
\examples{
gdata = data.frame(x2 = runif(nn <- 1000) - 0.5)
gdata = transform(gdata, x3 = runif(nn) - 0.5,
                         x4 = runif(nn) - 0.5)
gdata = transform(gdata, eta = 1.0 - 1.0 * x2 + 2.0 * x3)
gdata = transform(gdata, prob = logit(eta, inverse = TRUE))
gdata = transform(gdata, y = rgeom(nn, prob))
with(gdata, table(y))
fit = vglm(y ~ x2 + x3 + x4, geometric, gdata, trace = TRUE)
coef(fit, matrix = TRUE)
summary(fit)
}
\keyword{models}
\keyword{regression}

