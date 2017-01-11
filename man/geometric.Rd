\name{geometric}
\alias{geometric}
\alias{truncgeometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Geometric (Truncated and Untruncated) Distributions }
\description{
  Maximum likelihood estimation for the geometric
  and truncated geometric distributions.


}
\usage{
geometric(link = "logit", expected = TRUE, imethod = 1,
          iprob = NULL, zero = NULL)
truncgeometric(upper.limit = Inf,
               link = "logit", expected = TRUE, imethod = 1,
               iprob = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function applied to the
  probability parameter \eqn{p}{prob}, which lies in the unit interval.
  See \code{\link{Links}} for more choices.


  }
  \item{expected}{
  Logical.
  Fisher scoring is used if \code{expected = TRUE}, else Newton-Raphson.


  }
  \item{iprob, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for details.


  }

  \item{upper.limit}{
  Numeric.
  Upper values.
  As a vector, it is recycled across responses first.
  The default value means both family functions should give the same result.



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
  The geometric distribution is also a special case of the
  Borel distribution, which is a Lagrangian distribution.
  If \eqn{Y} has a geometric distribution with parameter \eqn{p}{prob} then
  \eqn{Y+1} has a positive-geometric distribution with the same parameter.
  Multiple responses are permitted.


  For \code{truncgeometric()},
  the (upper) truncated geometric distribution can have response integer
  values from 0 to \code{upper.limit}.
  It has density \code{prob * (1 - prob)^y / [1-(1-prob)^(1+upper.limit)]}.


  For a generalized truncated geometric distribution with
  integer values \eqn{L} to \eqn{U}, say, subtract \eqn{L}
  from the response and feed in \eqn{U-L} as the upper limit.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}

\author{
T. W. Yee.
Help from Viet Hoang Quoc is gratefully acknowledged.



}
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
  \code{\link{rbetageom}},
  \code{\link{simulate.vlm}}.


}
\examples{
gdata <- data.frame(x2 = runif(nn <- 1000) - 0.5)
gdata <- transform(gdata, x3 = runif(nn) - 0.5,
                          x4 = runif(nn) - 0.5)
gdata <- transform(gdata, eta  = -1.0 - 1.0 * x2 + 2.0 * x3)
gdata <- transform(gdata, prob = logit(eta, inverse = TRUE))
gdata <- transform(gdata, y1 = rgeom(nn, prob))
with(gdata, table(y1))
fit1 <- vglm(y1 ~ x2 + x3 + x4, geometric, data = gdata, trace = TRUE)
coef(fit1, matrix = TRUE)
summary(fit1)

# Truncated geometric (between 0 and upper.limit)
upper.limit <- 5
tdata <- subset(gdata, y1 <= upper.limit)
nrow(tdata)  # Less than nn
fit2 <- vglm(y1 ~ x2 + x3 + x4, truncgeometric(upper.limit),
             data = tdata, trace = TRUE)
coef(fit2, matrix = TRUE)

# Generalized truncated geometric (between lower.limit and upper.limit)
lower.limit <- 1
upper.limit <- 8
gtdata <- subset(gdata, lower.limit <= y1 & y1 <= upper.limit)
with(gtdata, table(y1))
nrow(gtdata)  # Less than nn
fit3 <- vglm(y1 - lower.limit ~ x2 + x3 + x4,
             truncgeometric(upper.limit - lower.limit),
             data = gtdata, trace = TRUE)
coef(fit3, matrix = TRUE)
}
\keyword{models}
\keyword{regression}

