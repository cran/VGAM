\name{gammahyp}
\alias{gammahyp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gamma Hyperbola Bivariate Distribution }
\description{
  Estimate the parameter of a gamma hyperbola bivariate distribution
  by maximum likelihood estimation.
}
\usage{
gammahyp(ltheta="loge", itheta=NULL, expected=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ltheta}{
  Link function applied to the (positive) parameter \eqn{\theta}{theta}.
  See \code{\link{Links}} for more choices.

  }
  \item{itheta}{
  Initial value for the parameter.
  The default is to estimate it internally.

  }
  \item{expected}{
  Logical. \code{FALSE} means the Newton-Raphson (using
  the observed information matrix) algorithm, otherwise the expected
  information matrix is used (Fisher scoring algorithm).

  }
}
\details{
  The joint probability density function is given by
  \deqn{f(y_1,y_2) = \exp( -e^{-\theta} y_1 / \theta - \theta y_2 )}{%
        f(y1,y2) = exp( -exp(-theta) * y1 / theta - theta * y2) }
  for \eqn{\theta > 0}{theta > 0}, \eqn{y_1 > 0}{y1 > 0}, \eqn{y_2 > 1}{y2 > 1}.
  The random variables \eqn{Y_1}{Y1}  and \eqn{Y_2}{Y2} are independent.
  The marginal distribution of \eqn{Y_1}{Y1} is an exponential distribution
  with rate parameter \eqn{\exp(-\theta)/\theta}{exp(-theta)/theta}.
  The marginal distribution of \eqn{Y_2}{Y2} is an exponential distribution
  that has been shifted to the right by 1 and with
  rate parameter \eqn{\theta}{theta}.
  The fitted values are stored in a two-column matrix with the marginal
  means, which are \eqn{\theta \exp(\theta)}{theta * exp(theta)} and
  \eqn{1 + 1/\theta}{1 + 1/theta}.

  The default algorithm is Newton-Raphson because Fisher scoring tends to
  be much slower for this distribution.
}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.
}

\references{

Reid, N. (2003)
Asymptotics and the theory of inference. 
\emph{Annals of Statistics},
\bold{31}, 1695--1731.

}
\author{ T. W. Yee }
\note{
  The response must be a two column matrix.

}
\seealso{
  \code{\link{exponential}}.
}
\examples{
x  = runif(n <- 1000)
theta = exp(-2+x)
y1 = rexp(n, rate=exp(-theta)/theta)
y2 = 1 + rexp(n, rate=theta)
fit = vglm(cbind(y1,y2) ~ x, fam=gammahyp(expected=TRUE), trace=TRUE)
fit = vglm(cbind(y1,y2) ~ x, fam=gammahyp, trace=TRUE, crit="coef")
coef(fit, matrix=TRUE)
Coef(fit)
fitted(fit)[1:4,]
summary(fit)
}
\keyword{models}
\keyword{regression}

