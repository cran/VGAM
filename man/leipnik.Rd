\name{leipnik}
\alias{leipnik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Leipnik Distribution Family Function}
\description{
  Estimates the two parameters of a (transformed) Leipnik distribution
  by maximum likelihood estimation.

}
\usage{
leipnik(lmu = "logit", llambda = "loge",
        emu=list(), elambda=list(),
        imu = NULL, ilambda = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, llambda}{
  Link function for the \eqn{\mu}{mu} and \eqn{\lambda}{lambda} parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{imu, ilambda}{
  Numeric. Optional initial values for \eqn{\mu}{mu} and
  \eqn{\lambda}{lambda}.

  }
  \item{emu, elambda}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
}
\details{
  The (transformed) Leipnik distribution has density function
  \deqn{f(y;\mu,\lambda) = \frac{ \{ y(1-y) \}^{-\frac12}}{
  \mbox{Beta}( \frac{\lambda+1}{2}, \frac12 )}
  \left[ 1 +  \frac{(y-\mu)^2 }{y(1-y)}
  \right]^{ -\frac{\lambda}{2}}}{%
  f(y;mu,lambda) = 
  (y(1-y))^(-1/2) * (1 + (y-mu)^2 / (y*(1-y)))^(-lambda/2) / 
  Beta((lambda+1)/2, 1/2)}
  where \eqn{0 < y < 1} and \eqn{\lambda > -1}{lambda > -1}.
  The mean is \eqn{\mu}{mu} (returned as the fitted values) and the variance is
  \eqn{1/\lambda}{1/lambda}.

  Jorgensen (1997) calls the above the \bold{transformed} Leipnik
  distribution, and if \eqn{y = (x+1)/2} and \eqn{\mu = (\theta+1)/2}{mu
  = (theta+1)/2}, then the distribution of \eqn{X} as a function of
  \eqn{x} and \eqn{\theta}{theta} is known as the the (untransformed)
  Leipnik distribution.  Here, both \eqn{x} and \eqn{\theta}{theta}
  are in \eqn{(-1,1)}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.

}
\references{
  Jorgensen, B. (1997)
  \emph{The Theory of Dispersion Models}.
  London: Chapman & Hall

  Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1995)
  \emph{Continuous Univariate Distributions},
  2nd edition,
  Volume 2,
  New York: Wiley.
  (pages 612--617).

}
\author{ T. W. Yee }
\note{
  Convergence may be slow or fail.
  Until better initial value estimates are forthcoming try assigning the
  argument \code{ilambda} some numerical value if it fails to converge.
  Currently, Newton-Raphson is implemented, not Fisher scoring.
  Currently, this family function probably only really works for
  intercept-only models, i.e., \code{y ~ 1} in the formula.

}

\section{Warning }{
  If \code{llambda="identity"} then it is possible that the
  \code{lambda} estimate becomes less than \eqn{-1}, i.e., out of
  bounds. One way to stop this is to choose \code{llambda="loge"},
  however, \code{lambda} is then constrained to be positive.

}

\seealso{ 
    \code{\link{mccullagh89}}.
}
\examples{
y = rnorm(n=2000, mean=0.5, sd=0.1)  # Not good data
fit = vglm(y ~ 1, leipnik(ilambda=1), tr=TRUE, checkwz=FALSE)
fit = vglm(y ~ 1, leipnik(ilambda=1,llam=logoff, elam=list(offset=1)),
           trace=TRUE, cri="coef")
fitted(fit)[1:5]
mean(y)
summary(fit)
coef(fit, matrix=TRUE)
Coef(fit)

sum(weights(fit))  # sum of the prior weights
sum(weights(fit, type="w")) # sum of the working weights
}
\keyword{models}
\keyword{regression}

%fit = vglm(y ~ 1, leipnik(ilambda=1), tr=TRUE, cri="c", checkwz=FALSE)