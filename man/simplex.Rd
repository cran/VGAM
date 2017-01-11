\name{simplex}
\alias{simplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simplex Distribution Family Function }
\description{
  The two parameters of the univariate standard simplex
  distribution are estimated by full maximum likelihood
  estimation.

}
\usage{
simplex(lmu = "logit", lsigma = "loge", imu = NULL, isigma = NULL,
        imethod = 1, ishrinkage = 0.95, zero = "sigma")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, lsigma}{
  Link function for \code{mu} and \code{sigma}.
  See \code{\link{Links}} for more choices.


  }
  \item{imu, isigma}{
  Optional initial values for \code{mu} and \code{sigma}.
  A \code{NULL} means a value is obtained internally.


  }
  \item{imethod, ishrinkage, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  The probability density function can be written
  \deqn{f(y; \mu, \sigma) = [2 \pi \sigma^2 (y (1-y))^3]^{-0.5}
  \exp[-0.5 (y-\mu)^2 / (\sigma^2 y (1-y) \mu^2 (1-\mu)^2)]
  }{%
  f(y; mu, sigma) = [2* pi * sigma^2 * (y*(1-y))^3]^(-0.5) *
   exp[-0.5 * (y-mu)^2 / (sigma^2 * y * (1-y) * mu^2 * (1-mu)^2)] }
  for \eqn{0 < y < 1},
      \eqn{0 < \mu < 1}{0 < mu < 1},
  and \eqn{\sigma > 0}{sigma > 0}.
  The mean of \eqn{Y} is \eqn{\mu}{mu} (called \code{mu}, and
  returned as the fitted values).



% This comes from Jorgensen but it is not confirmed by simulations:
% The variance of \eqn{Y} is \eqn{\mu (1 - \mu) - \sqrt{ \lambda / 2}
% \exp\{ \lambda / (\mu^2 (1 - \mu)^2) \}
% \Gamma(\lambda / (2 \mu^2 (1 - \mu)^2), 0.5)}{
% mu * (1 - mu) - sqrt(lambda / 2) *
% exp(lambda / (mu^2 * (1 - mu)^2)) *
% Gamma(lambda / (2 * mu^2 * (1 - mu)^2), 0.5)}.
% Here, \eqn{\Gamma(x, a)}{Gamma(x, a)} is the
% `upper' normalized incomplete gamma function given by
% \code{pgamma(x, a, lower = FALSE) * gamma(a)}.


  The second parameter, \code{sigma}, of this standard simplex
  distribution is known as the dispersion parameter.
  The unit variance function is
  \eqn{V(\mu) = \mu^3 (1-\mu)^3}{V(mu) = mu^3 (1-mu)^3}.
  Fisher scoring is applied to both parameters.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Jorgensen, B. (1997)
  \emph{The Theory of Dispersion Models}.
  London: Chapman & Hall


  Song, P. X.-K. (2007)
  \emph{Correlated Data Analysis: Modeling, Analytics, and Applications}.
  Springer.


}
\author{ T. W. Yee }
\note{
  This distribution is potentially useful for dispersion modelling.
  Numerical problems may occur when \code{mu} is very close to 0 or 1.


}

\seealso{
  \code{\link{dsimplex}},
  \code{\link{dirichlet}},
  \code{\link{rig}},
  \code{\link{binomialff}}.


}
\examples{
sdata <- data.frame(x2 = runif(nn <- 1000))
sdata <- transform(sdata, eta1 = 1 + 2 * x2,
                          eta2 = 1 - 2 * x2)
sdata <- transform(sdata, y = rsimplex(nn, mu = logit(eta1, inverse = TRUE),
                                       dispersion = exp(eta2)))
(fit <- vglm(y ~ x2, simplex(zero = NULL), data = sdata, trace = TRUE))
coef(fit, matrix = TRUE)
summary(fit)
}
\keyword{models}
\keyword{regression}






