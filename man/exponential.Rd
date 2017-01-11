\name{exponential}
\alias{exponential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exponential Distribution }
\description{
  Maximum likelihood estimation for the exponential distribution.


}
\usage{
exponential(link = "loge", location = 0, expected = TRUE,
            ishrinkage = 0.95, parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function applied to the positive parameter \eqn{rate}.
  See \code{\link{Links}} for more choices.


  }
  \item{location}{
  Numeric of length 1, the known location parameter, \eqn{A}, say.


  }
  \item{expected}{
  Logical. If \code{TRUE} Fisher scoring is used,
  otherwise Newton-Raphson. The latter is usually faster.


  }
  \item{ishrinkage, parallel, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


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

Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.



}

\author{ T. W. Yee }
\note{
  Suppose \eqn{A = 0}.
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
    \code{\link{amlexponential}},
    \code{\link{gpd}},
    \code{\link{laplace}},
    \code{\link{expgeometric}},
    \code{\link{explogff}},
    \code{\link{poissonff}},
    \code{\link{mix2exp}},
    \code{\link{freund61}},
    \code{\link{simulate.vlm}},
    \code{\link[stats]{Exponential}}.


%   \code{\link{cens.exponential}},


}

\examples{
edata <- data.frame(x2 = runif(nn <- 100) - 0.5)
edata <- transform(edata, x3 = runif(nn) - 0.5)
edata <- transform(edata, eta = 0.2 - 0.7 * x2 + 1.9 * x3)
edata <- transform(edata, rate = exp(eta))
edata <- transform(edata, y = rexp(nn, rate = rate))
with(edata, stem(y))

fit.slow <- vglm(y ~ x2 + x3, exponential, data = edata, trace = TRUE)
fit.fast <- vglm(y ~ x2 + x3, exponential(exp = FALSE), data = edata,
                 trace = TRUE, crit = "coef")
coef(fit.slow, mat = TRUE)
summary(fit.slow)


# Compare results with a GPD. Has a threshold.
threshold <- 0.5
gdata <- data.frame(y1 = threshold + rexp(n = 3000, rate = exp(1.5)))

fit.exp <- vglm(y1 ~ 1, exponential(location = threshold), data = gdata)
coef(fit.exp, matrix = TRUE)
Coef(fit.exp)
logLik(fit.exp)

fit.gpd <- vglm(y1 ~ 1, gpd(threshold =  threshold), data = gdata)
coef(fit.gpd, matrix = TRUE)
Coef(fit.gpd)
logLik(fit.gpd)
}
\keyword{models}
\keyword{regression}

