\name{inv.gaussianff}
\alias{inv.gaussianff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inverse Gaussian Distribution Family Function }
\description{
  Estimates the two parameters of the inverse Gaussian distribution
  by maximum likelihood estimation.


}
\usage{
inv.gaussianff(lmu = "loge", llambda = "loge",
               imethod = 1, ilambda = NULL,
               parallel = FALSE, ishrinkage = 0.99, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
%apply.parint = FALSE,
\arguments{
  \item{lmu, llambda}{
  Parameter link functions for the \eqn{\mu}{mu} and
  \eqn{\lambda}{lambda} parameters.
  See \code{\link{Links}} for more choices.


  }
  \item{ilambda, parallel}{
  See \code{\link{CommonVGAMffArguments}} for more information.
  If \code{parallel = TRUE} then the constraint is not applied to the
  intercept.


  }
  \item{imethod, ishrinkage, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  The standard (``canonical'') form of the
  inverse Gaussian distribution has a density
  that can be written as
  \deqn{f(y;\mu,\lambda) = \sqrt{\lambda/(2\pi y^3)}
       \exp\left(-\lambda (y-\mu)^2/(2 \mu^2 y)\right)}{%
       f(y;mu,lambda) = sqrt(lambda/(2*pi*y^3)) *
       exp(-lambda*(y-mu)^2/(2*mu^2*y))
  }
  where \eqn{y>0},
  \eqn{\mu>0}{mu>0}, and
  \eqn{\lambda>0}{lambda>0}.
  The mean of \eqn{Y} is \eqn{\mu}{mu} and its variance is
  \eqn{\mu^3/\lambda}{mu^3/lambda}.
  By default, \eqn{\eta_1=\log(\mu)}{eta1=log(mu)} and
  \eqn{\eta_2=\log(\lambda)}{eta2=log(lambda)}.
  The mean is returned as the fitted values.
  This \pkg{VGAM} family function can handle multiple
  responses (inputted as a matrix).



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{

Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
\emph{Continuous Univariate Distributions},
2nd edition, Volume 1, New York: Wiley.


Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}
\author{ T. W. Yee }
\note{
  The inverse Gaussian distribution can be fitted (to a certain extent)
  using the usual GLM framework involving a scale parameter. This family
  function is different from that approach in that it estimates both
  parameters by full maximum likelihood estimation.


}

\seealso{
  \code{\link{Inv.gaussian}},
  \code{\link{waldff}},
  \code{\link{bisa}}.


  The \R{} package \pkg{SuppDists} has several functions for evaluating
  the density, distribution function, quantile function and generating
  random numbers from the inverse Gaussian distribution.


}
\examples{
idata <- data.frame(x2 = runif(nn <- 1000))
idata <- transform(idata, mymu   = exp(2 + 1 * x2),
                          Lambda = exp(2 + 1 * x2))
idata <- transform(idata, y = rinv.gaussian(nn, mu = mymu, lambda = Lambda))
fit1 <-   vglm(y ~ x2, inv.gaussianff, data = idata, trace = TRUE)
rrig <- rrvglm(y ~ x2, inv.gaussianff, data = idata, trace = TRUE)
coef(fit1, matrix = TRUE)
coef(rrig, matrix = TRUE)
Coef(rrig)
summary(fit1)
}
\keyword{models}
\keyword{regression}

