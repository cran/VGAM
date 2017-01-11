\name{waldff}
\alias{waldff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Wald Distribution Family Function }
\description{
Estimates the parameter of the standard Wald distribution
by maximum likelihood estimation.


}
\usage{
waldff(llambda = "loge", ilambda = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llambda,ilambda}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  The standard Wald distribution is a special case of the
  inverse Gaussian distribution with \eqn{\mu=1}{mu=1}.
  It has a density that can be written as
  \deqn{f(y;\lambda) = \sqrt{\lambda/(2\pi y^3)}
       \; \exp\left(-\lambda (y-1)^2/(2 y)\right)}{%
       f(y;mu,lambda) = sqrt(lambda/(2*pi*y^3)) *
       exp(-lambda*(y-1)^2/(2*y))
 }
 where \eqn{y>0} and \eqn{\lambda>0}{lambda>0}.
 The mean of \eqn{Y} is \eqn{1}
 (returned as the fitted values) and its variance is
 \eqn{1/\lambda}{1/lambda}.
 By default, \eqn{\eta=\log(\lambda)}{eta=log(lambda)}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
\emph{Continuous Univariate Distributions},
2nd edition,
Volume 1,
New York: Wiley.


}
\author{ T. W. Yee }
\note{
  The \pkg{VGAM} family function \code{\link{inv.gaussianff}}
  estimates the location parameter \eqn{\mu}{mu} too.


}


\seealso{
  \code{\link{inv.gaussianff}},
  \code{\link{rinv.gaussian}}.


}
\examples{
wdata <- data.frame(y = rinv.gaussian(n = 1000, mu =  1, lambda = exp(1)))
wfit <- vglm(y ~ 1, waldff(ilambda = 0.2), data = wdata, trace = TRUE)
coef(wfit, matrix = TRUE)
Coef(wfit)
summary(wfit)
}
\keyword{models}
\keyword{regression}

