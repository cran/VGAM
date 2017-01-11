\name{riceff}
\alias{riceff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rice Distribution Family Function}
\description{
  Estimates the two parameters of a Rice distribution
  by maximum likelihood estimation.


}
\usage{
riceff(lsigma = "loge", lvee = "loge", isigma = NULL,
       ivee = NULL, nsimEIM = 100, zero = NULL, nowarning = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nowarning}{ Logical. Suppress a warning?
  Ignored for \pkg{VGAM} 0.9-7 and higher.


   }


  \item{lvee, lsigma}{
  Link functions for the \eqn{v} and \eqn{\sigma}{sigma} parameters.
  See \code{\link{Links}} for more choices and for general information.


  }
  \item{ivee, isigma}{
  Optional initial values for the parameters.
  If convergence failure occurs (this \pkg{VGAM} family function seems
  to require good initial values) try using these arguments.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  The Rician distribution has density function
  \deqn{f(y;v,\sigma) =
  \frac{ y }{\sigma^2} \, \exp(-(y^2+v^2) / (2\sigma^2)) \, I_0(y v / \sigma^2)
  }{%
  f(y;v,sigma) =
  (y/sigma^2) * exp(-(y^2+v^2) / (2*sigma^2)) * I_0(y*v/sigma^2)}
  where \eqn{y > 0},
  \eqn{v > 0},
  \eqn{\sigma > 0} and \eqn{I_0} is the modified Bessel function of the
  first kind with order zero.
  When \eqn{v = 0} the Rice distribution reduces to a Rayleigh distribution.
  The mean is
  \eqn{\sigma \sqrt{\pi/2} \exp(z/2)
       ((1-z) I_0(-z/2)-z I_1(-z/2))}{sigma*sqrt(pi/2) *
       exp(z/2)*((1-z) * I_0(-z/2)-z*I_1(-z/2))}
  (returned as the fitted values) where
  \eqn{z=-v^2/(2 \sigma^2)}{z=-v^2/(2*sigma^2)}.
  Simulated Fisher scoring is implemented.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

Rice, S. O. (1945)
Mathematical Analysis of Random Noise.
\emph{Bell System Technical Journal},
\bold{24}, 46--156.


}
\author{ T. W. Yee }
\note{
  Convergence problems may occur for data where \eqn{v=0}; if so, use
  \code{\link{rayleigh}} or possibly use an \code{\link{identity}} link.


  When \eqn{v} is large (greater than 3, say) then the mean is approximately
  \eqn{v} and the standard deviation is approximately \eqn{\sigma}{sigma}.


}

\seealso{
  \code{\link{drice}},
  \code{\link{rayleigh}},
  \code{\link[base:Bessel]{besselI}},
  \code{\link{simulate.vlm}}.


}
\examples{
\dontrun{ sigma <- exp(1); vee <- exp(2)
rdata <- data.frame(y = rrice(n <- 1000, sigma, vee = vee))
fit <- vglm(y ~ 1, riceff, data = rdata, trace = TRUE, crit = "coef")
c(with(rdata, mean(y)), fitted(fit)[1])
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
}
\keyword{models}
\keyword{regression}

