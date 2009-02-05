\name{riceff}
\alias{riceff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rice Distribution Family Function}
\description{
  Estimates the two parameters of a Rice distribution
  by maximum likelihood estimation.

}
\usage{
riceff(lvee="loge", lsigma="loge",
       evee=list(), esigma=list(),
       ivee=NULL, isigma=NULL,
       nsimEIM=100, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lvee, evee}{
  Link function and extra argument for the \eqn{v} parameter.
  See \code{\link{Links}} for more choices and for general information.

  }
  \item{lsigma, esigma}{
  Link function and extra argument for the \eqn{\sigma}{sigma} parameter.
  See \code{\link{Links}} for more choices and for general information.

  }
  \item{ivee, isigma}{
  Optional initial values for the parameters.
  See \code{\link{CommonVGAMffArguments}} for more information.
  If convergence failure occurs (this \pkg{VGAM} family function seems
  to require good initial values) try using these arguments.

  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
}
\details{
  The Rician distribution has density function
  \deqn{f(y;v,\sigma) =
  \frac{ y }{\sigma^2} \, \exp(-(y^2+v^2) / (2\sigma^2)) \, I_0(y v / \sigma^2)
  }{%
  f(y;v,sigma) = 
  (y/sigma^2) * exp(-(y^2+v^2) / (2*sigma^2)) * I_0(y*v/sigma^2)}
  where \eqn{y>0},
  \eqn{v > 0},
  \eqn{\sigma > 0} and \eqn{I_0} is the modified Bessel function of the
  first kind with order zero.
  When \eqn{v=0} the Rice distribution reduces to a Rayleigh distribution.
  The mean is
  \eqn{\sigma \sqrt{\pi/2} \exp(z/2) ((1-z) I_0(-z/2)-z I_1(-z/2))}{sigma*sqrt(pi/2)*exp(z/2)*((1-z)*I_0(-z/2)-z*I_1(-z/2))}
  (returned as the fitted values) where \eqn{z=-v^2/(2 \sigma^2)}{z=-v^2/(2*sigma^2)}.
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
  \code{\link[base:Bessel]{besselI}}.

}
\examples{
vee = exp(2); sigma = exp(1);
y = rrice(n <- 1000, vee, sigma)
fit = vglm(y ~ 1, riceff, trace=TRUE, crit="c")
c(mean(y), fitted(fit)[1])
coef(fit, matrix=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

