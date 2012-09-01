\name{hypersecant}
\alias{hypersecant}
\alias{hypersecant.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hyperbolic Secant Distribution Family Function }
\description{
  Estimation of the parameter of the hyperbolic secant
  distribution.

}
\usage{
hypersecant(link.theta = elogit(min = -pi/2, max = pi/2), init.theta = NULL)
hypersecant.1(link.theta = elogit(min = -pi/2, max = pi/2), init.theta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.theta}{
  Parameter link function applied to the parameter \eqn{\theta}{theta}.
  See \code{\link{Links}} for more choices.

  }
  \item{init.theta}{
  Optional initial value for \eqn{\theta}{theta}.
  If failure to converge occurs, try some other value.
  The default means an initial value is determined internally.

  }
}
\details{
  The probability density function of the hyperbolic secant distribution
  is given by
  \deqn{f(y)=\exp(\theta y + \log(\cos(\theta ))) / (2 \cosh(\pi y/2)),}{%
        f(y) =exp(theta*y + log(cos(theta ))) / (2*cosh(pi*y/2)),}
  for parameter \eqn{-\pi/2 < \theta < \pi/2}{pi/2 < theta < pi/2}
  and all real \eqn{y}.
  The mean of \eqn{Y} is \eqn{\tan(\theta)}{tan(theta)} (returned as
  the fitted values).


  Another parameterization is used for \code{hypersecant.1()}.
  This uses
  \deqn{f(y)=(\cos(\theta)/\pi) \times y^{-0.5+\theta/\pi} \times
             (1-y)^{-0.5-\theta/\pi},}{%
        f(y) =(cos(theta)/pi) * y^(-0.5+theta/pi) * (1-y)^(-0.5-theta/pi),}
  for parameter \eqn{-\pi/2 < \theta < \pi/2}{pi/2 < theta < pi/2}
  and \eqn{0 < y < 1}.
  Then the mean of \eqn{Y} is \eqn{0.5 + \theta/\pi}{0.5 + theta/pi}
  (returned as the fitted values) and the variance is
  \eqn{(\pi^2 - 4 \theta^2) / (8\pi^2)}{(pi^2 - 4*theta^2) / (8*pi^2)}.


  For both parameterizations Newton-Raphson is same as Fisher scoring.


}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Jorgensen, B. (1997)
  \emph{The Theory of Dispersion Models}.
  London: Chapman & Hall.
%  p.101, Eqn (3.37).


}

\author{ T. W. Yee }
%\note{ 

%}
\seealso{
  \code{\link{elogit}}.


}
\examples{
hdata <- data.frame(x2 = rnorm(nn <- 200))
hdata <- transform(hdata, y = rnorm(nn))  # Not very good data!
fit <- vglm(y ~ x2, hypersecant, hdata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
fit@misc$earg

# Not recommended:
fit <- vglm(y ~ x2, hypersecant(link = "identity"), hdata, trace = TRUE)
coef(fit, matrix = TRUE)
fit@misc$earg
}
\keyword{models}
\keyword{regression}
