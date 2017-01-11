\name{maxwell}
\alias{maxwell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maxwell Distribution Family Function }
\description{
 Estimating the parameter of the Maxwell  distribution by
 maximum likelihood estimation.

}
\usage{
maxwell(link = "loge", zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, zero}{
  Parameter link function applied to \eqn{a},
  which is called the parameter \code{rate}.
  See \code{\link{Links}} for more choices and information;
  a log link is the default because the parameter is positive.
  More information is at \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The Maxwell distribution, which is used in the area of
  thermodynamics,
  has a probability density function that can be written
  \deqn{f(y;a) = \sqrt{2/\pi} a^{3/2} y^2 \exp(-0.5 a y^2)}{%
    f(y;a) = sqrt(2/pi) * a^(3/2) * y^2 * exp(-0.5*a*y^2)}
  for \eqn{y>0} and \eqn{a>0}.
  The mean of \eqn{Y} is
  \eqn{\sqrt{8 / (a \pi)}}{sqrt(8 / (a * pi))}
  (returned as the fitted values), and its variance is
  \eqn{(3\pi - 8)/(\pi a)}{(3*pi - 8)/(pi*a)}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{

  von Seggern, D. H. (1993)
  \emph{CRC Standard Curves and Surfaces},
  Boca Raton, FL.: CRC Press.


}
\author{ T. W. Yee }
\note{
  Fisher-scoring and Newton-Raphson are the same here.
  A related distribution is the Rayleigh distribution.
  This \pkg{VGAM} family function handles multiple responses.
  This \pkg{VGAM} family function can be mimicked by
  \code{poisson.points(ostatistic = 1.5, dimension = 2)}.


}

\seealso{
  \code{\link{Maxwell}},
  \code{\link{rayleigh}},
  \code{\link{poisson.points}}.


}
\examples{
mdata <- data.frame(y = rmaxwell(1000, rate = exp(2)))
fit <- vglm(y ~ 1, maxwell, data = mdata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}
