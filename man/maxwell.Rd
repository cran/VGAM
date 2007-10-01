\name{maxwell}
\alias{maxwell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maxwell Distribution Family Function }
\description{
 Estimating the parameter of the Maxwell  distribution by
 maximum likelihood estimation.
}
\usage{
maxwell(link = "loge", earg=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function applied to the parameter \eqn{a}. 
  See \code{\link{Links}} for more choices.
  A log link is the default because the parameter is positive.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
}
\details{
  The Maxwell distribution, which is used in the area of
  thermodynamics,
  has a probability density function that can be written
  \deqn{f(y) = \sqrt{2/\pi} a^{3/2} y^2 \exp(-0.5 a y^2)}{%
    f(y) = sqrt(2/pi) * a^(3/2) * y^2 * exp(-0.5*a*y^2)}
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
  A related distribution is the Rayleigh distribution.
  Fisher-scoring and Newton-Raphson are the same here.

}

\seealso{ 
  \code{\link{Maxwell}},
  \code{\link{rayleigh}}.
}
\examples{
y = rmaxwell(n <- 1000, a=exp(2))
fit = vglm(y ~ 1, maxwell, trace=TRUE, crit="c")
coef(fit, matrix=TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}
