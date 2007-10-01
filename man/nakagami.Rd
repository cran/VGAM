\name{nakagami}
\alias{nakagami}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Nakagami Distribution Family Function }
\description{
  Estimation of the two parameters of the
  Nakagami distribution by maximum likelihood estimation.

}
\usage{
nakagami(lshape = "loge", lscale = "loge",
         eshape=list(), escale=list(), ishape = NULL, iscale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, lscale}{
  Parameter link function applied to the
  \emph{shape} and \emph{scale} parameters.
  Log links ensure they are positive.
  See \code{\link{Links}} for more choices.

  }
  \item{eshape, escale}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ishape, iscale}{
  Optional initial values for the shape and scale parameters.
  For \code{ishape}, a \code{NULL} value means it is obtained in the
  \code{initialize} slot based on the value of \code{iscale}.
  For \code{iscale}, assigning a \code{NULL} means a value is obtained in the
  \code{initialize} slot, however, setting another numerical
  value is recommended if convergence fails or is too slow.

  }
}
\details{
  The Nakagami distribution, which is useful for modelling
  wireless systems such as radio links, can be written
  \deqn{f(y) = 2 (shape/scale)^{shape} y^{2 \times shape-1}
        \exp(-shape \times y^2/scale) / \Gamma(shape)}{%
	2 * (shape/scale)^shape *
        y^(2*shape-1) * exp(-shape*y^2/scale) /
        gamma(shape)}
  for \eqn{y > 0}, \eqn{shape > 0}, \eqn{scale > 0}.
  The mean of \eqn{Y} is
  \eqn{\sqrt{scale/shape}  \times  \Gamma(shape+0.5) /
       \Gamma(shape)}{sqrt(scale/shape) * gamma(shape+0.5) / gamma(shape)} and
  these are returned as the fitted values.
  By default, the linear/additive predictors are
  \eqn{\eta_1=\log(shape)}{eta1=log(shape)} and
  \eqn{\eta_2=\log(scale)}{eta2=log(scale)}.  
  Fisher scoring is implemented.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{
  Nakagami, M. (1960)
  The  \emph{m}-distribution: a general  formula  of
  intensity  distribution  of  rapid  fading,
  pp.3--36 in:
  \emph{Statistical Methods in Radio Wave Propagation}.
  W. C. Hoffman, Ed., New York: Pergamon.

}
\author{ T. W. Yee }
\note{
  The Nakagami distribution is also known as the
  Nakagami-\emph{m} distribution, where \eqn{m=shape} here.
  Special cases: \eqn{m=0.5} is a one-sided Gaussian
  distribution and \eqn{m=1} is a Rayleigh distribution.
  The second moment is \eqn{E(Y^2)=m}.

  If \eqn{Y} has a Nakagami distribution with parameters \emph{shape}
  and \emph{scale} then \eqn{Y^2} has a gamma distribution with shape
  parameter \emph{shape} and scale parameter \emph{scale/shape}.

}

\seealso{
  \code{\link{rnaka}},
  \code{\link{gamma2}},
  \code{\link{rayleigh}}.
}
\examples{
n = 1000; shape = exp(0); Scale = exp(1)
y = sqrt(rgamma(n, shape=shape, scale=Scale/shape))
fit = vglm(y ~ 1, nakagami, trace=TRUE, crit="c")
y = rnaka(n, shape=shape, scale=Scale)
fit = vglm(y ~ 1, nakagami(iscale=3), trace=TRUE)
fitted(fit)[1:5]
mean(y)
coef(fit, matrix=TRUE)
(Cfit = Coef(fit))
\dontrun{
hist(sy <- sort(y), prob=TRUE, main="", xlab="y", ylim=c(0,0.6))
lines(sy, dnaka(sy, shape=Cfit[1], scale=Cfit[2]), col="red")
}
}
\keyword{models}
\keyword{regression}

