\name{nakagami}
\alias{nakagami}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Nakagami Distribution Family Function }
\description{
  Estimation of the two parameters of the
  Nakagami distribution by maximum likelihood estimation.

}
\usage{
nakagami(lscale = "loge", lshape = "loge", iscale = 1, ishape = NULL,
         nowarning = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nowarning}{ Logical. Suppress a warning? }


  \item{lscale, lshape}{
  Parameter link functions applied to the
  \emph{scale} and \emph{shape} parameters.
  Log links ensure they are positive.
  See \code{\link{Links}} for more choices
  and information.


  }
  \item{iscale, ishape}{
  Optional initial values for the shape and scale parameters.
  For \code{ishape}, a \code{NULL} value means it is obtained in the
  \code{initialize} slot based on the value of \code{iscale}.
  For \code{iscale}, assigning a \code{NULL} means a value
  is obtained in the \code{initialize} slot, however, setting
  another numerical value is recommended if convergence fails or
  is too slow.


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
       \Gamma(shape)}{sqrt(scale/shape) * gamma(shape+0.5) / gamma(shape)}
  and these are returned as the fitted values.
  By default, the linear/additive predictors are
  \eqn{\eta_1=\log(scale)}{eta1=log(scale)} and
  \eqn{\eta_2=\log(shape)}{eta2=log(shape)}.
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


  If \eqn{Y} has a Nakagami distribution with parameters
  \emph{shape} and \emph{scale} then \eqn{Y^2} has a gamma
  distribution with shape parameter \emph{shape} and scale
  parameter \emph{scale/shape}.


}

\seealso{
  \code{\link{rnaka}},
  \code{\link{gamma2}},
  \code{\link{rayleigh}}.


}
\examples{
nn <- 1000; shape <- exp(0); Scale <- exp(1)
ndata <- data.frame(y1 = sqrt(rgamma(nn, shape = shape, scale = Scale/shape)))
fit <- vglm(y1 ~ 1, nakagami, data = ndata, trace = TRUE, crit = "coef")
ndata <- transform(ndata, y2 = rnaka(nn, scale = Scale, shape = shape))
fit <- vglm(y2 ~ 1, nakagami(iscale = 3), data = ndata, trace = TRUE)
head(fitted(fit))
with(ndata, mean(y2))
coef(fit, matrix = TRUE)
(Cfit <- Coef(fit))
\dontrun{ sy <- with(ndata, sort(y2))
hist(with(ndata, y2), prob = TRUE, main = "", xlab = "y", ylim = c(0, 0.6),
     col = "lightblue")
lines(dnaka(sy, scale = Cfit["scale"], shape = Cfit["shape"]) ~ sy,
      data = ndata, col = "orange") }
}
\keyword{models}
\keyword{regression}

