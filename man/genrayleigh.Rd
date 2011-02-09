\name{genrayleigh}
\alias{genrayleigh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized Rayleigh Distribution Family Function}
\description{
  Estimates the two parameters of the generalized Rayleigh distribution
  by maximum likelihood estimation.

}
\usage{
genrayleigh(lshape = "loge", lscale = "loge",
            eshape = list(), escale = list(),
            ishape = NULL,   iscale = NULL,
            tol12 = 1e-05, nsimEIM = 300, zero = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, lscale}{
  Link function for the two positive parameters, shape and scale.
  See \code{\link{Links}} for more choices.

  }
  \item{eshape, escale}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ishape, iscale}{
  Numeric.
  Optional initial values for the shape and scale parameters.

  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}}.

}
  \item{tol12}{
  Numeric and positive. Tolerance for testing whether the
  second shape parameter is either 1 or 2. If so then the
  working weights need to handle these singularities.

}
}
\details{
  The generalized Rayleigh distribution has density function
  \deqn{f(y;a = shape,b = scale)  =
  (2 a y/b^{2}) e^{-(y/b)^{2}} (1 - e^{-(y/b)^{2}})^{a-1}}{%
  (2*a*y/b^2) * e^(-(y/b)^2) * (1 - e^(-(y/b)^2))^(a-1)}
  where \eqn{y > 0} and the two parameters,
  \eqn{a} and \eqn{b}, are positive.
  The mean cannot be expressed nicely so the median is returned as 
  the fitted values.
  Applications of the generalized Rayleigh distribution include modeling
  strength data and general lifetime data.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{
  Kundu, D., Raqab, M. C. (2005).
  Generalized Rayleigh distribution: different 
  methods of estimations.
  \emph{Computational Statistics and Data Analysis},
  \bold{49}, 187--200.

}
\author{ J. G. Lauder and T. W. Yee }
\note{
  We define \code{scale} as the reciprocal of the scale parameter 
  used by Kundu and Raqab (2005).

}

\seealso{ 
  \code{\link{dgenray}},
  \code{\link{rayleigh}}.

}
\examples{
shape = exp(1); scale = exp(2); 
rdata = data.frame(y = rgenray(n = 1000, shape, scale))
fit = vglm(y ~ 1, genrayleigh, rdata, trace = TRUE)
c(with(rdata, mean(y)), head(fitted(fit),1))
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

