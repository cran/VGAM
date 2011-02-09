\name{expgeometric}
\alias{expgeometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exponential geometric Distribution Family Function}
\description{
  Estimates the two parameters of the exponential geometric distribution
  by maximum likelihood estimation.

}
\usage{
expgeometric(lscale = "loge", lshape = "logit",
             escale = list(), eshape = list(),
             iscale = NULL,   ishape = NULL,
             zero = 1, nsimEIM = 400)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, lshape}{
  Link function for the two parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{escale, eshape}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{iscale, ishape}{
  Numeric.
  Optional initial values for the scale and shape parameters.

  }
  \item{zero, nsimEIM}{
  See \code{\link{CommonVGAMffArguments}}.

  }


}
\details{
  The exponential geometric distribution has density function
  \deqn{f(y; c = scale, s = shape)  = 
   c (1 - s) e^{-c y} (1 - s e^{-c y})^{-2}}{%
   c*(1-s)*e^(-c*y)*(1-s*e^(-c*y))^(-2)}
  where \eqn{y > 0}, \eqn{c > 0} and \eqn{s \in (0, 1)}{0 < s < 1}.
  The mean is \eqn{((s-1)/(sc)) \log(1-s)}{((s-1)/(sc)) * log(1-s)}
  is returned as the fitted values.
  Note the median is \eqn{(1/c) \log(2-s)}{(1/c) * log(2-s)}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{
  Adamidis, K., Loukas, S. (1998).
  A lifetime distribution with decreasing failure rate.
  \emph{Statistics and Probability Letters},
  \bold{39}, 35--42.

}
\author{ J. G. Lauder and T. W .Yee }
\note{
  Soon the scale term will be redefinined as the reciprocal of
  the scale parameter used by Adamidis and Loukas (1998).
  That is, division by the scale parameter will be used later.


}

\seealso{ 
  \code{dexpgeom},
  \code{\link{exponential}},
  \code{\link{geometric}}.

}
\examples{
scale = exp(2); shape = logit(-1, inverse = TRUE);
edata = data.frame(y = rexpgeom(n = 2000, scale=scale, shape=shape))
fit = vglm(y ~ 1, expgeometric, edata, trace = TRUE)
c(with(edata, mean(y)), head(fitted(fit), 1))
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

