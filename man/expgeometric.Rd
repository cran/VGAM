\name{expgeometric}
\alias{expgeometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exponential Geometric Distribution Family Function}
\description{
  Estimates the two parameters of the exponential geometric distribution
  by maximum likelihood estimation.

}
\usage{
expgeometric(lscale = "loge", lshape = "logit",
             iscale = NULL,   ishape = NULL,
             tol12 = 1e-05, zero = 1, nsimEIM = 400)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, lshape}{
  Link function for the two parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{iscale, ishape}{
  Numeric.
  Optional initial values for the scale and shape parameters.

  }
  \item{tol12}{
  Numeric.
  Tolerance for testing whether a parameter has value 1 or 2.

  }
  \item{zero, nsimEIM}{
  See \code{\link{CommonVGAMffArguments}}.

  }

}
\details{
  The exponential geometric distribution has density function
  \deqn{f(y; c = scale, s = shape)  = 
   (1/c) (1 - s) e^{-y/c} (1 - s e^{-y/c})^{-2}}{%
   (1/c) * (1 - s) * e^(-y/c) * (1 - s * e^(-y/c))^(-2)}
  where \eqn{y > 0}, \eqn{c > 0} and \eqn{s \in (0, 1)}{0 < s < 1}.
  The mean, \eqn{(c (s - 1)/ s) \log(1 - s)}{(c *
    (s - 1)/ s) * log(1 - s)}
  is returned as the fitted values.
  Note the median is \eqn{c \log(2 - s)}{c * log(2 - s)}.
  Simulated Fisher scoring is implemented.


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
\author{ J. G. Lauder and T. W. Yee }
\note{
  We define \code{scale} as the reciprocal of the scale parameter
  used by Adamidis and Loukas (1998).


}

\seealso{ 
  \code{\link{dexpgeom}},
  \code{\link{exponential}},
  \code{\link{geometric}}.


}
\examples{
\dontrun{
Scale <- exp(2); shape = logit(-1, inverse = TRUE);
edata <- data.frame(y = rexpgeom(n = 2000, scale = Scale, shape = shape))
fit <- vglm(y ~ 1, expgeometric, edata, trace = TRUE)
c(with(edata, mean(y)), head(fitted(fit), 1))
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
}
\keyword{models}
\keyword{regression}

