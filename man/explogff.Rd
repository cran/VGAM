\name{explogff}
\alias{explogff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exponential Logarithmic Distribution Family Function}
\description{
  Estimates the two parameters of the exponential logarithmic distribution
  by maximum likelihood estimation.

}
\usage{
explogff(lscale = "loge", lshape = "logit",
         iscale = NULL,   ishape = NULL,
         tol12 = 1e-05, zero = 1, nsimEIM = 400)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, lshape}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{tol12}{
  Numeric.
  Tolerance for testing whether a parameter has value 1 or 2.


  }
  \item{iscale, ishape, zero, nsimEIM}{
  See \code{\link{CommonVGAMffArguments}}.


  }

}
\details{
  The exponential logarithmic distribution has density function
  \deqn{f(y; c, s)  = 
   (1/(-\log p ))  (((1/c)   (1 - s)   e^{-y/c}) / (1 - (1 - s)   e^{-y/c}))}{%
   (1/(-log(p))) * (((1/c) * (1 - s) * e^(-y/c)) / (1 - (1 - s) * e^(-y/c)))}
  where \eqn{y > 0}, scale parameter \eqn{c > 0}, and
  shape parameter \eqn{s \in (0, 1)}{0 < s < 1}.
  The mean, \eqn{(-polylog(2, 1 - p) c) / \log(s)}{((-polylog(2, 1 - p) *
   c) / log(s)} is \emph{not} returned as the fitted values.
  Note the median is \eqn{c \log(1 + \sqrt{s})}{c * log(1 + sqrt(s))}
  and it is \emph{currently} returned as the fitted values.
  Simulated Fisher scoring is implemented.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
  Tahmasabi, R., Sadegh, R. (2008).
  A two-parameter lifetime distribution with decreasing failure rate.
  \emph{Computational Statistics and Data Analysis},
  \bold{52}, 3889--3901.
 
 
}
\author{ J. G. Lauder and T. W .Yee }
\note{
  We define \code{scale} as the reciprocal of the rate parameter
  used by Tahmasabi and Sadegh (2008).


  Yet to do: find a \code{polylog()} function.


}

\seealso{ 
  \code{\link{dexplog}},
  \code{\link{exponential}},


}
\examples{
\dontrun{ Scale <- exp(2); shape <- logit(-1, inverse = TRUE)
edata <- data.frame(y = rexplog(n = 2000, scale = Scale, shape = shape))
fit <- vglm(y ~ 1, explogff, data = edata, trace = TRUE)
c(with(edata, median(y)), head(fitted(fit), 1))
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
}
\keyword{models}
\keyword{regression}

