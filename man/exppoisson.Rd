\name{exppoisson}
\alias{exppoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exponential Poisson Distribution Family Function}
\description{
  Estimates the two parameters of the exponential Poisson distribution
  by maximum likelihood estimation.

}
\usage{
exppoisson(lrate = "loge", lshape = "loge",
           irate = 2, ishape = 1.1,  zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, lrate}{
  Link function for the two positive parameters.
  See \code{\link{Links}} for more choices.


  }
  \item{ishape, irate}{
  Numeric.
  Initial values for the \code{shape} and \code{rate} parameters.
  Currently this function is not intelligent enough to
  obtain better initial values.


  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The exponential Poisson distribution has density function
  \deqn{f(y; \beta = rate, \lambda = shape)  =
  \frac{\lambda \beta}{1 - e^{-\lambda}} \,
  e^{-\lambda - \beta y + \lambda \exp{(-\beta y)}}}{%
  f(y; a = shape, b = rate)  =
  (a*b/(1 - e^(-a))) * e^{-a - b*y + a * e^(-b*y)}}
  where \eqn{y > 0},
  and the parameters shape, \eqn{\lambda}{a},
  and rate, \eqn{\beta}{b}, are positive.
  The distribution implies a population facing discrete
  hazard rates which are multiples of a base hazard.
  This \pkg{VGAM} family function requires the \code{hypergeo} package
  (to use their \code{genhypergeo} function).
  The median is returned as the fitted value.






% This \pkg{VGAM} family function requires the \pkg{hypergeo} package
% (to use their \code{\link[hypergeo]{genhypergeo}} function).



}
\section{Warning }{
  This \pkg{VGAM} family function does not work properly!


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
  Kus, C., (2007).
  A new lifetime distribution.
  \emph{Computational Statistics and Data Analysis},
  \bold{51}, 4497--4509.


}
\author{ J. G. Lauder, jamesglauder@gmail.com }


\seealso{ 
  \code{\link{dexppois}},
  \code{\link{exponential}},
  \code{\link{poisson}}.


}
\examples{
\dontrun{
shape <- exp(1); rate <- exp(2)
rdata <- data.frame(y = rexppois(n = 1000, rate = rate, shape = shape))
library("hypergeo")  # Required!
fit <- vglm(y ~ 1, exppoisson, data = rdata, trace = FALSE, maxit = 1200)
c(with(rdata, median(y)), head(fitted(fit), 1))
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
}
\keyword{models}
\keyword{regression}

