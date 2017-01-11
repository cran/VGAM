\name{cardioid}
\alias{cardioid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cardioid Distribution Family Function }
\description{
  Estimates the two parameters of the
  cardioid distribution by maximum likelihood estimation.
}
\usage{
cardioid(lmu = extlogit(min = 0, max = 2*pi),
         lrho = extlogit(min = -0.5, max = 0.5),
         imu = NULL, irho = 0.3, nsimEIM = 100, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, lrho}{
  Parameter link functions applied to the \eqn{\mu}{mu}
  and \eqn{\rho}{rho} parameters, respectively.
  See \code{\link{Links}} for more choices.


  }
  \item{imu, irho}{
  Initial values.
  A \code{NULL} means an initial value is chosen internally.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }

}
\details{
  The two-parameter cardioid distribution
  has a density that can be written as
  \deqn{f(y;\mu,\rho) = \frac{1}{2\pi}
        \left(1 + 2\, \rho \cos(y - \mu) \right) }{%
        f(y;mu,rho) = (1 + 2*rho*cos(y-mu)) / (2*pi)}
  where \eqn{0 < y < 2\pi}{0 < y < 2*pi},
  \eqn{0 < \mu < 2\pi}{0 < mu < 2*pi}, and
  \eqn{-0.5 < \rho < 0.5}{-0.5 < rho < 0.5} is the concentration parameter.
 The default link functions enforce the range constraints
 of the parameters.


  For positive \eqn{\rho} the distribution is unimodal and symmetric about
  \eqn{\mu}{mu}.
  The mean of \eqn{Y} (which make up the fitted values) is
  \eqn{\pi + (\rho/\pi) ((2 \pi-\mu) \sin(2 \pi-\mu) +
       \cos(2 \pi-\mu) - \mu \sin(\mu) - \cos(\mu))}{
        pi + (rho/pi) ((2*pi-mu)*sin(2*pi-mu) +
       cos(2*pi-mu) - mu*sin(mu) - cos(mu))}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{

Jammalamadaka, S. R. and SenGupta, A. (2001)
\emph{Topics in Circular Statistics},
Singapore: World Scientific.

}
\author{ T. W. Yee }
\note{
 Fisher scoring using simulation is used.


}
\section{Warning }{
  Numerically, this distribution can be difficult to fit because of a
  log-likelihood having multiple maximums.
  The user is therefore encouraged to try different starting values,
  i.e., make use of \code{imu} and \code{irho}.


}

\seealso{
  \code{\link{rcard}},
  \code{\link{extlogit}},
  \code{\link{vonmises}}.


  \pkg{CircStats} and \pkg{circular} currently have a lot more
  R functions for circular data than the \pkg{VGAM} package.


}
\examples{
\dontrun{
cdata <- data.frame(y = rcard(n = 1000, mu = 4, rho = 0.45))
fit <- vglm(y ~ 1, cardioid, data = cdata, trace = TRUE)
coef(fit, matrix=TRUE)
Coef(fit)
c(with(cdata, mean(y)), head(fitted(fit), 1))
summary(fit)
}
}
\keyword{models}
\keyword{regression}

