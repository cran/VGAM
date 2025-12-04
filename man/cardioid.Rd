\name{cardioid}
\alias{cardioid}
\alias{cardioid2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cardioid Distribution Family Function }
\description{
  Estimates the two parameters of the
  cardioid distribution by maximum likelihood estimation.
}
\usage{
cardioid(lmu = "extlogitlink(min = 0, max = 2*pi)",
     lrho = "extlogitlink(min = -0.5, max = 0.5)",
     imu = NULL, irho = NULL,
     gmu  = ppoints(16) * 2 * pi,
     grho = ppoints(16) - 0.5, zero = NULL)
cardioid2(lmu = "extlogitlink(min = 0, max = 2*pi)",
     lrho2 = "logitlink", imu = NULL, irho2 = NULL,
     gmu   = ppoints(8) * 2 * pi,
     grho2 = ppoints(8), zero = "rho2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, lrho, lrho2}{
  Parameter link functions applied to the \eqn{\mu}{mu}
  and \eqn{\rho}{rho} and \eqn{\rho_2}{rho2}
  parameters, respectively.
  See \code{\link{Links}} for more choices.


  }
  \item{gmu, grho, grho2}{
  Grid of initial values.
  Because of a possible multimodal likelihood,
  it is a good idea to try a fine rectangular grid
  of values.
  The default may be too fine and made coarser when
  \eqn{n} is very large.
  Pewsey (2025) recommends that MLE be used
  for \eqn{n > 10} so a warning may occur if
  \eqn{n} is too small.
  See \code{\link{CommonVGAMffArguments}} for more
  information.


}

  \item{imu, irho, irho2}{
  Initial values.
  A \code{NULL} means an initial value is
  chosen internally.
  See \code{\link{CommonVGAMffArguments}} for
  information.


  }
  \item{zero}{
    See \code{\link{CommonVGAMffArguments}} for
    information.


  }

}
\details{
  The two-parameter cardioid distribution
  has a density that can be written as
  \deqn{f(y;\mu,\rho) = \frac{1}{2\pi}
        \left[1 + 2\, \rho \cos(y - \mu) \right] }{%
        f(y;mu,rho) = [1 + 2*rho*cos(y-mu)] / (2*pi)}
      where \eqn{0 \leq y < 2\pi}{0 <= y  < 2*pi},
      and \eqn{0 \leq \mu < 2\pi}{0 <= mu < 2*pi}
      is the location parameter (of the mode), and
  \eqn{-0.5 < \rho < 0.5}{-0.5 < rho < 0.5} is
  the concentration parameter.
  The default link functions enforce the range
  constraints of the parameters.


  For positive \eqn{\rho} the distribution is unimodal and
  symmetric about \eqn{\mu}{mu}.
  The mean of \eqn{Y} (which make up the fitted values) is
  \eqn{\pi + (\rho/\pi) ((2 \pi-\mu) \sin(2 \pi-\mu) +
       \cos(2 \pi-\mu) - \mu \sin(\mu) - \cos(\mu))}{
        pi + (rho/pi) ((2*pi-mu)*sin(2*pi-mu) +
       cos(2*pi-mu) - mu*sin(mu) - cos(mu))}.



  Pewsey (2025) considers another parameterization:
  since
  \eqn{f(y; \mu, -\rho) = f(y; \mu+\pi, \rho)}{
       f(y;  mu, - rho) = f(y; \mu+ pi,  rho)}
   then \eqn{\rho}{rho} could be restricted to
   \eqn{[0, 0.5]}.
   This suggests a variant parameterization:
   \eqn{\rho_2 = 2\rho}{rho2 = 2*rho}
   which lies in the unit interval and therefore
   \code{\link{logitlink}},
   \code{\link{probitlink}},
   \code{\link{clogloglink}}, etc. are options.
   This is implemented in
   \code{\link{dcard2}}, \code{\link{pcard2}},
   \code{\link{qcard2}}, \code{\link{rcard2}},
   \code{\link{cardioid2}} and is recommended
   above the original parameterization. Another
   reason this is better is that it 
   ensure identifiability because
   \code{dcard(theta, mu, -rho)} is the same as
   \code{dcard(theta, mu + pi, rho)}.
   

}
\value{
  An object of class \code{"vglmff"} (see
  \code{\link{vglmff-class}}).
  The object is used by modelling
  functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{


Pewsey, A. (2025).
On Jeffreys's cardioid distribution.
\emph{Computational Statistics and Data Analysis},
\bold{82}, in press.


Jammalamadaka, S. R. and SenGupta, A. (2001).
\emph{Topics in Circular Statistics},
Singapore: World Scientific.

}
\author{ T. W. Yee }
\note{
  Fisher scoring using simulation is no longer used.
  Instead, the EIM has been derived (Pewsey, 2025).


}
\section{Warning }{
  Numerically, this distribution can be difficult
  to fit because
  of a log-likelihood having multiple maximums.
  Indeed, Fisher scoring may not be suitable
  for estimation, especially for \eqn{n < 10}.
  Often the (global) solution occurs when
  \eqn{\rho}{rho}
  is on or near the boundary of the parameter space.
  The user is therefore encouraged to try different
  starting values, i.e.,
  make use of \code{imu} and \code{irho},
  even though \code{gmu} and
  \code{grho}/\code{grho2} should suffice.


% \eqn{\rho = 0.5}{rho = 0.5} or thereabouts, so



}

\seealso{
  \code{\link{Card}},
  \code{\link{extlogitlink}},
  \code{\link{vonmises}}.


  \pkg{CircStats}, \pkg{circular} and
  \pkg{Directional}
  currently have
  a lot more R functions for circular data than
  \pkg{VGAM}.


}
\examples{
\dontrun{
cdata <- data.frame(y = rcard2(1000, mu = 4, rho2 = 0.5))
fit <- vglm(y ~ 1, cardioid2, data = cdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
c(with(cdata, mean(y)), head(fitted(fit), 1))
summary(fit)
}
}
\keyword{models}
\keyword{regression}

