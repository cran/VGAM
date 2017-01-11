\name{betaff}
\alias{betaff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Two-parameter Beta Distribution Family Function }
\description{
  Estimation of the mean and precision parameters of the beta distribution.

}
\usage{
betaff(A = 0, B = 1, lmu = "logit", lphi = "loge",
       imu = NULL, iphi = NULL,
       gprobs.y = ppoints(8), gphi  = exp(-3:5)/4, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{


  \item{A, B}{
  Lower and upper limits of the distribution.
  The defaults correspond to the \emph{standard beta distribution}
  where the response lies between 0 and 1.
  }


  \item{lmu, lphi}{
  Link function for the mean and precision parameters.
  The values \eqn{A} and \eqn{B} are extracted from the
  \code{min} and \code{max} arguments of \code{\link{extlogit}}.
  Consequently, only \code{\link{extlogit}} is allowed.


% See below for more details.
% See \code{\link{Links}} for more choices.


  }
  \item{imu, iphi}{
  Optional initial value for the mean and precision parameters
  respectively. A \code{NULL} value means a value is obtained in the
  \code{initialize} slot.


  }
  \item{gprobs.y, gphi, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The two-parameter beta distribution can be written
  \eqn{f(y) =}
    \deqn{(y-A)^{\mu_1 \phi-1} \times (B-y)^{(1-\mu_1) \phi-1} / [beta(\mu_1
          \phi,(1-\mu_1) \phi) \times (B-A)^{\phi-1}]}{%
          (y-A)^(mu1*phi-1) * (B-y)^((1-mu1)*phi-1) / [beta(mu1*phi,(1-mu1)*phi) *
                 (B-A)^(phi-1)]}
  for \eqn{A < y < B}, and \eqn{beta(.,.)} is the beta function
  (see \code{\link[base:Special]{beta}}).
  The parameter \eqn{\mu_1}{mu1} satisfies
  \eqn{\mu_1 = (\mu - A) / (B-A)}{mu1 = (mu - A) / (B-A)}
  where \eqn{\mu}{mu} is the mean of \eqn{Y}.
  That is, \eqn{\mu_1}{mu1} is the mean of of a standard beta distribution:
  \eqn{E(Y) = A + (B-A) \times \mu_1}{E(Y) = A + (B-A)*mu1},
  and these are the fitted values of the object.
  Also, \eqn{\phi}{phi} is positive and \eqn{A < \mu < B}{A < mu < B}.
  Here, the limits \eqn{A} and \eqn{B} are \emph{known}.


  Another parameterization of the beta distribution involving the raw
  shape parameters is implemented in \code{\link{betaR}}.


  For general \eqn{A} and \eqn{B}, the variance of \eqn{Y} is
  \eqn{(B-A)^2 \times \mu_1 \times (1-\mu_1) / (1+\phi)}{(B-A)^2 *
       mu1 * (1-mu1) / (1+phi)}.
  Then \eqn{\phi}{phi} can be interpreted as a \emph{precision} parameter
  in the sense that, for fixed \eqn{\mu}{mu}, the larger the value of
  \eqn{\phi}{phi}, the smaller the variance of \eqn{Y}.
  Also, \eqn{\mu_1 = shape1/(shape1+shape2)}{mu1=shape1/(shape1+shape2)} and
  \eqn{\phi = shape1+shape2}{phi = shape1+shape2}.
  Fisher scoring is implemented.


% If \eqn{A} and \eqn{B} are unknown then the \pkg{VGAM} family function
% \code{beta4()} can be used to estimate these too.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Ferrari, S. L. P. and Francisco C.-N. (2004)
  Beta regression for modelling rates and proportions.
  \emph{Journal of Applied Statistics},
  \bold{31}, 799--815.


%  Documentation accompanying the \pkg{VGAM} package at
%  \url{https://www.stat.auckland.ac.nz/~yee}
%  contains further information and examples.


}
\author{ Thomas W. Yee }
\note{
  The response must have values in the interval (\eqn{A}, \eqn{B}).
  The user currently needs to manually choose \code{lmu} to match
  the input of arguments \code{A} and \code{B}, e.g.,
  with \code{\link{extlogit}}; see the example below.


}

\seealso{
  \code{\link{betaR}},
% \code{\link{zoibetaR}},
  \code{\link[stats:Beta]{Beta}},
  \code{\link{dzoabeta}},
  \code{\link{genbetaII}},
  \code{\link{betaII}},
  \code{\link{betabinomialff}},
  \code{\link{betageometric}},
  \code{\link{betaprime}},
  \code{\link{rbetageom}},
  \code{\link{rbetanorm}},
  \code{\link{kumar}},
  \code{\link{extlogit}},
  \code{\link{simulate.vlm}}.


}
\examples{
bdata <- data.frame(y = rbeta(nn <- 1000, shape1 = exp(0), shape2 = exp(1)))
fit1 <- vglm(y ~ 1, betaff, data = bdata, trace = TRUE)
coef(fit1, matrix = TRUE)
Coef(fit1)  # Useful for intercept-only models

# General A and B, and with a covariate
bdata <- transform(bdata, x2 = runif(nn))
bdata <- transform(bdata, mu   = logit(0.5 - x2, inverse = TRUE),
                          prec =   exp(3.0 + x2))  # prec == phi
bdata <- transform(bdata, shape2 = prec * (1 - mu),
                         shape1 = mu * prec)
bdata <- transform(bdata,
                   y = rbeta(nn, shape1 = shape1, shape2 = shape2))
bdata <- transform(bdata, Y = 5 + 8 * y)  # From 5 to 13, not 0 to 1
fit <- vglm(Y ~ x2, data = bdata, trace = TRUE,
            betaff(A = 5, B = 13, lmu = extlogit(min = 5, max = 13)))
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}
%      imu = NULL, iphi = NULL, imethod = 1, zero = NULL)



