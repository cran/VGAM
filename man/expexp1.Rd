\name{expexp1}
\alias{expexp1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exponentiated Exponential Distribution }
\description{
  Estimates the two parameters of the exponentiated exponential
  distribution by maximizing a profile (concentrated) likelihood.

}
\usage{
expexp1(lscale = "loge", escale = list(), iscale = NULL, ishape = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale}{
  Parameter link function for the (positive) \eqn{\lambda}{scale} parameter.
  See \code{\link{Links}} for more choices.

  }
  \item{escale}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{iscale}{
  Initial value for the \eqn{\lambda}{scale} parameter.
  By default, an initial value is chosen internally using \code{ishape}.

  }
  \item{ishape}{
  Initial value for the \eqn{\alpha}{shape} parameter. If convergence
  fails try setting a different value for this argument.

  }
}
\details{
  See \code{\link{expexp}} for details about the exponentiated
  exponential distribution. This family function uses a different
  algorithm for fitting the model. Given \eqn{\lambda}{scale},
  the MLE of \eqn{\alpha}{shape} can easily be solved in terms of
  \eqn{\lambda}{scale}. This family function maximizes a profile
  (concentrated) likelihood with respect to \eqn{\lambda}{scale}.
  Newton-Raphson is used, which compares with Fisher scoring with
  \code{\link{expexp}}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.
}
\references{
  Gupta, R. D. and Kundu, D. (2001)
  Exponentiated exponential family: an alternative to
  gamma and Weibull distributions,
  \emph{Biometrical Journal},
  \bold{43},
  117--130.

}

\author{ T. W. Yee }
\note{
  This family function works only for intercept-only models,
  i.e., \code{y ~ 1} where \code{y} is the response.

  The estimate of  \eqn{\alpha}{shape} is attached to the
  \code{misc} slot of the object, which is a list and contains
  the component \code{shape}.
  
  As Newton-Raphson is used, the working weights are sometimes
  negative, and some adjustment is made to these to make them
  positive.

  Like \code{\link{expexp}}, good initial
  values are needed. Convergence may be slow.

}

\section{Warning }{The standard errors produced by a
  \code{summary} of the model may be wrong.

}

\seealso{
  \code{\link{expexp}},
  \code{\link{CommonVGAMffArguments}}.

}
\examples{
# Ball bearings data (number of million revolutions before failure)
bbearings = data.frame(y = c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60,
48.80, 51.84, 51.96, 54.12, 55.56, 67.80, 68.64, 68.64,
68.88, 84.12, 93.12, 98.64, 105.12, 105.84, 127.92,
128.04, 173.40))
fit = vglm(y ~ 1, expexp1(ishape = 4), bbearings, trace = TRUE,
           maxit = 50, checkwz = FALSE)
coef(fit, matrix = TRUE)
Coef(fit) # Authors get c(0.0314, 5.2589) with log-lik -112.9763
fit@misc$shape    # Estimate of shape
logLik(fit)


# Failure times of the airconditioning system of an airplane
acplane = data.frame(y = c(23, 261, 87, 7, 120, 14, 62, 47,
225, 71, 246, 21, 42, 20, 5, 12, 120, 11, 3, 14,
71, 11, 14, 11, 16, 90, 1, 16, 52, 95))
fit = vglm(y ~ 1, expexp1(ishape = 0.8), acplane, trace = TRUE,
           maxit = 50, checkwz = FALSE)
coef(fit, matrix = TRUE)
Coef(fit) # Authors get c(0.0145, 0.8130) with log-lik -152.264
fit@misc$shape    # Estimate of shape
logLik(fit)
}
\keyword{models}
\keyword{regression}

