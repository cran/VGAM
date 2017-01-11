\name{skellam}
\alias{skellam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Skellam Distribution Family Function}
\description{
  Estimates the two parameters of a Skellam distribution
  by maximum likelihood estimation.

}
\usage{
skellam(lmu1 = "loge", lmu2 = "loge", imu1 = NULL, imu2 = NULL,
        nsimEIM = 100, parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu1, lmu2}{
  Link functions for the \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2} parameters.
  See \code{\link{Links}} for more choices and for general information.


  }
  \item{imu1, imu2}{
  Optional initial values for the parameters.
  See \code{\link{CommonVGAMffArguments}} for more information.
  If convergence failure occurs (this \pkg{VGAM} family function seems
  to require good initial values) try using these arguments.


  }
  \item{nsimEIM, parallel, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  In particular, setting \code{parallel=TRUE} will constrain the
  two means to be equal.


  }
}
\details{
  The Skellam distribution models the difference between two
  independent Poisson distributions
  (with means \eqn{\mu_{j}}{mu_j}, say).
  It has density function
  \deqn{f(y;\mu_1,\mu_2) =
  \left( \frac{ \mu_1 }{\mu_2} \right)^{y/2} \,
  \exp(-\mu_1-\mu_2 ) \, I_{|y|}( 2 \sqrt{ \mu_1 \mu_2})
  }{%
f(y;mu1,mu2) =
  ( \mu1 / mu_2 )^(y/2) *
  exp(-mu1-mu2 ) * I_(|y|)( 2 * sqrt(mu1*mu2))
  }
  where \eqn{y} is an integer,
  \eqn{\mu_1 > 0}{mu1 > 0},
  \eqn{\mu_2 > 0}{mu2 > 0}.
  Here, \eqn{I_v} is the modified Bessel function of the
  first kind with order \eqn{v}.


  The mean is \eqn{\mu_1 - \mu_2}{mu1 - mu2}
  (returned as the fitted values),
  and the variance is \eqn{\mu_1 + \mu_2}{mu1 + mu2}.
  Simulated Fisher scoring is implemented.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\section{Warning }{
  This \pkg{VGAM} family function seems fragile and very sensitive to
  the initial values.
  Use very cautiously!!


}

\references{

Skellam, J. G. (1946)
The frequency distribution of the difference between
two Poisson variates belonging to different populations.
\emph{Journal of the Royal Statistical Society, Series A},
\bold{109}, 296.


}
%\author{ T. W. Yee }
\note{
  Numerical problems may occur for data if \eqn{\mu_1}{mu1} and/or
  \eqn{\mu_2}{mu2} are large.


}

\seealso{
  \code{\link{dskellam}},
  \code{\link[stats:Poisson]{dpois}},
  \code{\link{poissonff}}.


}
\examples{
\dontrun{
sdata <- data.frame(x2 = runif(nn <- 1000))
sdata <- transform(sdata, mu1 = exp(1 + x2), mu2 = exp(1 + x2))
sdata <- transform(sdata, y = rskellam(nn, mu1, mu2))
fit1 <- vglm(y ~ x2, skellam, data = sdata, trace = TRUE, crit = "coef")
fit2 <- vglm(y ~ x2, skellam(parallel = TRUE), data = sdata, trace = TRUE)
coef(fit1, matrix = TRUE)
coef(fit2, matrix = TRUE)
summary(fit1)
# Likelihood ratio test for equal means:
pchisq(2 * (logLik(fit1) - logLik(fit2)),
       df = df.residual(fit2) - df.residual(fit1), lower.tail = FALSE)
lrtest(fit1, fit2)  # Alternative
}
}
\keyword{models}
\keyword{regression}

