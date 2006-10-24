\name{mix2normal1}
\alias{mix2normal1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mixture of Two Univariate Normal Distributions }
\description{
  Estimates the five parameters of a mixture of two univariate 
  normal distributions by maximum likelihood estimation.

}
\usage{
mix2normal1(lphi="logit", lmu="identity", lsd="loge",
            iphi=0.5, imu1=NULL, imu2=NULL, isd1=NULL, isd2=NULL,
            qmu=c(0.2, 0.8), esd=FALSE, zero=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi}{
    Link function for the parameter \eqn{\phi}{phi}.
    See below for more details.
    See \code{\link{Links}} for more choices.

  }
  \item{lmu}{
    Link function applied to each \eqn{\mu}{mu} parameter.
    See \code{\link{Links}} for more choices.

  }
  \item{lsd}{
    Link function applied to each \eqn{\sigma}{sd} parameter.
    See \code{\link{Links}} for more choices.

  }
  \item{iphi}{
  Initial value for \eqn{\phi}{phi}, whose value must lie
  between 0 and 1.

  }
  \item{imu1, imu2}{
    Optional initial value for \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2}.
    The default is to compute initial values internally using
    the argument \code{qmu}.

  }
  \item{isd1, isd2}{
    Optional initial value for \eqn{\sigma_1}{sd1} and \eqn{\sigma_2}{sd2}.
    The default is to compute initial values internally based on
    the argument \code{qmu}.

  }
  \item{qmu}{
    Vector with two values giving the probabilities relating to the sample
    quantiles for obtaining initial values for \eqn{\mu_1}{mu1}
    and \eqn{\mu_2}{mu2}.
    The two values are fed in as the \code{probs} argument into
    \code{\link[stats]{quantile}}.

  }
  \item{esd}{
    Logical indicating whether the two standard deviations should be 
    constrained to be equal. If set \code{TRUE}, the appropriate
    constraint matrices will be used.

  }
\item{zero}{
  An integer specifying which linear/additive predictor is modelled as
  intercepts only.  If given, the value or values must be from the
  set \eqn{\{1,2,\ldots,5\}}{1,2,...,5}. 
  The default is the first one only, meaning \eqn{\phi}{phi}
  is a single parameter even when there are explanatory variables.
  Set \code{zero=NULL} to model all linear/additive predictors as
  functions of the explanatory variables.

}
}
\details{
  The probability function can be loosely written as 
  \deqn{f(y) = \phi \, N(\mu_1,\sigma_1^2) + (1-\phi) \, N(\mu_2, \sigma_2^2)}{%
        f(y) = phi * N(mu1, sd1^2) + (1-phi) * N(mu2, sd2^2)}
  where \eqn{\phi}{phi} is the probability an observation belongs
  to the first group.
  The parameters \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2} are the means, and 
  \eqn{\sigma_1}{sd1} and \eqn{\sigma_2}{sd2} are the standard deviations.
  The parameter \eqn{\phi}{phi} satisfies \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is
  \eqn{\phi \mu_1 + (1-\phi) \mu_2}{phi*mu1 + (1-phi)*mu2}
  and this is returned as the fitted values.
  By default, the five linear/additive predictors are
  \eqn{(logit(\phi), \mu_1, \log(\sigma_1), \mu_2, \log(\sigma_2))^T}{(logit(phi),
  mu1, log(sd1), mu2, log(sd2))^T}.
  If \code{esd=TRUE} then \eqn{\sigma_1 = \sigma_2}{sd1=sd2} is enforced.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{

Everitt, B. S. and Hand, D. J. (1981)
\emph{Finite Mixture Distributions}.
London: Chapman & Hall.

}
\section{Warning }{
  Numerical problems can occur.
  Half-stepping is not uncommon.
  If failure to converge occurs, try obtaining better initial values,
  e.g., by using \code{iphi} and \code{qmu} etc.

  This function uses a quasi-Newton update for the working weight matrices
  (BFGS variant). It builds up approximations to the weight matrices,
  and currently the code is not fully tested.
  In particular, results based on the weight matrices (e.g., from
  \code{vcov} and \code{summary}) may be quite incorrect, especially when
  the arguments \code{weights} is used to input prior weights.

}

\author{ T. W. Yee }
\note{
  Fitting this model successfully to data can be difficult due to
  numerical problems and ill-conditioned data.  It pays to fit the
  model several times with different initial values, and check that
  the best fit looks reasonable. Plotting the results is recommended.
  This function works better as \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2}
  become more different.

  Convergence is often slow, especially when the two component
  distributions are not well separated. The control argument \code{maxit}
  should be set to a higher value, e.g., 200, and use \code{trace=TRUE}
  to monitor convergence.  If appropriate in the first place, setting
  \code{esd=TRUE} often makes the optimization problem much easier
  in general.

}

\seealso{
  \code{\link{normal1}},
  \code{\link[stats:Normal]{Normal}},
  \code{\link{mix2poisson}}.
}

\examples{
n = 1000
mu1 =  99   # Mean IQ of geography professors
mu2 = 150   # Mean IQ of mathematics professors
sd1 = sd2 = 16
phi = 0.3
y = ifelse(runif(n) < phi, rnorm(n, mu1, sd1), rnorm(n, mu2, sd2))

# Good idea to have trace=TRUE:
fit = vglm(y ~ 1, mix2normal1(esd=TRUE), maxit=200)
coef(fit, matrix=TRUE) # the estimates
c(phi, mu1, sd1, mu2, sd2) # the truth

\dontrun{# Plot the results
xx = seq(min(y), max(y), len=200)
plot(xx, (1-phi)*dnorm(xx, mu2, sd2), type="l", xlab="IQ",
     main="Red=estimate, blue=truth", col="blue")
phi.est = logit(coef(fit)[1], inverse=TRUE)
sd.est = exp(coef(fit)[3])
lines(xx, phi*dnorm(xx, mu1, sd1), col="blue")
lines(xx, phi.est * dnorm(xx, Coef(fit)[2], sd.est), col="red")
lines(xx, (1-phi.est) * dnorm(xx, Coef(fit)[4], sd.est), col="red")
abline(v=Coef(fit)[c(2,4)], lty=2, col="red")
abline(v=c(mu1, mu2), lty=2, col="blue")
}
}
\keyword{models}
\keyword{regression}


