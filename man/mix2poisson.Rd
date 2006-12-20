\name{mix2poisson}
\alias{mix2poisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mixture of Two Poisson Distributions }
\description{
  Estimates the three parameters of a mixture of two Poisson distributions
  by maximum likelihood estimation.

}
\usage{
mix2poisson(lphi = "logit", llambda = "loge",
            ephi=list(), el1=list(), el2=list(),
            iphi = 0.5, il1 = NULL, il2 = NULL,
            qmu = c(0.2, 0.8), zero = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi}{
    Link function for the parameter \eqn{\phi}{phi}.
    See below for more details.
    See \code{\link{Links}} for more choices.

  }
  \item{llambda}{
    Link function applied to each \eqn{\lambda}{lambda} parameter.
    See \code{\link{Links}} for more choices.

  }
  \item{ephi, el1, el2}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{iphi}{
  Initial value for \eqn{\phi}{phi}, whose value must lie
  between 0 and 1.

  }
  \item{il1, il2}{
    Optional initial value for \eqn{\lambda_1}{lambda1} and
    \eqn{\lambda_2}{lambda2}. These values must be positive.
    The default is to compute initial values internally using
    the argument \code{qmu}.

  }
  \item{qmu}{
    Vector with two values giving the probabilities relating to the sample
    quantiles for obtaining initial values for \eqn{\lambda_1}{lambda1}
    and \eqn{\lambda_2}{lambda2}.
    The two values are fed in as the \code{probs} argument into
    \code{\link[stats]{quantile}}.

  }
  \item{zero}{
  An integer specifying which linear/additive predictor is modelled as
  intercepts only.  If given, the value must be either 1 and/or 2 and/or
  3, and the default is the first one only, meaning \eqn{\phi}{phi}
  is a single parameter even when there are explanatory variables.
  Set \code{zero=NULL} to model all linear/additive predictors as
  functions of the explanatory variables.

  }
}
\details{
  The probability function can be loosely written as 
  \deqn{P(Y=y) = \phi \, Poisson(\lambda_1) + (1-\phi) \, Poisson(\lambda_2)}{%
        P(Y=y) = phi * Poisson(lambda1) + (1-phi) * Poisson(lambda2)}
  where \eqn{\phi}{phi} is the probability an observation belongs
  to the first group, and \eqn{y=0,1,2,\ldots}{y=0,1,2,...}.
  The parameter \eqn{\phi}{phi} satisfies \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is
  \eqn{\phi \lambda_1 + (1-\phi) \lambda_2}{phi*lambda1 + (1-phi)*lambda2}
  and this is returned as the fitted values.
  By default, the three linear/additive predictors are
  \eqn{(logit(\phi), \log(\lambda_1), \log(\lambda_2))^T}{(logit(phi),
  log(lambda1), log(lambda2))^T}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
% \references{ ~put references to the literature/web site here ~ }
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

  This \pkg{VGAM} family function should be used with caution.

}

\author{ T. W. Yee }
\note{
  Fitting this model successfully to data can be difficult due to
  numerical problems and ill-conditioned data.  It pays to fit the model
  several times with different initial values, and check that the best fit
  looks reasonable. Plotting the results is recommended.  This function
  works better as \eqn{\lambda_1}{lambda1} and \eqn{\lambda_2}{lambda2}
  become more different.

  Convergence is often slow, especially when the two component
  distributions are not well separated. The control argument \code{maxit}
  should be set to a higher value, e.g., 200, and use \code{trace=TRUE}
  to monitor convergence.

}

\seealso{
  \code{\link[stats:Poisson]{rpois}},
  \code{\link{mix2normal1}}.
}

\examples{
n = 3000
mu1 = exp(2.4) # also known as lambda1
mu2 = exp(3.1)
phi = 0.3
y = ifelse(runif(n) < phi, rpois(n, mu1), rpois(n, mu2))

fit = vglm(y ~ 1, mix2poisson, maxit=200) # good idea to have trace=TRUE
coef(fit, matrix=TRUE)
Coef(fit) # the estimates
c(phi, mu1, mu2) # the truth

\dontrun{# Plot the results
ty = table(y)
plot(names(ty), ty, type="h", main="Red=estimate, blue=truth")
abline(v=Coef(fit)[-1], lty=2, col="red")
abline(v=c(mu1, mu2), lty=2, col="blue")
}
}
\keyword{models}
\keyword{regression}


