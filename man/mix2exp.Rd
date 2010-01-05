\name{mix2exp}
\alias{mix2exp}
%- Also NEED an '\alias' for EACH other topic documented here.
%- Adapted from mix2poisson.Rd
\title{ Mixture of Two Exponential Distributions }
\description{
  Estimates the three parameters of a mixture of two exponential
  distributions by maximum likelihood estimation.

}
\usage{
mix2exp(lphi="logit", llambda="loge", ephi=list(), el1=list(), el2=list(),
        iphi=0.5, il1=NULL, il2=NULL, qmu=c(0.8, 0.2), nsimEIM=100, zero=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi}{
    Link function for the parameter \eqn{\phi}{phi}.
    See \code{\link{Links}} for more choices.

  }
  \item{llambda}{
    Link function applied to each \eqn{\lambda}{lambda} (rate) parameter.
    See \code{\link{Links}} for more choices.
    Note that the mean of an ordinary exponential distribution is
    \eqn{1 / \lambda}.

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
  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}}.

  }
  \item{zero}{
  An integer specifying which linear/additive predictor is modelled as
  intercepts only.  If given, the value must be either 1 and/or 2 and/or
  3, and the default is the first one only, meaning \eqn{\phi}{phi}
  is a single parameter even when there are explanatory variables.
  Set \code{zero=NULL} to model all linear/additive predictors as
  functions of the explanatory variables.
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
}
\details{
  The probability function can be loosely written as 
  \deqn{P(Y=y) = \phi\,Exponential(\lambda_1) + (1-\phi)\,Exponential(\lambda_2)}{%
        P(Y=y) = phi * Exponential(lambda1) + (1-phi) * Exponential(lambda2)}
  where \eqn{\phi}{phi} is the probability an observation belongs
  to the first group, and \eqn{y>0}.
  The parameter \eqn{\phi}{phi} satisfies \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is
  \eqn{\phi / \lambda_1 + (1-\phi) / \lambda_2}{phi/lambda1 + (1-phi)/lambda2}
  and this is returned as the fitted values.
  By default, the three linear/additive predictors are
  \eqn{(logit(\phi), \log(\lambda_1), \log(\lambda_2))^T}{(logit(phi),
  log(lambda1), log(lambda2))^T}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
% \references{ ~put references to the literature/web site here ~ }
\section{Warning }{
  This \pkg{VGAM} family function requires care for a successful
  application.
  In particular, good initial values are required because of the presence
  of local solutions. Therefore running this function with several
  different combinations of arguments such as \code{iphi}, \code{il1},
  \code{il2}, \code{qmu} is highly recommended. Graphical methods such
  as \code{\link[graphics]{hist}} can be used as an aid.

}

\author{ T. W. Yee }
\note{
  Fitting this model successfully to data can be difficult due to
  local solutions, uniqueness problems and ill-conditioned data. It
  pays to fit the model several times with different initial values
  and check that the best fit looks reasonable. Plotting the results is
  recommended. This function works better as \eqn{\lambda_1}{lambda1}
  and \eqn{\lambda_2}{lambda2} become more different.
  The default control argument \code{trace=TRUE} is to encourage
  monitoring convergence.

}

\seealso{
  \code{\link[stats:Exponential]{rexp}},
  \code{\link{exponential}},
  \code{\link{mix2poisson}}.
}

\examples{
lambda1 = exp(1); lambda2 = exp(3)
(phi = logit(-1, inverse=TRUE))
mdata = data.frame(y1 = rexp(nn <- 1000, lambda1))
mdata = transform(mdata, y2 = rexp(nn, lambda2))
mdata = transform(mdata, y  = ifelse(runif(nn) < phi, y1, y2))
fit = vglm(y ~ 1, mix2exp, mdata, trace=TRUE)
coef(fit, matrix=TRUE)

# Compare the results with the truth
round(rbind('Estimated'=Coef(fit),
            'Truth'=c(phi, lambda1, lambda2)), dig=2)

\dontrun{# Plot the results
with(mdata, hist(y, prob=TRUE, main="Red=estimate, blue=truth"))
abline(v=1/Coef(fit)[c(2,3)], lty=2, col="red", lwd=2)
abline(v=1/c(lambda1, lambda2), lty=2, col="blue", lwd=2) }
}
\keyword{models}
\keyword{regression}


