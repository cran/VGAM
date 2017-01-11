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
            iphi = 0.5, il1 = NULL, il2 = NULL,
            qmu = c(0.2, 0.8), nsimEIM = 100, zero = "phi")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi, llambda}{
    Link functions for the parameter \eqn{\phi}{phi} and
    \eqn{\lambda}{lambda}.
    See \code{\link{Links}} for more choices.

  }

% \item{ephi, el1, el2}{
%           ephi = list(), el1 = list(), el2 = list(),
% List. Extra argument for each of the links.
% See \code{earg} in \code{\link{Links}} for general information.
% }

  \item{iphi}{
  Initial value for \eqn{\phi}{phi}, whose value must lie
  between 0 and 1.

  }
  \item{il1, il2}{
    Optional initial value for \eqn{\lambda_1}{lambda1} and
    \eqn{\lambda_2}{lambda2}. These values must be positive.
    The default is to compute initial values internally using
    the argument \code{qmu}.
%   If these arguments are supplied then practical experience
%   suggests they should be quite well-separated.

  }
  \item{qmu}{
    Vector with two values giving the probabilities relating to the sample
    quantiles for obtaining initial values for \eqn{\lambda_1}{lambda1}
    and \eqn{\lambda_2}{lambda2}.
    The two values are fed in as the \code{probs} argument into
    \code{\link[stats]{quantile}}.


  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}}.


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


  With grouped data (i.e., using the \code{weights} argument)
  one has to use a large value of \code{nsimEIM};
  see the example below.




  This \pkg{VGAM} family function is experimental and
  should be used with care.


}

\author{ T. W. Yee }
\note{
  The response must be integer-valued since \code{\link[stats]{dpois}}
  is invoked.


  Fitting this model successfully to data can be difficult due to local
  solutions and ill-conditioned data. It pays to fit the model several
  times with different initial values, and check that the best fit
  looks reasonable. Plotting the results is recommended. This function
  works better as \eqn{\lambda_1}{lambda1} and \eqn{\lambda_2}{lambda2}
  become more different.
  The default control argument \code{trace = TRUE} is to encourage
  monitoring convergence.


}

\seealso{
  \code{\link[stats:Poisson]{rpois}},
  \code{\link{poissonff}},
  \code{\link{mix2normal}}.


}

\examples{
\dontrun{ # Example 1: simulated data
nn <- 1000
mu1 <- exp(2.5)  # Also known as lambda1
mu2 <- exp(3)
(phi <- logit(-0.5, inverse = TRUE))
mdata <- data.frame(y = rpois(nn, ifelse(runif(nn) < phi, mu1, mu2)))
mfit <- vglm(y ~ 1, mix2poisson, data = mdata)
coef(mfit, matrix = TRUE)

# Compare the results with the truth
round(rbind('Estimated' = Coef(mfit), 'Truth' = c(phi, mu1, mu2)), digits = 2)

ty <- with(mdata, table(y))
plot(names(ty), ty, type = "h", main = "Orange=estimate, blue=truth",
     ylab = "Frequency", xlab = "y")
abline(v = Coef(mfit)[-1], lty = 2, col = "orange", lwd = 2)
abline(v = c(mu1, mu2), lty = 2, col = "blue", lwd = 2)

# Example 2: London Times data (Lange, 1997, p.31)
ltdata1 <- data.frame(deaths = 0:9,
                      freq = c(162, 267, 271, 185, 111, 61, 27, 8, 3, 1))
ltdata2 <- data.frame(y = with(ltdata1, rep(deaths, freq)))

# Usually this does not work well unless nsimEIM is large
Mfit <- vglm(deaths ~ 1, weight = freq, data = ltdata1,
            mix2poisson(iphi = 0.3, il1 = 1, il2 = 2.5, nsimEIM = 5000))

# This works better in general
Mfit <- vglm(y ~ 1, mix2poisson(iphi = 0.3, il1 = 1, il2 = 2.5), data = ltdata2)
coef(Mfit, matrix = TRUE)
Coef(Mfit)
}
}
\keyword{models}
\keyword{regression}


