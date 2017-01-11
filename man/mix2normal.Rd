\name{mix2normal}
\alias{mix2normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mixture of Two Univariate Normal Distributions }
\description{
  Estimates the five parameters of a mixture of two univariate
  normal distributions by maximum likelihood estimation.

}
\usage{
mix2normal(lphi = "logit", lmu = "identitylink", lsd = "loge",
           iphi = 0.5, imu1 = NULL, imu2 = NULL, isd1 = NULL, isd2 = NULL,
           qmu = c(0.2, 0.8), eq.sd = TRUE, nsimEIM = 100, zero = "phi")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi,lmu,lsd}{
    Link functions for the parameters \eqn{\phi}{phi},
    \eqn{\mu}{mu}, and
    \eqn{\sigma}{sd}.
    See \code{\link{Links}} for more choices.


  }

% \item{ephi, emu1, emu2, esd1, esd2}{
% List. Extra argument for each of the links.
% See \code{earg} in \code{\link{Links}} for general information.
% If \code{eq.sd = TRUE} then \code{esd1} must equal \code{esd2}.
% }

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
    Currently these are not great, therefore using these arguments
    where practical is a good idea.


  }
  \item{qmu}{
    Vector with two values giving the probabilities relating to the sample
    quantiles for obtaining initial values for \eqn{\mu_1}{mu1}
    and \eqn{\mu_2}{mu2}.
    The two values are fed in as the \code{probs} argument into
    \code{\link[stats]{quantile}}.


  }
  \item{eq.sd}{
    Logical indicating whether the two standard deviations should be
    constrained to be equal. If \code{TRUE} then the appropriate
    constraint matrices will be used.


  }
  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}}.


  }
  \item{zero}{
  May be an integer vector
  specifying which linear/additive predictors are modelled as
  intercept-only.  If given, the value or values can be from the
  set \eqn{\{1,2,\ldots,5\}}{1,2,...,5}.
  The default is the first one only, meaning \eqn{\phi}{phi}
  is a single parameter even when there are explanatory variables.
  Set \code{zero = NULL} to model all linear/additive predictors as
  functions of the explanatory variables.
  See \code{\link{CommonVGAMffArguments}} for more information.


}
}
\details{
  The probability density function can be loosely written as
  \deqn{f(y) = \phi \, N(\mu_1,\sigma_1) + (1-\phi) \, N(\mu_2, \sigma_2)}{%
        f(y) = phi * N(mu1, sd1) + (1-phi) * N(mu2, sd2)}
  where \eqn{\phi}{phi} is the probability an observation belongs
  to the first group.
  The parameters \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2} are the means, and
  \eqn{\sigma_1}{sd1} and \eqn{\sigma_2}{sd2} are the standard deviations.
  The parameter \eqn{\phi}{phi} satisfies \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is
  \eqn{\phi \mu_1 + (1-\phi) \mu_2}{phi*mu1 + (1-phi)*mu2}
  and this is returned as the fitted values.
  By default, the five linear/additive predictors are
  \eqn{(logit(\phi), \mu_1, \log(\sigma_1), \mu_2, \log(\sigma_2))^T}{
       (logit(phi),   mu1,   log(sd1), mu2, log(sd2))^T}.
  If \code{eq.sd = TRUE} then \eqn{\sigma_1 = \sigma_2}{sd1=sd2}
  is enforced.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


McLachlan, G. J. and Peel, D. (2000)
\emph{Finite Mixture Models}.
New York: Wiley.


Everitt, B. S. and Hand, D. J. (1981)
\emph{Finite Mixture Distributions}.
London: Chapman & Hall.


}
\section{Warning }{
  Numerical problems can occur and
  half-stepping is not uncommon.
  If failure to converge occurs, try inputting better initial values,
  e.g., by using \code{iphi},
  \code{qmu},
  \code{imu1},
  \code{imu2},
  \code{isd1},
  \code{isd2},
  etc.


% This function uses a quasi-Newton update for the working weight matrices
% (BFGS variant). It builds up approximations to the weight matrices,
% and currently the code is not fully tested.
% In particular, results based on the weight matrices (e.g., from
% \code{vcov} and \code{summary}) may be quite incorrect, especially when
% the arguments \code{weights} is used to input prior weights.


  This \pkg{VGAM} family function is experimental and
  should be used with care.


}

\author{ T. W. Yee }
\note{
  Fitting this model successfully to data can be difficult due to
  numerical problems and ill-conditioned data.  It pays to fit the
  model several times with different initial values and check that
  the best fit looks reasonable. Plotting the results is recommended.
  This function works better as \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2}
  become more different.


  Convergence can be slow, especially when the two component
  distributions are not well separated.
  The default control argument \code{trace = TRUE} is to encourage
  monitoring convergence.
  Having \code{eq.sd = TRUE} often makes the overall optimization problem
  easier.


}

\seealso{
  \code{\link{uninormal}},
  \code{\link[stats:Normal]{Normal}},
  \code{\link{mix2poisson}}.


}

\examples{
\dontrun{ mu1 <-  99; mu2 <- 150; nn <- 1000
sd1 <- sd2 <- exp(3)
(phi <- logit(-1, inverse = TRUE))
mdata <- data.frame(y = ifelse(runif(nn) < phi, rnorm(nn, mu1, sd1),
                                                rnorm(nn, mu2, sd2)))
fit <- vglm(y ~ 1, mix2normal(eq.sd = TRUE), data = mdata)

# Compare the results
cfit <- coef(fit)
round(rbind('Estimated' = c(logit(cfit[1], inverse = TRUE),
            cfit[2], exp(cfit[3]), cfit[4]),
            'Truth' = c(phi, mu1, sd1, mu2)), digits = 2)

# Plot the results
xx <- with(mdata, seq(min(y), max(y), len = 200))
plot(xx, (1-phi) * dnorm(xx, mu2, sd2), type = "l", xlab = "y",
     main = "Orange = estimate, blue = truth",
     col = "blue", ylab = "Density")
phi.est <- logit(coef(fit)[1], inverse = TRUE)
sd.est <- exp(coef(fit)[3])
lines(xx, phi*dnorm(xx, mu1, sd1), col = "blue")
lines(xx, phi.est * dnorm(xx, Coef(fit)[2], sd.est), col = "orange")
lines(xx, (1-phi.est) * dnorm(xx, Coef(fit)[4], sd.est), col = "orange")
abline(v = Coef(fit)[c(2,4)], lty = 2, col = "orange")
abline(v = c(mu1, mu2), lty = 2, col = "blue")
}
}
\keyword{models}
\keyword{regression}


