\name{AR1}
\alias{AR1}
% \alias{AR1.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Autoregressive Process with Order-1 Family Function }
\description{
  Maximum likelihood estimation of the three-parameter AR-1 model


}
\usage{
AR1(ldrift = "identitylink", lsd  = "loge", lvar = "loge", lrho = "rhobit",
    idrift  = NULL, isd  = NULL, ivar = NULL, irho = NULL, imethod = 1,
    ishrinkage = 0.95, type.likelihood = c("exact", "conditional"),
    type.EIM  = c("exact", "approximate"), var.arg = FALSE, nodrift = FALSE,
    print.EIM = FALSE, zero = c(if (var.arg) "var" else "sd", "rho"))

}
%   zero = c(-2, -3)
% AR1.control(epsilon = 1e-6, maxit = 30, stepsize = 1, ...)

%     deviance.arg = FALSE,

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ldrift, lsd, lvar, lrho}{
    Link functions applied to the scaled mean, standard deviation
    or variance, and correlation parameters.
    The parameter \code{drift} is known as the \emph{drift}, and
    it is a scaled mean.
    See \code{\link{Links}} for more choices.


  }
  \item{idrift, isd, ivar, irho}{
  Optional initial values for the parameters.
  If failure to converge occurs then try different values
  and monitor convergence by using \code{trace = TRUE}.
  For a \eqn{S}-column response, these arguments can be of length
  \eqn{S}, and they are recycled by the columns first.
  A value \code{NULL} means an initial value for each response is
  computed internally.



  }

  \item{ishrinkage, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.
  The default for \code{zero} assumes there is a drift parameter to
  be estimated (the default for that argument), so if a drift parameter
  is suppressed and there are covariates, then \code{zero} will need
  to be assigned the value 1 or 2 or \code{NULL}.


  }
  \item{var.arg}{
  Same meaning as \code{\link{uninormal}}.


  }
  \item{nodrift}{
  Logical, for determining whether to estimate the drift parameter.
  The default is to estimate it.
  If \code{TRUE}, the drift parameter is set to 0 and not estimated.


  }
  \item{type.EIM}{
   What type of expected information matrix (EIM) is used in
   Fisher scoring. By default, this family function calls
   \code{\link[VGAM:AR1EIM]{AR1EIM}}, which recursively
   computes the exact EIM for the AR process with Gaussian
   white noise. See Porat and Friedlander (1986) for further
   details on the exact EIM.

   If \code{type.EIM = "approximate"} then
   approximate expression for the EIM of Autoregressive processes
   is used; this approach holds when the number of observations
   is large enough. Succinct details about the approximate EIM
   are delineated at Porat and Friedlander (1987).


  }
  \item{print.EIM}{
  Logical. If \code{TRUE}, then the first few EIMs are printed.
  Here, the result shown is the sum of each EIM.


  }
  \item{type.likelihood}{
    What type of likelihood function is maximized.
    The first choice (default) is the sum of the marginal likelihood
    and the conditional likelihood.
    Choosing the conditional likelihood means that the first observation is
    effectively ignored (this is handled internally by setting
    the value of the first prior weight to be some small
    positive number, e.g., \code{1.0e-6}).
    See the note below.


  }
%  \item{epsilon, maxit, stepsize,...}{
%  Same as \code{\link[VGAM:vglm.control]{vglm.control}}.
%
%
%  }


}
\details{
  The AR-1 model implemented here has
  \deqn{Y_1 \sim N(\mu, \sigma^2 / (1-\rho^2)), }{%
    Y(1) ~ N(mu, sigma^2 / (1-rho^2), }
  and
  \deqn{Y_i  = \mu^* + \rho Y_{i-1} + e_i, }{%
    Y(i) =  mu^* +  rho * Y(i-1) + e(i) }
  where the \eqn{e_i}{e(i)} are i.i.d. Normal(0, sd = \eqn{\sigma}{sigma})
  random variates.


  Here are a few notes:
  (1). A test for weak stationarity might be to verify whether
  \eqn{1/\rho}{1/rho} lies outside the unit circle.
  (2). The mean of all the \eqn{Y_i}{Y(i)}
  is \eqn{\mu^* /(1-\rho)}{mu^* / (1-rho)} and
  these are returned as the fitted values.
  (3). The correlation of all the \eqn{Y_i}{Y(i)} with \eqn{Y_{i-1}}{Y(i-1)}
  is \eqn{\rho}{rho}.
  (4). The default link function ensures that
                         \eqn{-1 < \rho < 1}{-1 < rho < 1}.



% (1). ...  whether \eqn{\mu^*}{mu^*} is intercept-only.


}
\section{Warning}{
  Monitoring convergence is urged, i.e., set \code{trace = TRUE}.

  Moreover, if the exact EIMs are used, set \code{print.EIM = TRUE}
  to compare the computed exact to the approximate EIM.

  Under the VGLM/VGAM approach, parameters can be modelled in terms
  of covariates. Particularly, if the standard deviation of
  the white noise is modelled in this way, then
  \code{type.EIM = "exact"} may certainly lead to unstable
  results. The reason is that white noise is a stationary
  process, and consequently, its variance must remain as a constant.
  Consequently, the use of variates to model
  this parameter contradicts the assumption of
  stationary random components to compute the exact EIMs proposed
  by Porat and Friedlander (1987).

  To prevent convergence issues in such cases, this family function
  internally verifies whether the variance of the white noise remains
  as a constant at each Fisher scoring iteration.
  If this assumption is violated and \code{type.EIM = "exact"} is set,
  then \code{AR1} automatically shifts to
  \code{type.EIM = "approximate"}.
  Also, a warning is accordingly displayed.


  %Thus, if modelling the standard deviation of the white noise
  %is required, the use of \code{type.EIM = "approximate"} is
  %highly recommended.

  %Yet to do: add an argument that allows the scaled mean parameter
  %to be deleted, i.e, a 2-parameter model is fitted.
  %Yet to do: \code{ARff(p.lag = 1)} should hopefully be written soon.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Porat, B. and Friedlander, B. (1987)
  The Exact Cramer-Rao Bond for Gaussian Autoregressive Processes.
  \emph{IEEE Transactions on Aerospace and Electronic Systems},
  \bold{AES-23(4)}, 537--542.


}
\author{ Victor Miranda (exact method) and
         Thomas W. Yee (approximate method).}


\note{
  %For \code{type.likelihood = "conditional"},
  %the prior weight for the first observation is set to
  %some small positive number, which has the effect of deleting
  %that observation.
  %However, \eqn{n} is still the original \eqn{n} so that
  %statistics such as the residual degrees of freedom are
  %unchanged (uncorrected possibly).


  Multiple responses are handled.
  The mean is returned as the fitted values.


% Argument \code{zero} can be either a numeric or a character
% vector. It must specify the position(s) or name(s) of the
% parameters to be modeled as intercept-only. If names are used,
% notice that parameter names in this family function are


% \deqn{c("drift", "var" or "sd", "rho").}


  %Practical experience has shown that half-stepping is a very
  %good idea. The default options use step sizes that are
  %about one third the usual step size. Consequently,
  %\code{maxit} is increased to about 100, by default.



}

\seealso{
  \code{\link{AR1EIM}},
  \code{\link{vglm.control}},
  \code{\link{dAR1}},
  \code{\link[stats]{arima.sim}}.


}
\examples{
### Example 1: using  arima.sim() to generate a 0-mean stationary time series.
nn <- 500
tsdata <- data.frame(x2 =  runif(nn))
ar.coef.1 <- rhobit(-1.55, inverse = TRUE)  # Approx -0.65
ar.coef.2 <- rhobit( 1.0, inverse = TRUE)   # Approx  0.50
set.seed(1)
tsdata  <- transform(tsdata,
              index = 1:nn,
              TS1 = arima.sim(nn, model = list(ar = ar.coef.1),
                              sd = exp(1.5)),
              TS2 = arima.sim(nn, model = list(ar = ar.coef.2),
                              sd = exp(1.0 + 1.5 * x2)))

### An autoregressive intercept--only model.   ###
### Using the exact EIM, and "nodrift = TRUE"  ###
fit1a <- vglm(TS1 ~ 1, data = tsdata, trace = TRUE,
              AR1(var.arg = FALSE, nodrift = TRUE,
                  type.EIM = "exact",
                  print.EIM = FALSE),
              crit = "coefficients")
Coef(fit1a)
summary(fit1a)

\dontrun{
### Two responses. Here, the white noise standard deviation of TS2   ###
### is modelled in terms of 'x2'. Also, 'type.EIM = exact'.  ###
fit1b <- vglm(cbind(TS1, TS2) ~ x2,
              AR1(zero = NULL, nodrift = TRUE,
                  var.arg = FALSE,
                  type.EIM = "exact"),
              constraints = list("(Intercept)" = diag(4),
                                 "x2" = rbind(0, 0, 1, 0)),
              data = tsdata, trace = TRUE, crit = "coefficients")
coef(fit1b, matrix = TRUE)
summary(fit1b)

### Example 2: another stationary time series
nn     <- 500
my.rho <- rhobit(1.0, inverse = TRUE)
my.mu  <- 1.0
my.sd  <- exp(1)
tsdata  <- data.frame(index = 1:nn, TS3 = runif(nn))

set.seed(2)
for (ii in 2:nn)
  tsdata$TS3[ii] <- my.mu/(1 - my.rho) +
                    my.rho * tsdata$TS3[ii-1] + rnorm(1, sd = my.sd)
tsdata <- tsdata[-(1:ceiling(nn/5)), ]  # Remove the burn-in data:

### Fitting an AR(1). The exact EIMs are used.
fit2a <- vglm(TS3 ~ 1, AR1(type.likelihood = "exact",  # "conditional",
                                type.EIM = "exact"),
              data = tsdata, trace = TRUE, crit = "coefficients")

Coef(fit2a)
summary(fit2a)      # SEs are useful to know

Coef(fit2a)["rho"]    # Estimate of rho, for intercept-only models
my.rho                # The 'truth' (rho)
Coef(fit2a)["drift"]  # Estimate of drift, for intercept-only models
my.mu /(1 - my.rho)   # The 'truth' (drift)
}
}
\keyword{models}
\keyword{regression}



%fit1a <- vglm(cbind(TS1, TS2) ~ x2, AR1(zero = c(1:4, 6)),
%             data = tsdata, trace = TRUE)


%head(weights(fit2a, type = "prior"))    # First one is effectively deleted
%head(weights(fit2a, type = "working"))  # Ditto



