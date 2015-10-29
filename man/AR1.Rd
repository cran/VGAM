\name{AR1}
\alias{AR1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Autoregressive Process with Order-1 Family Function }
\description{
  Maximum likelihood estimation of the three-parameter AR-1 model


}
\usage{
AR1(ldrift = "identitylink", lsd  = "loge", lvar = "loge",
    lrho = "rhobit", idrift  = NULL,
    isd  = NULL, ivar = NULL, irho = NULL,
    ishrinkage = 0.9, type.likelihood = c("exact", "conditional"),
    var.arg = FALSE, nodrift = FALSE, almost1 = 0.99, zero = c(-2, -3))
}

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

  \item{ishrinkage, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{var.arg}{
  Same meaning as \code{\link{uninormal}}.


  }
  \item{nodrift}{
  Logical, for determining whether to estimate the drift parameter.
  The default is to estimate it.
  If \code{TRUE}, the drift parameter is set to 0 and not estimated.


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


\item{almost1}{
  A value close to 1 but slightly smaller. One of the off-diagonal
  EIM elements is multiplied by this, to ensure that the EIM is
  positive-definite.



}

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
  1. A test for stationarity might be to test
  whether \eqn{\mu^*}{mu^*} is intercept-only.
  2. The mean of all the \eqn{Y_i}{Y(i)}
  is \eqn{\mu^* /(1-\rho)}{mu^* / (1-rho)} and
  these are returned as the fitted values.
  3. The correlation of all the \eqn{Y_i}{Y(i)} with \eqn{Y_{i-1}}{Y(i-1)}
  is \eqn{\rho}{rho}.
  4. The default link function ensures that \eqn{-1 < \rho < 1}{-1 < rho < 1}.

  
}
\section{Warning}{
  Monitoring convergence is urged: set \code{trace = TRUE}.

  

  Yet to do: add an argument that allows the scaled mean parameter
  to be deleted, i.e, a 2-parameter model is fitted.
  Yet to do: \code{ARff(p.lag = 1)} should hopefully be written soon.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
%\references{


%}
\author{ Thomas W. Yee and Victor Miranda }
\note{
  For \code{type.likelihood = "conditional"},
  the prior weight for the first observation is set to
  some small positive number, which has the effect of deleting
  that observation.
  However, \eqn{n} is still the original \eqn{n} so that
  statistics such as the residual degrees of freedom are
  unchanged (uncorrected possibly).


  Multiple responses are handled.
  The mean is returned as the fitted values.


  Practical experience has shown that half-stepping is a very
  good idea. The default options use step sizes that are
  about one third the usual step size. Consequently,
  \code{maxit} is increased to about 100, by default.


}
  


\seealso{
  \code{\link{vglm.control}},
  \code{\link{dAR1}},
  \code{\link{uninormal}},
  \code{\link[stats]{arima.sim}},


}
\examples{
# Example 1: using  arima.sim() to generate a stationary time series
nn <- 100; set.seed(1)
tsdata <- data.frame(x2 =  runif(nn))
tsdata  <- transform(tsdata,
              index = 1:nn,
              TS1 = arima.sim(nn, model = list(ar = -0.80),
                              sd = exp(1.0)),
              TS2 = arima.sim(nn, model = list(ar =  0.50),
                              sd = exp(1.0 + 2 * x2)))
fit1a <- vglm(cbind(TS1, TS2) ~ x2, AR1(zero = c(1:4, 6)),
             data = tsdata, trace = TRUE)
rhobit(-0.8)
rhobit( 0.5)
coef(fit1a, matrix = TRUE)
summary(fit1a)  # SEs are useful to know

# Example 2: another stationary time series
nn <- 1000
my.rho <- rhobit(-1.0, inverse = TRUE)
my.mu  <- 2.5
my.sd  <- exp(1)
tsdata  <- data.frame(index = 1:nn, TS3 = runif(nn))
for (ii in 2:nn)
  tsdata$TS3[ii] <- my.mu + my.rho * tsdata$TS3[ii-1] + rnorm(1, sd = my.sd)
tsdata <- tsdata[-(1:ceiling(nn/5)), ]  # Remove the burn-in data:
fit2a <- vglm(TS3 ~ 1, AR1(type.likelihood = "conditional"),
             data = tsdata, trace = TRUE)
coef(fit2a, matrix = TRUE)
summary(fit2a)  # SEs are useful to know
Coef(fit2a)["rho"]  # Estimate of rho for intercept-only models
my.rho
coef(fit2a)[1]  # drift
my.mu           # Should be the same
head(weights(fit2a, type= "prior"))    # First one is effectively deleted
head(weights(fit2a, type= "working"))  # Ditto
}
\keyword{models}
\keyword{regression}



