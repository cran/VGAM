\name{dAR1}
\alias{dAR1}
\alias{dAR1}
%\alias{pAR1}
%\alias{qAR1}
%\alias{rAR1}
\title{The AR-1 Autoregressive Process}
\description{
  Density for the AR-1 model.

}
\usage{
dAR1(x, drift = 0, var.error = 1, ARcoef1 = 0.0,
     type.likelihood = c("exact", "conditional"), log = FALSE)
}
\arguments{
  \item{x,}{vector of quantiles.}

  \item{drift}{
    the scaled mean (also known as the \emph{drift} parameter),
    \eqn{\mu^*}{mu^*}.
    Note that the mean is \eqn{\mu^* /(1-\rho)}{mu^* / (1-rho)}.
    The default corresponds to observations that have mean 0.



  }
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm of the density is returned.


  }
  \item{type.likelihood, var.error, ARcoef1}{
    See \code{\link{AR1}}.
    The argument \code{ARcoef1} is \eqn{\rho}{rho}.
    The argument \code{var.error} is the variance of  the
    i.i.d. random noise, i.e., \eqn{\sigma^2}{sigma^2}.
    If \code{type.likelihood = "conditional"} then the
    first element or row of the result is currently
    assigned \code{NA}---this
    is because the density of the first observation is effectively
    ignored.


  }
}
\value{
  \code{dAR1} gives the density.



% \code{pAR1} gives the distribution function, and
% \code{qAR1} gives the quantile function, and
% \code{rAR1} generates random deviates.


}

\author{ T. W. Yee and Victor Miranda }

\details{
  Most of the background to this function is given
  in \code{\link{AR1}}.
  All the arguments are converted into matrices, and then
  all their dimensions are obtained. They are then coerced
  into the same size: the number of rows is the maximum
  of all the single rows, and ditto for the number of columns.



}
%\note{
%}
\seealso{
  \code{\link{AR1}}.


}
\examples{
nn <- 100; set.seed(1)
tdata <- data.frame(index = 1:nn,
                    TS1 = arima.sim(nn, model = list(ar = -0.50),
                                    sd = exp(1)))
fit1 <- vglm(TS1 ~ 1, AR1, data = tdata, trace = TRUE)
rhobit(-0.5)
coef(fit1, matrix = TRUE)
(Cfit1 <- Coef(fit1))
summary(fit1)  # SEs are useful to know
logLik(fit1)
sum(dAR1(depvar(fit1), drift = Cfit1[1], var.error = (Cfit1[2])^2,
         ARcoef1 = Cfit1[3], log = TRUE))

fit2 <- vglm(TS1 ~ 1, AR1(type.likelihood = "cond"), data = tdata, trace = TRUE)
(Cfit2 <- Coef(fit2))  # Okay for intercept-only models
logLik(fit2)
head(keep <- dAR1(depvar(fit2), drift = Cfit2[1], var.error = (Cfit2[2])^2,
                  ARcoef1 = Cfit2[3], type.likelihood = "cond", log = TRUE))
sum(keep[-1])
}
\keyword{distribution}


