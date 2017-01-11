\name{tobit}
\alias{tobit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tobit Model }
\description{
  Fits a Tobit model.

}
\usage{
tobit(Lower = 0, Upper = Inf, lmu = "identitylink", lsd = "loge",
      imu = NULL, isd = NULL,
      type.fitted = c("uncensored", "censored", "mean.obs"),
      byrow.arg = FALSE, imethod = 1, zero = "sd")
}
% 20151024 yettodo: maybe add a new option to 'type.fitted':
%     type.fitted = c("uncensored", "censored", "mean.obs", "truncated"),
% where "truncated" is only concerned with values of y > Lower;
% values of y <= Lower are ignored.
%
%
%
%
%
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lower}{
  Numeric. It is the value \eqn{L} described below.
  Any value of the linear model
  \eqn{x_i^T \beta}{x_i^T beta} that
  is less than this lowerbound is assigned this value.
  Hence this should be the smallest possible value in the response
  variable.
  May be a vector (see below for more information).


  }
  \item{Upper}{
  Numeric. It is the value \eqn{U} described below.
  Any value of the linear model
  \eqn{x_i^T \beta}{x_i^T beta} that
  is greater than this upperbound is assigned this value.
  Hence this should be the largest possible value in the response
  variable.
  May be a vector (see below for more information).


  }
  \item{lmu, lsd}{
  Parameter link functions for the mean and standard deviation parameters.
  See \code{\link{Links}} for more choices.
  The standard deviation is a positive quantity, therefore a log link
  is its default.


  }

  \item{imu, isd, byrow.arg}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{type.fitted}{
  Type of fitted value returned.
  The first choice is default and is the ordinary uncensored or
  unbounded linear model.
  If \code{"censored"} then the fitted values in the interval \eqn{[L, U]}.
  If \code{"mean.obs"} then the mean of the observations is returned;
  this is a doubly truncated normal distribution
  augmented by point masses at the truncation points
  (see \code{\link{dtobit}}).
  See \code{\link{CommonVGAMffArguments}} for more information.



  }
  \item{imethod}{
  Initialization method. Either 1 or 2 or 3, this specifies
  some methods for obtaining initial values for the parameters.
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{zero}{
  A vector, e.g., containing the value 1 or 2. If so,
  the mean or standard deviation respectively are modelled
  as an intercept-only.
  Setting \code{zero = NULL} means both linear/additive predictors
  are modelled as functions of the explanatory variables.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The Tobit model can be written
  \deqn{y_i^* = x_i^T \beta + \varepsilon_i}{%
        y_i^* = x_i^T  beta +  e_i}
  where the \eqn{e_i \sim N(0,\sigma^2)}{e_i ~ N(0,sigma^2)}
  independently and \eqn{i=1,\ldots,n}{i=1,...,n}.
  However, we measure \eqn{y_i = y_i^*} only if \eqn{y_i^* > L} and
  \eqn{y_i^* < U} for some
  cutpoints \eqn{L} and \eqn{U}. Otherwise we let \eqn{y_i=L} or
  \eqn{y_i=U}, whatever is closer.
  The Tobit model is thus a multiple linear regression but with censored
  responses if it is below or above certain cutpoints.



  The defaults for \code{Lower} and \code{Upper} and
  \code{lmu} correspond to the \emph{standard} Tobit model.
  Fisher scoring is used for the standard and nonstandard
  models.
  By default, the mean \eqn{x_i^T \beta}{x_i^T beta} is
  the first linear/additive predictor, and the log of
  the standard deviation is the second linear/additive
  predictor. The Fisher information matrix for uncensored
  data is diagonal. The fitted values are the estimates
  of \eqn{x_i^T \beta}{x_i^T beta}.



}
\section{Warning }{
  If values of the response and \code{Lower} and/or \code{Upper}
  are not integers then there is the danger that the value is
  wrongly interpreted as uncensored.
  For example, if the first 10 values of the response were
  \code{runif(10)} and \code{Lower} was assigned these value then
  testing \code{y[1:10] == Lower[1:10]} is numerically fraught.
  Currently, if any \code{y < Lower} or \code{y > Upper} then
  a warning is issued.


}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Tobin, J. (1958)
  Estimation of relationships for limited dependent variables.
  \emph{Econometrica} \bold{26}, 24--36.


}

\author{ Thomas W. Yee }
\note{
  The response can be a matrix.
  If so, then \code{Lower} and \code{Upper}
  are recycled into a matrix with the number of columns equal
  to the number of responses,
  and the recycling is done row-wise \emph{if} \code{byrow.arg = TRUE}.
  The default order is as \code{\link[base]{matrix}}, which
  is \code{byrow.arg = FALSE}.
  For example, these are returned in \code{fit4@misc$Lower} and
  \code{fit4@misc$Upper} below.



  If there is no censoring then
  \code{\link{uninormal}} is recommended instead. Any value of the
  response less than \code{Lower} or greater than \code{Upper} will
  be assigned the value \code{Lower} and \code{Upper} respectively,
  and a warning will be issued.
  The fitted object has components \code{censoredL} and \code{censoredU}
  in the \code{extra} slot which specifies whether observations
  are censored in that direction.
  The function \code{\link{cens.normal}} is an alternative
  to \code{tobit()}.



% 20150417; McClelland Kemp bug:


  When obtaining initial values, if the algorithm would
  otherwise want to fit an underdetermined system of
  equations, then it uses the entire data set instead.
  This might result in rather poor quality initial values,
  and consequently, monitoring convergence is advised.



}
\seealso{
  \code{\link{rtobit}},
  \code{\link{cens.normal}},
  \code{\link{uninormal}},
  \code{\link{double.cens.normal}},
  \code{\link{posnormal}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link[stats:Normal]{rnorm}}.



}
\examples{
# Here, fit1 is a standard Tobit model and fit2 is a nonstandard Tobit model
tdata <- data.frame(x2 = seq(-1, 1, length = (nn <- 100)))
set.seed(1)
Lower <- 1; Upper <- 4  # For the nonstandard Tobit model
tdata <- transform(tdata,
                   Lower.vec = rnorm(nn, Lower, 0.5),
                   Upper.vec = rnorm(nn, Upper, 0.5))
meanfun1 <- function(x) 0 + 2*x
meanfun2 <- function(x) 2 + 2*x
meanfun3 <- function(x) 2 + 2*x
meanfun4 <- function(x) 3 + 2*x
tdata <- transform(tdata,
  y1 = rtobit(nn, mean = meanfun1(x2)),  # Standard Tobit model
  y2 = rtobit(nn, mean = meanfun2(x2), Lower = Lower, Upper = Upper),
  y3 = rtobit(nn, mean = meanfun3(x2), Lower = Lower.vec, Upper = Upper.vec),
  y4 = rtobit(nn, mean = meanfun3(x2), Lower = Lower.vec, Upper = Upper.vec))
with(tdata, table(y1 == 0))  # How many censored values?
with(tdata, table(y2 == Lower | y2 == Upper))  # How many censored values?
with(tdata, table(attr(y2, "cenL")))
with(tdata, table(attr(y2, "cenU")))

fit1 <- vglm(y1 ~ x2, tobit, data = tdata, trace = TRUE)
coef(fit1, matrix = TRUE)
summary(fit1)

fit2 <- vglm(y2 ~ x2, tobit(Lower = Lower, Upper = Upper, type.f = "cens"),
             data = tdata, trace = TRUE)
table(fit2@extra$censoredL)
table(fit2@extra$censoredU)
coef(fit2, matrix = TRUE)

fit3 <- vglm(y3 ~ x2, tobit(Lower = with(tdata, Lower.vec),
                            Upper = with(tdata, Upper.vec), type.f = "cens"),
             data = tdata, trace = TRUE)
table(fit3@extra$censoredL)
table(fit3@extra$censoredU)
coef(fit3, matrix = TRUE)

# fit4 is fit3 but with type.fitted = "uncen".
fit4 <- vglm(cbind(y3, y4) ~ x2,
             tobit(Lower = rep(with(tdata, Lower.vec), each = 2),
                   Upper = rep(with(tdata, Upper.vec), each = 2),
                   byrow.arg = TRUE),
             data = tdata, crit = "coeff", trace = TRUE)
head(fit4@extra$censoredL)  # A matrix
head(fit4@extra$censoredU)  # A matrix
head(fit4@misc$Lower)       # A matrix
head(fit4@misc$Upper)       # A matrix
coef(fit4, matrix = TRUE)

\dontrun{ # Plot fit1--fit4
par(mfrow = c(2, 2))

plot(y1 ~ x2, tdata, las = 1, main = "Standard Tobit model",
     col = as.numeric(attr(y1, "cenL")) + 3,
     pch = as.numeric(attr(y1, "cenL")) + 1)
legend(x = "topleft", leg = c("censored", "uncensored"),
       pch = c(2, 1), col = c("blue", "green"))
legend(-1.0, 2.5, c("Truth", "Estimate", "Naive"),
       col = c("purple", "orange", "black"), lwd = 2, lty = c(1, 2, 2))
lines(meanfun1(x2) ~ x2, tdata, col = "purple", lwd = 2)
lines(fitted(fit1) ~ x2, tdata, col = "orange", lwd = 2, lty = 2)
lines(fitted(lm(y1 ~ x2, tdata)) ~ x2, tdata, col = "black",
      lty = 2, lwd = 2)  # This is simplest but wrong!

plot(y2 ~ x2, data = tdata, las = 1, main = "Tobit model",
     col = as.numeric(attr(y2, "cenL")) + 3 +
           as.numeric(attr(y2, "cenU")),
     pch = as.numeric(attr(y2, "cenL")) + 1 +
           as.numeric(attr(y2, "cenU")))
legend(x = "topleft", leg = c("censored", "uncensored"),
       pch = c(2, 1), col = c("blue", "green"))
legend(-1.0, 3.5, c("Truth", "Estimate", "Naive"),
       col = c("purple", "orange", "black"), lwd = 2, lty = c(1, 2, 2))
lines(meanfun2(x2) ~ x2, tdata, col = "purple", lwd = 2)
lines(fitted(fit2) ~ x2, tdata, col = "orange", lwd = 2, lty = 2)
lines(fitted(lm(y2 ~ x2, tdata)) ~ x2, tdata, col = "black",
      lty = 2, lwd = 2)  # This is simplest but wrong!

plot(y3 ~ x2, data = tdata, las = 1,
     main = "Tobit model with nonconstant censor levels",
     col = as.numeric(attr(y3, "cenL")) + 2 +
           as.numeric(attr(y3, "cenU") * 2),
     pch = as.numeric(attr(y3, "cenL")) + 1 +
           as.numeric(attr(y3, "cenU") * 2))
legend(x = "topleft", leg = c("censored", "uncensored"),
       pch = c(2, 1), col = c("blue", "green"))
legend(-1.0, 3.5, c("Truth", "Estimate", "Naive"),
       col = c("purple", "orange", "black"), lwd = 2, lty = c(1, 2, 2))
lines(meanfun3(x2) ~ x2, tdata, col = "purple", lwd = 2)
lines(fitted(fit3) ~ x2, tdata, col = "orange", lwd = 2, lty = 2)
lines(fitted(lm(y3 ~ x2, tdata)) ~ x2, tdata, col = "black",
      lty = 2, lwd = 2)  # This is simplest but wrong!

plot(y3 ~ x2, data = tdata, las = 1,
     main = "Tobit model with nonconstant censor levels",
     col = as.numeric(attr(y3, "cenL")) + 2 +
           as.numeric(attr(y3, "cenU") * 2),
     pch = as.numeric(attr(y3, "cenL")) + 1 +
           as.numeric(attr(y3, "cenU") * 2))
legend(x = "topleft", leg = c("censored", "uncensored"),
       pch = c(2, 1), col = c("blue", "green"))
legend(-1.0, 3.5, c("Truth", "Estimate", "Naive"),
       col = c("purple", "orange", "black"), lwd = 2, lty = c(1, 2, 2))
lines(meanfun3(x2) ~ x2, data = tdata, col = "purple", lwd = 2)
lines(fitted(fit4)[, 1] ~ x2, tdata, col = "orange", lwd = 2, lty = 2)
lines(fitted(lm(y3 ~ x2, tdata)) ~ x2, data = tdata, col = "black",
      lty = 2, lwd = 2)  # This is simplest but wrong!
}
}
\keyword{models}
\keyword{regression}
