\name{lms.bcn}
\alias{lms.bcn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LMS Quantile Regression with a Box-Cox Transformation to Normality }
\description{
  LMS quantile regression with the Box-Cox transformation
  to normality.

}
\usage{
lms.bcn(percentiles = c(25, 50, 75), zero = c("lambda", "sigma"),
        llambda = "identitylink", lmu = "identitylink", lsigma = "loge",
        idf.mu = 4, idf.sigma = 2, ilambda = 1,
        isigma = NULL, tol0 = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{percentiles}{
  A numerical vector containing values between 0 and 100,
  which are the quantiles.
  They will be returned as `fitted values'.



% or expectiles.
% 20140624; withdrawn 'expectiles'.
%       isigma = NULL, tol0 = 0.001, expectiles = FALSE




  }
  \item{zero}{
  Can be an integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,3\}.
  The default value usually increases the chance of successful convergence.
  Setting \code{zero = NULL} means they all are
  functions of the covariates.
  For more information see \code{\link{CommonVGAMffArguments}}.


  }
  \item{llambda, lmu, lsigma}{
  Parameter link functions applied to the first, second and third
  linear/additive predictors.
  See \code{\link{Links}} for more choices,
  and \code{\link{CommonVGAMffArguments}}.


  }
  \item{idf.mu}{
  Degrees of freedom for the cubic smoothing spline fit applied to
  get an initial estimate of mu.
  See \code{\link{vsmooth.spline}}.


  }
  \item{idf.sigma}{
  Degrees of freedom for the cubic smoothing spline fit applied to
  get an initial estimate of sigma.
  See \code{\link{vsmooth.spline}}.
  This argument may be assigned \code{NULL} to get an initial value
  using some other algorithm.


  }
  \item{ilambda}{
  Initial value for lambda.
  If necessary, it is recycled to be a vector of length \eqn{n}
  where \eqn{n} is the number of (independent) observations.


  }
  \item{isigma}{
  Optional initial value for sigma.
  If necessary, it is recycled to be a vector of length \eqn{n}.
  The default value, \code{NULL}, means an initial value is computed
  in the \code{@initialize} slot of the family function.


  }
  \item{tol0}{
  Small positive number,
  the tolerance for testing if lambda is equal to zero.


  }
% \item{expectiles}{
% Experimental; please do not use.

% A single logical. If \code{TRUE} then the method is LMS-expectile
% regression; \emph{expectiles} are returned rather than quantiles.
% The default is LMS quantile regression based on the normal distribution.


% }

}
\details{

  Given a value of the covariate, this function applies
  a Box-Cox transformation to the response to best obtain
  normality. The parameters chosen to do this are estimated
  by maximum likelihood or penalized maximum likelihood.


In more detail,
the basic idea behind this method is that, for a fixed
value of \eqn{x}, a Box-Cox transformation of the response \eqn{Y}
is applied to obtain standard normality. The 3 parameters
(\eqn{\lambda}{lambda}, \eqn{\mu}{mu}, \eqn{\sigma}{sigma},
which start with the letters ``L-M-S''
respectively, hence its name) are chosen to maximize a penalized
log-likelihood (with \code{\link{vgam}}). Then the
appropriate quantiles of the standard normal distribution
are back-transformed onto the original scale to get the
desired quantiles.
The three parameters may vary as a smooth function of \eqn{x}.


The Box-Cox power transformation here of the \eqn{Y}, given \eqn{x}, is
\deqn{Z = [(Y / \mu(x))^{\lambda(x)} - 1] / ( \sigma(x) \, \lambda(x) )}{
      Z = [(Y / mu(x))^{lambda(x)} - 1] / (sigma(x) * lambda(x))}
for \eqn{\lambda(x) \neq 0}{lambda(x) != 0}.
(The singularity at \eqn{\lambda(x) = 0}{lambda(x) = 0}
is handled by a simple function involving a logarithm.)
Then \eqn{Z} is assumed to have a standard normal distribution.
The parameter \eqn{\sigma(x)}{sigma(x)} must be positive, therefore
\pkg{VGAM} chooses
\eqn{\eta(x)^T = (\lambda(x), \mu(x), \log(\sigma(x)))}{eta(x)^T =
(lambda(x), mu(x), log(sigma(x)))}
by default.
The parameter \eqn{\mu}{mu} is also positive, but while
\eqn{\log(\mu)}{log(mu)} is
available, it is not the default because \eqn{\mu}{mu} is
more directly interpretable.
Given the estimated linear/additive predictors, the
\eqn{100\alpha}{100*alpha} percentile can be estimated
by inverting the Box-Cox power transformation at the
\eqn{100\alpha}{100*alpha} percentile of the standard
normal distribution.


Of the three functions, it is often a good idea to allow
\eqn{\mu(x)}{mu(x)} to be more flexible because the functions
\eqn{\lambda(x)}{lambda(x)} and \eqn{\sigma(x)}{sigma(x)}
usually vary more smoothly with \eqn{x}. This is somewhat
reflected in the default value for the argument \code{zero},
viz. \code{zero = c(1, 3)}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{
Cole, T. J. and Green, P. J. (1992)
Smoothing Reference Centile Curves: The LMS Method and
Penalized Likelihood.
\emph{Statistics in Medicine},
\bold{11}, 1305--1319.


Green, P. J. and Silverman, B. W. (1994)
\emph{Nonparametric Regression and Generalized Linear Models: A
Roughness Penalty Approach},
London: Chapman & Hall.


Yee, T. W. (2004)
Quantile regression via vector generalized additive models.
\emph{Statistics in Medicine},
\bold{23}, 2295--2315.


%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
  The response must be positive because the Box-Cox
  transformation cannot handle negative values. The
  LMS-Yeo-Johnson-normal method can handle both positive
  and negative values.


% LMS-BCN expectile regression is a \emph{new} methodology proposed
% by myself!


  In general, the lambda and sigma functions should be more smoother
  than the mean function.
  Having \code{zero = 1}, \code{zero = 3} or \code{zero = c(1, 3)}
  is often a good idea. See the example below.



% While it is usual to regress the response against a single
% covariate, it is possible to add other explanatory variables,
% e.g., gender.
% See
% \url{http://www.stat.auckland.ac.nz/~yee}
% for further information and examples about this feature.


}

\section{Warning }{
  The computations are not simple, therefore convergence may fail.
  Set \code{trace = TRUE} to monitor convergence if it isn't set already.
  Convergence failure will occur if, e.g., the response is bimodal
  at any particular value of \eqn{x}.
  In case of convergence failure, try different starting values.
  Also, the estimate may diverge quickly near the solution,
  in which case try prematurely
  stopping the iterations by assigning \code{maxits} to be the iteration
  number corresponding to the highest likelihood value.


  One trick is to fit a simple model and use it to provide
  initial values for a more complex model; see in the
  examples below.


}
\seealso{
  \code{\link{lms.bcg}},
  \code{\link{lms.yjn}},
  \code{\link{qtplot.lmscreg}},
  \code{\link{deplot.lmscreg}},
  \code{\link{cdf.lmscreg}},
  \code{\link{alaplace1}},
  \code{\link{amlnormal}},
  \code{\link{denorm}},
  \code{\link{CommonVGAMffArguments}}.


% \code{\link{bmi.nz}},


}

\examples{
\dontrun{ require("VGAMdata")
mysubset <- subset(xs.nz, sex == "M" & ethnicity == "Maori" & study1)
mysubset <- transform(mysubset, BMI = weight / height^2)
BMIdata <- na.omit(mysubset)
BMIdata <- subset(BMIdata, BMI < 80 & age < 65,
                   select = c(age, BMI))  # Delete an outlier
summary(BMIdata)

fit <- vgam(BMI ~ s(age, df = c(4, 2)), lms.bcn(zero = 1), data = BMIdata)

par(mfrow = c(1, 2))
plot(fit, scol = "blue", se = TRUE)  # The two centered smooths

head(predict(fit))
head(fitted(fit))
head(BMIdata)
head(cdf(fit))  # Person 46 is probably overweight, given his age
100 * colMeans(depvar(fit, drop = TRUE) < fitted(fit))  # Empirical proportions

# Convergence problems? Try this trick: fit0 is a simpler model used for fit1
fit0 <- vgam(BMI ~ s(age, df = 4), lms.bcn(zero = c(1, 3)), data = BMIdata)
fit1 <- vgam(BMI ~ s(age, df = c(4, 2)), lms.bcn(zero = 1), data = BMIdata,
            etastart = predict(fit0))
}

\dontrun{
# Quantile plot
par(bty = "l", mar = c(5, 4, 4, 3) + 0.1, xpd = TRUE)
qtplot(fit, percentiles = c(5, 50, 90, 99), main = "Quantiles",
       xlim = c(15, 66), las = 1, ylab = "BMI", lwd = 2, lcol = 4)

# Density plot
ygrid <- seq(15, 43, len = 100)  # BMI ranges
par(mfrow = c(1, 1), lwd = 2)
(aa <- deplot(fit, x0 = 20, y = ygrid, xlab = "BMI", col = "black",
  main = "Density functions at Age = 20 (black), 42 (red) and 55 (blue)"))
aa <- deplot(fit, x0 = 42, y = ygrid, add = TRUE, llty = 2, col = "red")
aa <- deplot(fit, x0 = 55, y = ygrid, add = TRUE, llty = 4, col = "blue",
             Attach = TRUE)
aa@post$deplot  # Contains density function values
}
}
\keyword{models}
\keyword{regression}

% BMIdata <- subset(mysubset, select = c(age, BMI))
% BMIdata <- mysubset[, c("age", "BMI")]

