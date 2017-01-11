\name{gev}
\alias{gev}
\alias{gevff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Extreme Value Distribution Family Function }
\description{
  Maximum likelihood estimation of the
  3-parameter generalized extreme value (GEV) distribution.

}
\usage{
gev(llocation = "identitylink", lscale = "loge",
    lshape = logoff(offset = 0.5), percentiles = c(95, 99),
    ilocation = NULL, iscale = NULL, ishape = NULL, imethod = 1,
    gprobs.y = (1:9)/10, gscale.mux = exp((-5:5)/6),
    gshape = (-5:5) / 11 + 0.01,
    iprobs.y = NULL, tolshape0 = 0.001,
    type.fitted = c("percentiles", "mean"),
    zero = c("scale", "shape"))
gevff(llocation = "identitylink", lscale = "loge",
    lshape = logoff(offset = 0.5), percentiles = c(95, 99),
    ilocation = NULL, iscale = NULL, ishape = NULL, imethod = 1,
    gprobs.y = (1:9)/10, gscale.mux = exp((-5:5)/6),
    gshape = (-5:5) / 11 + 0.01,
    iprobs.y = NULL, tolshape0 = 0.001,
    type.fitted = c("percentiles", "mean"), zero = c("scale", "shape"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale, lshape}{
  Parameter link functions for \eqn{\mu}{mu}, \eqn{\sigma}{sigma} and
  \eqn{\xi}{xi} respectively.
  See \code{\link{Links}} for more choices.


  For the shape parameter,
  the default \code{\link{logoff}} link has an offset
  called \eqn{A} below; and then the linear/additive predictor is
  \eqn{\log(\xi+A)}{log(xi+A)} which means that
  \eqn{\xi > -A}{xi > -A}.
  For technical reasons (see \bold{Details}) it is a good idea
  for \eqn{A = 0.5}.


  }



%  \item{Offset}{
%  Numeric, of length 1.
%  Called \eqn{A} below.
%  Offset value if \code{lshape = "logoff"}.
%  Then the linear/additive predictor is
%  \eqn{\log(\xi+A)}{log(xi+A)} which means that
%  \eqn{\xi > -A}{xi > -A}.
%  For technical reasons (see \bold{Details}) it is a good idea for
%  \code{Offset = 0.5}.

%  }


  \item{percentiles}{
  Numeric vector of percentiles used for the fitted values.
  Values should be between 0 and 100.
  This argument is ignored if \code{type.fitted = "mean"}.



% 20140912: this is still true, but using 'type.fitted' is better.
% However, if \code{percentiles = NULL}, then the mean
% \eqn{\mu + \sigma (\Gamma(1-\xi)-1) / \xi}{mu + sigma * (gamma(1-xi)-1)/xi}
% is returned, and this is only defined if \eqn{\xi<1}{xi<1}.


  }
  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The default is to use the \code{percentiles} argument.
  If \code{"mean"} is chosen, then the mean
  \eqn{\mu + \sigma (\Gamma(1-\xi)-1) / \xi}{mu + sigma * (gamma(1-xi)-1)/xi}
  is returned as the fitted values,
  and these are only defined for \eqn{\xi<1}{xi<1}.



  }
  \item{ilocation, iscale, ishape}{
  Numeric. Initial value for the location parameter, \eqn{\sigma}{sigma} and
  \eqn{\xi}{xi}. A \code{NULL} means a value is computed internally.
  The argument \code{ishape} is more important than the other two.
  If a failure to converge occurs, or even to obtain initial values occurs,
  try assigning \code{ishape} some value
  (positive or negative; the sign can be very important).
  Also, in general, a larger value of \code{iscale} tends to be better than a
  smaller value.


% because they are initialized from the initial \eqn{\xi}{xi}.



  }
%  \item{rshape}{
%  Numeric, of length 2.
%  Range of \eqn{\xi}{xi} if \code{lshape = "extlogit"} is chosen.
%  The rationale for the default values is given below.
%  }


% \item{mean}{
% Logical. If \code{TRUE}, the mean is computed and returned
% as the fitted values. This argument overrides the
% \code{percentiles} argument.
% See \bold{Details} for more details.
% }


  \item{imethod}{
  Initialization method. Either the value 1 or 2.
  If both methods fail then try using \code{ishape}.
  See \code{\link{CommonVGAMffArguments}} for information.


% Method 1 involves choosing the best \eqn{\xi}{xi} on the grid values
% given by \code{gshape}.
% Method 2 is similar to the method of moments.


  }
  \item{gshape}{
  Numeric vector.
  The values are used for a grid search for an initial value
  for \eqn{\xi}{xi}.
  See \code{\link{CommonVGAMffArguments}} for information.


% Used only if \code{imethod} equals 1.


  }
  \item{gprobs.y, gscale.mux, iprobs.y}{
  Numeric vectors, used for the initial values.
  See \code{\link{CommonVGAMffArguments}} for information.



  }

  \item{tolshape0}{
  Passed into \code{\link{dgev}} when computing the log-likelihood.


  }
  \item{zero}{
  A specifying which
  linear/additive predictors are modelled as intercepts only.
  The values can be from the set \{1,2,3\} corresponding
  respectively to \eqn{\mu}{mu}, \eqn{\sigma}{sigma}, \eqn{\xi}{xi}.
  If \code{zero = NULL} then all linear/additive predictors are modelled as
  a linear combination of the explanatory variables.
  For many data sets having \code{zero = 3} is a good idea.
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  The GEV distribution function can be written
  \deqn{G(y) = \exp( -[ (y-\mu)/ \sigma  ]_{+}^{- 1/ \xi}) }{%
        G(y) =  exp( -[ (y- mu)/  sigma  ]_{+}^{- 1/  xi}) }
  where \eqn{\sigma > 0}{sigma > 0},
  \eqn{-\infty < \mu < \infty}{-Inf < mu < Inf},
  and \eqn{1 + \xi(y-\mu)/\sigma > 0}{1 + xi*(y-mu)/sigma > 0}.
  Here, \eqn{x_+ = \max(x,0)}{x_+ = max(x,0)}.
  The \eqn{\mu}{mu}, \eqn{\sigma}{sigma}, \eqn{\xi}{xi} are known as the
  \emph{location}, \emph{scale} and \emph{shape} parameters respectively.
  The cases
  \eqn{\xi>0}{xi>0},
  \eqn{\xi<0}{xi<0},
  \eqn{\xi = 0}{xi = 0}
  correspond to the Frechet, Weibull, and Gumbel types respectively.
  It can be noted that the Gumbel (or Type I) distribution accommodates
  many commonly-used distributions such as the normal, lognormal,
  logistic, gamma, exponential and Weibull.


  For the GEV distribution, the \eqn{k}th moment about the mean exists
  if \eqn{\xi < 1/k}{xi < 1/k}.
  Provided they exist, the mean and variance are given by
  \eqn{\mu+\sigma\{ \Gamma(1-\xi)-1\}/ \xi}{mu + sigma \{ Gamma(1-xi)-1\} / xi}
  and
  \eqn{\sigma^2 \{ \Gamma(1-2\xi) - \Gamma^2(1-\xi) \} / \xi^2}{sigma^2
  \{ Gamma(1-2 xi) - Gamma^2 (1- xi) \} / xi^2}
  respectively,
  where \eqn{\Gamma}{Gamma} is the gamma function.


  Smith (1985) established that when \eqn{\xi > -0.5}{xi > -0.5},
  the maximum likelihood estimators are completely regular.
  To have some control over the estimated \eqn{\xi}{xi} try
  using \code{lshape = logoff(offset = 0.5)}, say,
  or \code{lshape = extlogit(min = -0.5, max = 0.5)}, say.


% and when \eqn{-1 < \xi < -0.5}{-1 < xi < -0.5} they exist but are
% non-regular; and when \eqn{\xi < -1}{xi < -1} then the maximum
% likelihood estimators do not exist. In most environmental data
% sets \eqn{\xi > -1}{xi > -1} so maximum likelihood works fine.


}
\section{Warning }{
  Currently, if an estimate of \eqn{\xi}{xi} is too close to 0 then
  an error may occur for \code{gev()} with multivariate responses.
  In general, \code{gevff()} is more reliable than \code{gev()}.


  Fitting the GEV by maximum likelihood estimation can be numerically
  fraught. If \eqn{1 + \xi (y-\mu)/ \sigma \leq 0}{1 + xi*(y-mu)/sigma <=
  0} then some crude evasive action is taken but the estimation process
  can still fail. This is particularly the case if \code{\link{vgam}}
  with \code{\link{s}} is used; then smoothing is best done with
  \code{\link{vglm}} with regression splines (\code{\link[splines]{bs}}
  or \code{\link[splines]{ns}}) because \code{\link{vglm}} implements
  half-stepsizing whereas \code{\link{vgam}} doesn't (half-stepsizing
  helps handle the problem of straying outside the parameter space.)



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Yee, T. W. and Stephenson, A. G. (2007)
  Vector generalized linear and additive extreme value models.
  \emph{Extremes}, \bold{10}, 1--19.


  Tawn, J. A. (1988)
  An extreme-value theory model for dependent observations.
  \emph{Journal of Hydrology}, \bold{101}, 227--250.


  Prescott, P. and Walden, A. T. (1980)
  Maximum likelihood estimation of the parameters of the
  generalized extreme-value distribution.
  \emph{Biometrika}, \bold{67}, 723--724.


  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of nonregular cases.
  \emph{Biometrika}, \bold{72}, 67--90.


}
\author{ T. W. Yee }

\note{
  The \pkg{VGAM} family function \code{gev} can handle a multivariate
  (matrix) response, cf. multiple responses.
  If so, each row of the matrix is sorted into
  descending order and \code{NA}s are put last.
  With a vector or one-column matrix response using
  \code{gevff} will give the same result but be faster and it handles
  the \eqn{\xi = 0}{xi = 0} case.
  The function \code{gev} implements Tawn (1988) while
  \code{gevff} implements Prescott and Walden (1980).


  Function \code{egev()} has been replaced by the
  new family function \code{gevff()}. It now
  conforms to the usual \pkg{VGAM} philosophy of
  having \code{M1} linear predictors per (independent) response.
  This is the usual way multiple responses are handled.
  Hence \code{vglm(cbind(y1, y2)\ldots, gevff, \ldots)} will have
  6 linear predictors and it is possible to constrain the
  linear predictors so that the answer is similar to \code{gev()}.
  Missing values in the response of \code{gevff()} will be deleted;
  this behaviour is the same as with almost every other
  \pkg{VGAM} family function.


  The shape parameter \eqn{\xi}{xi} is difficult to estimate
  accurately unless there is a lot of data.
  Convergence is slow when \eqn{\xi}{xi} is near \eqn{-0.5}.
  Given many explanatory variables, it is often a good idea
  to make sure \code{zero = 3}.
  The range restrictions of the parameter \eqn{\xi}{xi} are not
  enforced; thus it is possible for a violation to occur.


  Successful convergence often depends on having a reasonably good initial
  value for \eqn{\xi}{xi}. If failure occurs try various values for the
  argument \code{ishape}, and if there are covariates,
  having \code{zero = 3} is advised.


}

\seealso{
  \code{\link{rgev}},
  \code{\link{gumbel}},
  \code{\link{gumbelff}},
  \code{\link{guplot}},
  \code{\link{rlplot.gevff}},
  \code{\link{gpd}},
  \code{\link{weibullR}},
  \code{\link{frechet}},
  \code{\link{extlogit}},
  \code{\link{oxtemp}},
  \code{\link{venice}},
  \code{\link{CommonVGAMffArguments}}.


  %\code{\link{gevff}},
  %\code{\link{ogev}},


}

\examples{
\dontrun{
# Multivariate example
fit1 <- vgam(cbind(r1, r2) ~ s(year, df = 3), gev(zero = 2:3),
             data = venice, trace = TRUE)
coef(fit1, matrix = TRUE)
head(fitted(fit1))
par(mfrow = c(1, 2), las = 1)
plot(fit1, se = TRUE, lcol = "blue", scol = "forestgreen",
     main = "Fitted mu(year) function (centered)", cex.main = 0.8)
with(venice, matplot(year, depvar(fit1)[, 1:2], ylab = "Sea level (cm)",
     col = 1:2, main = "Highest 2 annual sea levels", cex.main = 0.8))
with(venice, lines(year, fitted(fit1)[,1], lty = "dashed", col = "blue"))
legend("topleft", lty = "dashed", col = "blue", "Fitted 95 percentile")

# Univariate example
(fit <- vglm(maxtemp ~ 1, gevff, data = oxtemp, trace = TRUE))
head(fitted(fit))
coef(fit, matrix = TRUE)
Coef(fit)
vcov(fit)
vcov(fit, untransform = TRUE)
sqrt(diag(vcov(fit)))  # Approximate standard errors
rlplot(fit)
}
}
\keyword{models}
\keyword{regression}




%   type.fitted = c("percentiles", "mean"), giveWarning = TRUE,

% \item{gshape}{
% Numeric, of length 2.
% Range of \eqn{\xi}{xi} used for a grid search for a good initial value
% for \eqn{\xi}{xi}.
% Used only if \code{imethod} equals 1.
% }



