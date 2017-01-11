\name{gpd}
\alias{gpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Pareto Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter
  generalized  Pareto distribution (GPD).

}
\usage{
gpd(threshold = 0, lscale = "loge", lshape = logoff(offset = 0.5),
    percentiles = c(90, 95), iscale = NULL, ishape = NULL,
    tolshape0 = 0.001, type.fitted = c("percentiles", "mean"),
    imethod = 1, zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{threshold}{
  Numeric, values are recycled if necessary.
  The threshold value(s), called \eqn{\mu}{mu} below.


  }
  \item{lscale}{
  Parameter link function for the scale parameter \eqn{\sigma}{sigma}.
  See \code{\link{Links}} for more choices.


  }
  \item{lshape}{
  Parameter link function for the shape parameter \eqn{\xi}{xi}.
  See \code{\link{Links}} for more choices.
  The default constrains the parameter to be greater than \eqn{-0.5}
  because if \eqn{\xi \leq -0.5}{xi <= -0.5} then Fisher
  scoring does not work.
  See the Details section below for more information.


  For the shape parameter,
  the default \code{\link{logoff}} link has an offset
  called \eqn{A} below; and then the second linear/additive predictor is
  \eqn{\log(\xi+A)}{log(xi+A)} which means that
  \eqn{\xi > -A}{xi > -A}.
  The working weight matrices are positive definite if \eqn{A = 0.5}.


  }

% \item{Offset}{
% Numeric, of length 1.
% Called \eqn{A} below.
% Offset value if \code{lshape = "logoff"}.
% Then the second linear/additive predictor is
% \eqn{\log(\xi+A)}{log(xi+A)} which means that
% \eqn{\xi > -A}{xi > -A}.
% The working weight matrices are positive definite if \code{Offset = 0.5}.

% }
  \item{percentiles}{
  Numeric vector of percentiles used
  for the fitted values. Values should be between 0 and 100.
  See the example below for illustration.
  This argument is ignored if \code{type.fitted = "mean"}.



% However, if \code{percentiles = NULL} then the mean
% \eqn{\mu + \sigma / (1-\xi)}{mu + sigma / (1-xi)} is returned;
% this is only defined if \eqn{\xi<1}{xi<1}.


  }
  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The default is to use the \code{percentiles} argument.
  If \code{"mean"} is chosen, then the mean
  \eqn{\mu + \sigma / (1-\xi)}{mu + sigma / (1-xi)}
  is returned as the fitted values,
  and these are only defined for \eqn{\xi<1}{xi<1}.



  }

  \item{iscale, ishape}{
  Numeric. Optional initial values for \eqn{\sigma}{sigma}
  and \eqn{\xi}{xi}.
  The default is to use \code{imethod} and compute a value internally for
  each parameter.
  Values of \code{ishape} should be between \eqn{-0.5} and \eqn{1}.
  Values of \code{iscale} should be positive.


  }
% \item{rshape}{
% Numeric, of length 2.
% Range of \eqn{\xi}{xi} if \code{lshape = "extlogit"} is chosen.
% The default values ensures the algorithm works (\eqn{\xi > -0.5}{xi > -0.5})
% and the variance exists (\eqn{\xi < 0.5}{xi < 0.5}).


% }
  \item{tolshape0}{
  Passed into \code{\link{dgpd}} when computing the log-likelihood.


  }

% \item{tolshape0}{
% Positive numeric.
% Threshold/tolerance value for resting whether \eqn{\xi}{xi} is zero.
% If the absolute value of the estimate of \eqn{\xi}{xi} is less than
% this value then it will be assumed zero and exponential distribution
% derivatives etc. will be used.


% }
  \item{imethod}{
   Method of initialization, either 1 or 2. The first is the method of
   moments, and the second is a variant of this.  If neither work, try
   assigning values to arguments \code{ishape} and/or \code{iscale}.


  }
  \item{zero}{
  Can be an integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  For one response, the value should be from the set \{1,2\} corresponding
  respectively to \eqn{\sigma}{sigma} and \eqn{\xi}{xi}.
  It is often a good idea for the \eqn{\sigma}{sigma} parameter only
  to be modelled through
  a linear combination of the explanatory variables because the
  shape parameter is probably best left as an intercept only:
  \code{zero = 2}.
  Setting \code{zero = NULL} means both parameters are modelled with
  explanatory variables.
  See \code{\link{CommonVGAMffArguments}} for more details.




  }
}
\details{
  The distribution function of the GPD can be written
    \deqn{G(y) = 1 - [1 + \xi (y-\mu) / \sigma  ]_{+}^{- 1/ \xi} }{%
          G(y) = 1 - [1 +  xi (y-mu)/  sigma  ]_{+}^{- 1/  xi}  }
  where
  \eqn{\mu}{mu} is the location parameter
  (known, with value \code{threshold}),
  \eqn{\sigma > 0}{sigma > 0} is the scale parameter,
    \eqn{\xi}{xi} is the shape parameter, and
  \eqn{h_+ = \max(h,0)}{h_+ = max(h,0)}.
  The function \eqn{1-G} is known as the \emph{survivor function}.
  The limit \eqn{\xi \rightarrow 0}{xi --> 0}
  gives the \emph{shifted exponential} as a special case:
    \deqn{G(y) = 1 - \exp[-(y-\mu)/ \sigma]. }{%
          G(y) = 1 -  exp[-(y-mu)/  sigma]. }
  The support is \eqn{y>\mu}{y>mu} for \eqn{\xi>0}{xi>0},
  and
  \eqn{\mu < y <\mu-\sigma / \xi}{mu < y <mu-sigma / xi} for \eqn{\xi<0}{xi<0}.


  Smith (1985) showed that if \eqn{\xi <= -0.5}{xi <= -0.5} then
  this is known as the nonregular case and problems/difficulties
  can arise both theoretically and numerically. For the (regular)
  case \eqn{\xi > -0.5}{xi > -0.5} the classical asymptotic
  theory of maximum likelihood estimators is applicable; this is
  the default.


  Although for \eqn{\xi < -0.5}{xi < -0.5} the usual asymptotic properties
  do not apply, the maximum likelihood estimator generally exists and
  is superefficient for \eqn{-1 < \xi < -0.5}{-1 < xi < -0.5}, so it is
  ``better'' than normal.
  When \eqn{\xi < -1}{xi < -1} the maximum
  likelihood estimator generally does not exist as it effectively becomes
  a two parameter problem.


  The mean of \eqn{Y} does not exist unless \eqn{\xi < 1}{xi < 1}, and
  the variance does not exist unless \eqn{\xi < 0.5}{xi < 0.5}.  So if
  you want to fit a model with finite variance use \code{lshape = "extlogit"}.


}
\note{
  The response in the formula of \code{\link{vglm}}
  and \code{\link{vgam}} is \eqn{y}.
  Internally, \eqn{y-\mu}{y-mu} is computed.
  This \pkg{VGAM} family function can handle a multiple
  responses, which is inputted as a matrix.
  The response stored on the object is the original uncentred data.



  With functions \code{\link{rgpd}}, \code{\link{dgpd}}, etc., the
  argument \code{location} matches with the argument \code{threshold}
  here.


}
\section{Warning}{
  Fitting the GPD by maximum likelihood estimation can be numerically
  fraught. If \eqn{1 + \xi (y-\mu)/ \sigma \leq 0}{1 + xi*(y-mu)/sigma <=
  0} then some crude evasive action is taken but the estimation process
  can still fail. This is particularly the case if \code{\link{vgam}}
  with \code{\link{s}} is used. Then smoothing is best done with
  \code{\link{vglm}} with regression splines (\code{\link[splines]{bs}}
  or \code{\link[splines]{ns}}) because \code{\link{vglm}} implements
  half-stepsizing whereas \code{\link{vgam}} doesn't. Half-stepsizing
  helps handle the problem of straying outside the parameter space.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.
  However, for this \pkg{VGAM} family function, \code{\link{vglm}}
  is probably preferred over \code{\link{vgam}} when there is smoothing.

}
\references{
  Yee, T. W. and Stephenson, A. G. (2007)
  Vector generalized linear and additive extreme value models.
  \emph{Extremes}, \bold{10}, 1--19.


  Coles, S. (2001)
  \emph{An Introduction to Statistical Modeling of Extreme Values}.
  London: Springer-Verlag.


  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of nonregular cases.
  \emph{Biometrika}, \bold{72}, 67--90.


}
\author{ T. W. Yee }

\seealso{
  \code{\link{rgpd}},
  \code{\link{meplot}},
  \code{\link{gev}},
  \code{\link{paretoff}},
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{s}}.

}

\examples{
# Simulated data from an exponential distribution (xi = 0)
Threshold <- 0.5
gdata <- data.frame(y1 = Threshold + rexp(n = 3000, rate = 2))
fit <- vglm(y1 ~ 1, gpd(threshold = Threshold), data = gdata, trace = TRUE)
head(fitted(fit))
summary(depvar(fit))  # The original uncentred data
coef(fit, matrix = TRUE)  # xi should be close to 0
Coef(fit)
summary(fit)

head(fit@extra$threshold)  # Note the threshold is stored here

# Check the 90 percentile
ii <- depvar(fit) < fitted(fit)[1, "90\%"]
100 * table(ii) / sum(table(ii))  # Should be 90%

# Check the 95 percentile
ii <- depvar(fit) < fitted(fit)[1, "95\%"]
100 * table(ii) / sum(table(ii))  # Should be 95%

\dontrun{ plot(depvar(fit), col = "blue", las = 1,
               main = "Fitted 90\% and 95\% quantiles")
matlines(1:length(depvar(fit)), fitted(fit), lty = 2:3, lwd = 2) }


# Another example
gdata <- data.frame(x2 = runif(nn <- 2000))
Threshold <- 0; xi <- exp(-0.8) - 0.5
gdata <- transform(gdata, y2 = rgpd(nn, scale = exp(1 + 0.1*x2), shape = xi))
fit <- vglm(y2 ~ x2, gpd(Threshold), data = gdata, trace = TRUE)
coef(fit, matrix = TRUE)


\dontrun{ # Nonparametric fits
# Not so recommended:
fit1 <- vgam(y2 ~ s(x2), gpd(Threshold), data = gdata, trace = TRUE)
par(mfrow = c(2, 1))
plot(fit1, se = TRUE, scol = "blue")
# More recommended:
fit2 <- vglm(y2 ~ sm.bs(x2), gpd(Threshold), data = gdata, trace = TRUE)
plot(as(fit2, "vgam"), se = TRUE, scol = "blue") }
}
\keyword{models}
\keyword{regression}

%
% # gdata <- transform(gdata, yy = y2 + rnorm(nn, sd = 0.1))
%
%   giveWarning = TRUE, imethod = 1, zero = "shape"
